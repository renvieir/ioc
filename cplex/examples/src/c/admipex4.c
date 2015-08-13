/* --------------------------------------------------------------------------
 * File: admipex4.c
 * Version 12.6.1
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
 * Copyright IBM Corporation 1997, 2014. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 * --------------------------------------------------------------------------
 */

/* Examples admipex4.c and admipex5.c both solve the MIPLIB
   3.0 model noswot.mps by adding user cuts.  admipex4.c adds
   these cuts to the cut table before the branch-and-cut
   process begins, while admipex5.c adds them through the
   user cut callback during the branch-and-cut process. 

   When this example is run the program reads a problem from a
   file named "noswot.mps" either from the directory ../../data
   if no argument is passed to the executable or from the directory
   that is specified as the first (and only) argument to the executable. */


/* Bring in the CPLEX function declarations and the C library 
   header file stdio.h with the following single include */

#include <ilcplex/cplex.h>

/* Bring in the declarations for the string and character functions, 
   malloc, and fabs */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Declarations for the functions in this program */

static int
   addusercuts   (CPXENVptr env, CPXLPptr lp);

static void
   free_and_null (char **ptr);

int
main (int argc, char *argv[])
{
   int status = 0;

   /* Declare and allocate space for the variables and arrays where
      we will store the optimization results, including the status, 
      objective value, and variable values */
   
   int    solstat;
   double objval;
   double *x = NULL;
   
   CPXENVptr env = NULL;
   CPXLPptr  lp = NULL;

   int j;
   int cur_numcols;

   const char * datadir = argc <= 1 ? "../../../examples/data" : argv[1];
   char *noswot = NULL;

   noswot = (char *) malloc (strlen (datadir) + 1 + strlen("noswot.mps") + 1);
   sprintf (noswot, "%s/noswot.mps", datadir);

   /* Initialize the CPLEX environment */

   env = CPXopenCPLEX (&status);

   /* If an error occurs, the status value indicates the reason for
      failure.  A call to CPXgeterrorstring will produce the text of
      the error message.  Note that CPXopenCPLEX produces no
      output, so the only way to see the cause of the error is to use
      CPXgeterrorstring.  For other CPLEX routines, the errors will
      be seen if the CPXPARAM_ScreenOutput parameter is set to CPX_ON */

   if ( env == NULL ) {
      char errmsg[CPXMESSAGEBUFSIZE];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      goto TERMINATE;
   }

   /* Turn on output to the screen */

   status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON);
   if ( status != 0 ) {
      fprintf (stderr, 
               "Failure to turn on screen indicator, error %d.\n",
               status);
      goto TERMINATE;
   }
   CPXsetintparam (env, CPXPARAM_MIP_Interval, 1000);

   /* Create the problem, using the filename as the problem name */

   lp = CPXcreateprob (env, &status, "noswot");

   /* A returned pointer of NULL may mean that not enough memory
      was available or there was some other problem.  In the case of
      failure, an error message will have been written to the error
      channel from inside CPLEX.  In this example, the setting of
      the parameter CPXPARAM_ScreenOutput causes the error message to
      appear on stdout.  Note that most CPLEX routines return
      an error code to indicate the reason for failure */

   if ( lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      goto TERMINATE;
   }

   /* Now read the file, and copy the data into the created lp */

   status = CPXreadcopyprob (env, lp, noswot, NULL);
   if ( status ) {
      fprintf (stderr,
               "Failed to read and copy the problem data.\n");
      goto TERMINATE;
   }

   /* Set parameters */

   /* Assure linear mappings between the presolved and original
      models */

   status = CPXsetintparam (env, CPXPARAM_Preprocessing_Linear, 0);
   if ( status )  goto TERMINATE;


   /* Create user cuts for noswot problem */

   status = addusercuts (env, lp); 
   if ( status )  goto TERMINATE;

   /* Optimize the problem and obtain solution */

   status = CPXmipopt (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize MIP.\n");
      goto TERMINATE;
   }

   solstat = CPXgetstat (env, lp);
   printf ("Solution status %d.\n", solstat);

   status = CPXgetobjval (env, lp, &objval);
   if ( status ) {
      fprintf (stderr,"Failed to obtain objective value.\n");
      goto TERMINATE;
   }

   printf ("Objective value %.10g\n", objval);

   cur_numcols = CPXgetnumcols (env, lp);

   /* Allocate space for solution */

   x = (double *) malloc (cur_numcols * sizeof (double));
   if ( x == NULL ) {
      fprintf (stderr, "No memory for solution values.\n");
      goto TERMINATE;
   }

   status = CPXgetx (env, lp, x, 0, cur_numcols-1);
   if ( status ) {
      fprintf (stderr, "Failed to obtain solution.\n");
      goto TERMINATE;
   }

   /* Write out the solution */

   for (j = 0; j < cur_numcols; j++) {
      if ( fabs (x[j]) > 1e-10 ) {
         printf ("Column %d:  Value = %17.10g\n", j, x[j]);
      }
   }


TERMINATE:

   /* Free the filename */

   free_and_null ((char **) &noswot);

   /* Free the solution vector */

   free_and_null ((char **) &x);

   /* Free the problem as allocated by CPXcreateprob and
      CPXreadcopyprob, if necessary */

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n",
                  status);
      }
   }

   /* Free the CPLEX environment, if necessary */

   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      /* Note that CPXcloseCPLEX produces no output, so the only 
         way to see the cause of the error is to use
         CPXgeterrorstring.  For other CPLEX routines, the errors 
         will be seen if the CPXPARAM_ScreenOutput parameter is set to 
         CPX_ON */

      if ( status ) {
         char errmsg[CPXMESSAGEBUFSIZE];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
   }
     
   return (status);

} /* END main */


/* This simple routine frees up the pointer *ptr, and sets *ptr 
   to NULL */

static void
free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */ 


/* Valid cuts for noswot.mps: 
   cut1: X21 - X22 <= 0
   cut2: X22 - X23 <= 0
   cut3: X23 - X24 <= 0
   cut4: 2.08 X11 + 2.98 X21 + 3.47 X31 + 2.24 X41 + 2.08 X51 
       + 0.25 W11 + 0.25 W21 + 0.25 W31 + 0.25 W41 + 0.25 W51
         <= 20.25
   cut5: 2.08 X12 + 2.98 X22 + 3.47 X32 + 2.24 X42 + 2.08 X52 
       + 0.25 W12 + 0.25 W22 + 0.25 W32 + 0.25 W42 + 0.25 W52
         <= 20.25
   cut6: 2.08 X13 + 2.98 X23 + 3.4722 X33 + 2.24 X43 + 2.08 X53
       + 0.25 W13 + 0.25 W23 + 0.25 W33 + 0.25 W43 + 0.25 W53
         <= 20.25
   cut7: 2.08 X14 + 2.98 X24 + 3.47 X34 + 2.24 X44 + 2.08 X54 
       + 0.25 W14 + 0.25 W24 + 0.25 W34 + 0.25 W44 + 0.25 W54
         <= 20.25
   cut8: 2.08 X15 + 2.98 X25 + 3.47 X35 + 2.24 X45 + 2.08 X55 
       + 0.25 W15 + 0.25 W25 + 0.25 W35 + 0.25 W45 + 0.25 W55
         <= 16.25
*/

static int
addusercuts (CPXENVptr env,
             CPXLPptr  lp)
{
   int status = 0;

   int cutbeg[] = {0, 2, 4, 6, 16, 26, 36, 46, 56};

   double cutval[] = 
      {1, -1, 
       1, -1, 
       1, -1, 
       2.08, 2.98, 3.47, 2.24, 2.08, 0.25, 0.25, 0.25, 0.25, 0.25,
       2.08, 2.98, 3.47, 2.24, 2.08, 0.25, 0.25, 0.25, 0.25, 0.25,
       2.08, 2.98, 3.47, 2.24, 2.08, 0.25, 0.25, 0.25, 0.25, 0.25,
       2.08, 2.98, 3.47, 2.24, 2.08, 0.25, 0.25, 0.25, 0.25, 0.25,
       2.08, 2.98, 3.47, 2.24, 2.08, 0.25, 0.25, 0.25, 0.25, 0.25};

   char *varname[] = 
      {"X21", "X22", 
       "X22", "X23", 
       "X23", "X24",
       "X11", "X21", "X31", "X41", "X51",
       "W11", "W21", "W31", "W41", "W51",
       "X12", "X22", "X32", "X42", "X52",
       "W12", "W22", "W32", "W42", "W52",
       "X13", "X23", "X33", "X43", "X53",
       "W13", "W23", "W33", "W43", "W53",
       "X14", "X24", "X34", "X44", "X54",
       "W14", "W24", "W34", "W44", "W54",
       "X15", "X25", "X35", "X45", "X55",
       "W15", "W25", "W35", "W45", "W55"};
   double cutrhs[] = {0, 0, 0, 20.25, 20.25, 20.25, 20.25, 16.25};

   char cutsens[] = "LLLLLLLL";
   int  cutind[56];

   int i, varind;
   int nz = 56;
   int cuts = 8;

   for (i = 0; i < nz; i++) {
      status = CPXgetcolindex (env, lp, varname[i], &varind);
      if ( status )  {
         fprintf (stderr,
                  "Failed to get index from variable name.\n");
         goto TERMINATE;
      }
      cutind[i] = varind;
   }

   status = CPXaddusercuts (env, lp, cuts, nz, cutrhs, 
                            cutsens, cutbeg, cutind, cutval, NULL);

   /* Alternative approach:  Replace CPXaddusercuts() by
      the following CPXaddrows() call, and delete the call to 
      CPXsetintparam (env, CPXPARAM_Preprocessing_Linear, 0) in main().

      Note that in this example, we're adding a small number
      of cuts, so it is reasonable (preferable, in fact) to
      to add them directly to the model.  If you have
      a large number of cuts, however, they can slow down
      the linear programming relaxations significantly, so
      placing them in the user cut table rather than in the
      model itself will usually be preferable.

   status = CPXaddrows (env, lp, 0, cuts, nz, cutrhs, 
                             cutsens, cutbeg, cutind, cutval,
                             NULL, NULL);
   */ 

TERMINATE:
 
   return (status);

} /* END addusercuts */

