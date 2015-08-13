/* --------------------------------------------------------------------------
 * File: admipex5.c
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

/* Example admipex5.c solves the MIPLIB 3.0 model noswot.mps by adding
   user cuts via a user cut callback during the branch-and-cut process.

   Then it modifies the problem by specifying a lazy constraint generator that
   tests against one single cut that the original optimal solution does not
   satisfy.

   Finally the modified problem is solved again, this time without user cuts
   (but the lazy constraint generator is still active).

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

/* The following macro defines the smallest improvement 
   on the value of the objective function that is required
   when adding user cuts from within a callback.
   If the improvement on the value of the ojective function
   is not large enough, the callback will abort the cut loop. */

#define EPSOBJ 0.1

/* The following structure will hold the information we need to 
   pass to the cut callback function */

struct cutinfo {
   CPXLPptr lp;
   int      numcols;
   int      num;
   double   *x;
   int      *beg;
   int      *ind; 
   double   *val;
   double   *rhs;
   int      nodeid;
   double   nodeobjval;
   int      objsen;
};
typedef struct cutinfo CUTINFO, *CUTINFOptr;

/* Declarations for functions in this program */

static int
   myoptimize         (CPXENVptr env, CPXLPptr lp, double true_objval);

static int CPXPUBLIC
   mycutcallback      (CPXCENVptr env, void *cbdata, int wherefrom,
                       void *cbhandle, int *useraction_p);

static int
   makeusercuts       (CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo);

static int
   makelazyconstraint (CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo);

static void
   initnodeobjvalinfo (CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo);

static void
   free_and_null      (char **ptr);

int
main (int argc, char *argv[])
{
   int status = 0;

   /* Declare and allocate space for the variables and arrays where
      we will store the optimization results, including the status, 
      objective value, and variable values */
   
   CPXENVptr env = NULL;
   CPXLPptr  lp = NULL;

   CUTINFO usercutinfo;
   CUTINFO lazyconinfo;

   const char * datadir = argc <= 1 ? "../../../examples/data" : argv[1];
   char *noswot = NULL;

   usercutinfo.x   = NULL;
   usercutinfo.beg = NULL;
   usercutinfo.ind = NULL; 
   usercutinfo.val = NULL;
   usercutinfo.rhs = NULL;
  
   lazyconinfo.x   = NULL;
   lazyconinfo.beg = NULL;
   lazyconinfo.ind = NULL; 
   lazyconinfo.val = NULL;
   lazyconinfo.rhs = NULL;

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

   /* The problem will be solved several times, so turn off advanced start */
   
   status = CPXsetintparam (env, CPXPARAM_Advance, CPX_OFF);
   if ( status )  goto TERMINATE;

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


   /* Turn on traditional search for use with control callbacks */

   status = CPXsetintparam (env, CPXPARAM_MIP_Strategy_Search,
                            CPX_MIPSEARCH_TRADITIONAL);
   if ( status )  goto TERMINATE;

   /* Let MIP callbacks work on the original model */

   status = CPXsetintparam (env, CPXPARAM_MIP_Strategy_CallbackReducedLP,
                            CPX_OFF);
   if ( status )  goto TERMINATE;

   /* Create user cuts for noswot problem */

   status = makeusercuts (env, lp, &usercutinfo);
   if ( status )  goto TERMINATE;

   /* Init information on the node objval for the user cut callback */

   initnodeobjvalinfo (env, lp, &usercutinfo);

   /* Set up to use MIP usercut callback */

   status = CPXsetusercutcallbackfunc (env, mycutcallback, &usercutinfo);
   if ( status )  goto TERMINATE;

   status = myoptimize (env, lp, -41.0);
   if ( status )  goto TERMINATE;

   /*=======================================================================*/

   /* Create a lazy constraint to alter the optimum */

   status = makelazyconstraint (env, lp, &lazyconinfo);
   if ( status )  goto TERMINATE;

   /* Init information on the node objval for the user cut callback.
      No need to initialize the information on the node objval,
      for the lazy constraint callback, because those information are
      used only in the user cut callback. */

   initnodeobjvalinfo (env, lp, &usercutinfo);

   /* Set up to use MIP lazyconstraint callback. The callback funtion
    * registered is the same, but the data will be different. */

   status = CPXsetlazyconstraintcallbackfunc (env, mycutcallback, &lazyconinfo);
   if ( status )  goto TERMINATE;

   status = myoptimize (env, lp, -39.0);
   if ( status )  goto TERMINATE;

   /*=======================================================================*/

   /* Now solve the problem without usercut callback */
   
   status = CPXsetusercutcallbackfunc (env, NULL, NULL);
   if ( status )  goto TERMINATE;

   status = myoptimize (env, lp, -39.0);
   if ( status )  goto TERMINATE;

   /*=======================================================================*/


TERMINATE:

   /* Free the filename */

   free_and_null ((char **) &noswot);

   /* Free the allocated vectors */

   free_and_null ((char **) &usercutinfo.x);
   free_and_null ((char **) &usercutinfo.beg);
   free_and_null ((char **) &usercutinfo.ind);
   free_and_null ((char **) &usercutinfo.val);
   free_and_null ((char **) &usercutinfo.rhs);
   free_and_null ((char **) &lazyconinfo.x);
   free_and_null ((char **) &lazyconinfo.beg);
   free_and_null ((char **) &lazyconinfo.ind);
   free_and_null ((char **) &lazyconinfo.val);
   free_and_null ((char **) &lazyconinfo.rhs);

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

int
myoptimize (CPXENVptr env, CPXLPptr lp, double true_objval)
{
   int status = 0;
   int solstat = 0;

   int j;
   double objval;
   double *x = NULL;
   int cur_numcols = CPXgetnumcols (env, lp);

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
         char *colname[1];
         char namestore[6];
         int surplus = 0;
         status = CPXgetcolname (env, lp, colname, namestore, 6,
                                 &surplus, j, j);
         if ( status ) {
            namestore[0] = 0;
            colname[0] = namestore;
         }
         printf ( "Column %3d (%5s):  Value = %17.10g\n", j, colname[0], x[j]);
      }
   }


TERMINATE:

   free_and_null ((char **) &x);

   return (status);
}

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


static int CPXPUBLIC 
mycutcallback (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               int        *useraction_p)
{
   int status = 0;

   CUTINFOptr cutinfo = (CUTINFOptr) cbhandle;

   int      numcols  = cutinfo->numcols;
   int      numcuts  = cutinfo->num;
   double   *x       = cutinfo->x;
   int      *beg     = cutinfo->beg;
   int      *ind     = cutinfo->ind;
   double   *val     = cutinfo->val;
   double   *rhs     = cutinfo->rhs;
   int      *cutind  = NULL;
   double   *cutval  = NULL;
   double   cutvio;
   int      addedcuts = 0;
   int      i, j, k, cutnz;

   *useraction_p = CPX_CALLBACK_DEFAULT; 

   /* If we are called as a user cut callback, decide
      first if we want to add cuts or abort the cut loop.
      When adding user cuts with purgeable flag set to 
      CPX_USECUT_PURGE or CPX_USECUT_FILTER, it is important 
      to avoid the possibility of an infinite cut loop, where 
      the same cuts are added to the LP and then immediately 
      purged at every cut pass. Such a situation can be avoided,
      for instance, by applying a tailing off criterion and aborting 
      the cut loop where no progress in the objval is observed.
      Note, however, that the same approach must not be applied
      with lazy constraints. In this case, if lazy constraints are
      added with purgeable flag set to CPX_USECUT_PURGE, adding
      the same lazy constraint more than once could be required
      to ensure the correctness of the final result. */

   if ( wherefrom == CPX_CALLBACK_MIP_CUT_LOOP ||
        wherefrom == CPX_CALLBACK_MIP_CUT_LAST   ) {
      int    oldnodeid     = cutinfo->nodeid;
      double oldnodeobjval = cutinfo->nodeobjval;

      /* Retrieve nodeid and node objval of the current node */

      status = CPXgetcallbacknodeinfo (env, cbdata, wherefrom, 0,
                                       CPX_CALLBACK_INFO_NODE_SEQNUM,
                                       &cutinfo->nodeid);
      if ( status ) {
         fprintf(stderr, "Failed to get node id.\n");
         goto TERMINATE;
      }

      status = CPXgetcallbacknodeinfo (env, cbdata, wherefrom, 0,
                                       CPX_CALLBACK_INFO_NODE_OBJVAL,
                                       &cutinfo->nodeobjval);
      if ( status ) {
         fprintf(stderr, "Failed to get node objval.\n");
         goto TERMINATE;
      }

      /* Abort the cut loop if we are stuck at the same node
         as before and there is no progress in the node objval */

      if ( oldnodeid == cutinfo->nodeid ) {
         double objchg = (cutinfo->nodeobjval - oldnodeobjval);
         /* Multiply objchg by objsen to normalize 
            the change in the objective function to 
            the case of a minimization problem */
         objchg *= cutinfo->objsen;
         if ( objchg <= EPSOBJ ) {
            *useraction_p = CPX_CALLBACK_ABORT_CUT_LOOP;
            goto TERMINATE;
         }
      }
   }

   /* If we reached this point, we are 
      .. in a lazyconstraint callback, or 
      .. in a user cut callback, and cuts seem to help 
         improving the node objval. 
      In both cases, we retrieve the x solution and 
      look for violated cuts. */

   status = CPXgetcallbacknodex (env, cbdata, wherefrom, x,
                                 0, numcols-1); 
   if ( status ) {
      fprintf(stderr, "Failed to get node solution.\n");
      goto TERMINATE;
   }

   for (i = 0; i < numcuts; i++) {
      cutvio = -rhs[i];
      k = beg[i];
      cutnz = beg[i+1] - k;
      cutind = ind + k;
      cutval = val + k;
      for (j = 0; j < cutnz; j++) {
         cutvio += x[cutind[j]] * cutval[j];
      }

      /* Use a cut violation tolerance of 0.01 */

      if ( cutvio > 0.01 ) { 
         status = CPXcutcallbackadd (env, cbdata, wherefrom,
                                     cutnz, rhs[i], 'L',
                                     cutind, cutval, CPX_USECUT_PURGE);
         if ( status ) {
            fprintf (stderr, "Failed to add cut.\n");
            goto TERMINATE;
         }
         addedcuts++;
      }
   }

   /* Tell CPLEX that cuts have been created */ 
   if ( addedcuts > 0 ) {
      *useraction_p = CPX_CALLBACK_SET; 
   }

TERMINATE:

   return (status);

} /* END mycutcallback */


/* Valid cuts for noswot 
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
makeusercuts (CPXENVptr  env,
              CPXLPptr   lp,
              CUTINFOptr usercutinfo)
{
   int status = 0;

   int beg[] = {0, 2, 4, 6, 16, 26, 36, 46, 56};

   double val[] = 
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

   double rhs[] = {0, 0, 0, 20.25, 20.25, 20.25, 20.25, 16.25};

   int    *cutbeg = NULL;
   int    *cutind = NULL;
   double *cutval = NULL;
   double *cutrhs = NULL; 

   int i, varind;
   int nz   = 56;
   int cuts = 8;

   int cur_numcols = CPXgetnumcols (env, lp);

   usercutinfo->lp = lp;
   usercutinfo->numcols = cur_numcols;

   usercutinfo->x = (double *) malloc (cur_numcols * sizeof (double));
   if ( usercutinfo->x == NULL ) {
      fprintf (stderr, "No memory for solution values.\n");
      goto TERMINATE;
   }

   cutbeg = (int *)    malloc ((cuts+1) * sizeof (int));
   cutind = (int *)    malloc (nz * sizeof (int));
   cutval = (double *) malloc (nz * sizeof (double));
   cutrhs = (double *) malloc (cuts * sizeof (double));

   if ( cutbeg == NULL ||
        cutind == NULL ||
        cutval == NULL ||
        cutrhs == NULL   ) {
      fprintf (stderr, "No memory.\n");
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   } 
      
   for (i = 0; i < nz; i++) {
      status = CPXgetcolindex (env, lp, varname[i], &varind);
      if ( status )  {
         fprintf (stderr,
                  "Failed to get index from variable name.\n");
         goto TERMINATE;
      }
      cutind[i] = varind;
      cutval[i] = val[i];
   }

   for (i = 0; i < cuts; i++) {
      cutbeg[i] = beg[i];
      cutrhs[i] = rhs[i];
   }
   cutbeg[cuts] = beg[cuts];

   usercutinfo->num      = cuts;
   usercutinfo->beg      = cutbeg;
   usercutinfo->ind      = cutind;
   usercutinfo->val      = cutval;
   usercutinfo->rhs      = cutrhs;

TERMINATE:

   if ( status ) {
      free_and_null ((char **) &cutbeg);
      free_and_null ((char **) &cutind);
      free_and_null ((char **) &cutval);
      free_and_null ((char **) &cutrhs);
   }
 
   return (status);

} /* END makeusercuts */

/* A constraint that cuts off the optimal solution of noswot:
   W11 + W12 + W13 + W14 + W15 <= 3
 */

static int
makelazyconstraint (CPXENVptr  env,
                    CPXLPptr   lp,
                    CUTINFOptr lazyconinfo)
{
   int status = 0;

   int beg[] = {0, 5};
   double val[] = {1, 1, 1, 1, 1};
   char *varname[] = { "W11", "W12", "W13", "W14", "W15" };
   double rhs[] = {3};

   int    *cutbeg = NULL;
   int    *cutind = NULL;
   double *cutval = NULL;
   double *cutrhs = NULL; 

   int i, varind;
   int nz   = 5;
   int cuts = 1;

   int cur_numcols = CPXgetnumcols (env, lp);

   lazyconinfo->lp = lp;
   lazyconinfo->numcols = cur_numcols;

   lazyconinfo->x = (double *) malloc (cur_numcols * sizeof (double));
   if ( lazyconinfo->x == NULL ) {
      fprintf (stderr, "No memory for solution values.\n");
      goto TERMINATE;
   }

   cutbeg = (int *)    malloc ((cuts+1) * sizeof (int));
   cutind = (int *)    malloc (nz * sizeof (int));
   cutval = (double *) malloc (nz * sizeof (double));
   cutrhs = (double *) malloc (cuts * sizeof (double));

   if ( cutbeg == NULL ||
        cutind == NULL ||
        cutval == NULL ||
        cutrhs == NULL   ) {
      fprintf (stderr, "No memory.\n");
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   } 
      
   for (i = 0; i < nz; i++) {
      status = CPXgetcolindex (env, lp, varname[i], &varind);
      if ( status )  {
         fprintf (stderr,
                  "Failed to get index from variable name.\n");
         goto TERMINATE;
      }
      cutind[i] = varind;
      cutval[i] = val[i];
   }

   for (i = 0; i < cuts; i++) {
      cutbeg[i] = beg[i];
      cutrhs[i] = rhs[i];
   }
   cutbeg[cuts] = beg[cuts];

   lazyconinfo->num      = cuts;
   lazyconinfo->beg      = cutbeg;
   lazyconinfo->ind      = cutind;
   lazyconinfo->val      = cutval;
   lazyconinfo->rhs      = cutrhs;

TERMINATE:

   if ( status ) {
      free_and_null ((char **) &cutbeg);
      free_and_null ((char **) &cutind);
      free_and_null ((char **) &cutval);
      free_and_null ((char **) &cutrhs);
   }
 
   return (status);

} /* END makelazyconstraint */

/* Init information on the node objval for the user cut callback */

static void 
initnodeobjvalinfo (CPXENVptr env, CPXLPptr lp, CUTINFOptr cutinfo)
{
   cutinfo->nodeid = -1;
   cutinfo->nodeobjval = 0.0;
   cutinfo->objsen = CPXgetobjsen (env, lp);
   if ( cutinfo->objsen == CPX_MIN )
      cutinfo->objsen = 1;
   else
      cutinfo->objsen = -1;

} /* END initnodeobjvalinfo */

   
