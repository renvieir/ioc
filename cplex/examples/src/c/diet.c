/* --------------------------------------------------------------------------
 * File: diet.c
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

/* diet.c - reading data for a dietary problem, building the model
   and solving it.  methods for creating the model:

      diet -r <datafile>   generates the problem by adding rows
      diet -c <datafile>   generates the problem by adding columns
 */

/* Bring in the CPLEX function declarations and the C library 
   header file stdio.h with the following single include. */

#include <ilcplex/cplex.h>

/* Bring in the declarations for the string functions */

#include <stdlib.h>
#include <string.h>

/* Include declaration for functions at end of program */


static int
   readarray         (FILE *in, int *num_p, double **data_p),
   readdata          (char* file,
                      int *nfoods_p, double **cost_p,
                      double **lb_p, double **ub_p, 
                      int *nnutr_p, double **nutrmin_p, double **nutrmax_p,
                      double ***nutrper_p),
   populatebyrow     (CPXENVptr env, CPXLPptr lp,
                      int nfoods, double *cost, double *lb, double *ub, 
                      int nnutr, double *nutrmin, double *nutrmax,
                      double **nutrper),
   populatebycolumn  (CPXENVptr env, CPXLPptr lp,
                      int nfoods, double *cost, double *lb, double *ub, 
                      int nnutr, double *nutrmin, double *nutrmax,
                      double **nutrper);

static void
   free_and_null     (char **ptr),
   usage             (char *progname);



int
main (int argc, char **argv)
{
   int status = 0;

   int    nfoods;
   int    nnutr;
   double *cost     = NULL;
   double *lb       = NULL;
   double *ub       = NULL;
   double *nutrmin  = NULL;
   double *nutrmax  = NULL;
   double **nutrper = NULL;

   double *x = NULL;
   double objval;
   int    solstat;

   /* Declare and allocate space for the variables and arrays where we
      will store the optimization results including the status, objective
      value, variable values, dual values, row slacks and variable
      reduced costs. */

   CPXENVptr     env = NULL;
   CPXLPptr      lp = NULL;
   int           i, j;

   /* Check the command line arguments */

   if (( argc != 3 )                                        ||
       ( argv[1][0] != '-' )                                ||
       ( strchr ("rc", argv[1][1]) == NULL )  ) {
      usage (argv[0]);
      goto TERMINATE;
   }

   status = readdata(argv[2], &nfoods, &cost, &lb, &ub,
                     &nnutr, &nutrmin, &nutrmax, &nutrper);
   if ( status ) goto TERMINATE;


   /* Initialize the CPLEX environment */

   env = CPXopenCPLEX (&status);

   /* If an error occurs, the status value indicates the reason for
      failure.  A call to CPXgeterrorstring will produce the text of
      the error message.  Note that CPXopenCPLEX produces no output,
      so the only way to see the cause of the error is to use
      CPXgeterrorstring.  For other CPLEX routines, the errors will
      be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */

   if ( env == NULL ) {
      char  errmsg[CPXMESSAGEBUFSIZE];
      fprintf (stderr, "Could not open CPLEX environment.\n");
      CPXgeterrorstring (env, status, errmsg);
      fprintf (stderr, "%s", errmsg);
      goto TERMINATE;
   }

   /* Turn on output to the screen */

   status = CPXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON);
   if ( status ) {
      fprintf (stderr, 
               "Failure to turn on screen indicator, error %d.\n", status);
      goto TERMINATE;
   }

   /* Turn on data checking */

   status = CPXsetintparam (env, CPXPARAM_Read_DataCheck, CPX_ON);
   if ( status ) {
      fprintf (stderr, 
               "Failure to turn on data checking, error %d.\n", status);
      goto TERMINATE;
   }

   /* Create the problem. */

   lp = CPXcreateprob (env, &status, "diet");

   /* A returned pointer of NULL may mean that not enough memory
      was available or there was some other problem.  In the case of 
      failure, an error message will have been written to the error 
      channel from inside CPLEX.  In this example, the setting of
      the parameter CPXPARAM_ScreenOutput causes the error message to
      appear on stdout.  */

   if ( lp == NULL ) {
      fprintf (stderr, "Failed to create LP.\n");
      goto TERMINATE;
   }

   /* Now populate the problem with the data.  For building large
      problems, consider setting the row, column and nonzero growth
      parameters before performing this task. */

   switch (argv[1][1]) {
      case 'r':
         status = populatebyrow (env, lp, nfoods, cost, lb, ub, 
                                 nnutr, nutrmin, nutrmax, nutrper);
         break;
      case 'c':
         status = populatebycolumn (env, lp, nfoods, cost, lb, ub, 
                                    nnutr, nutrmin, nutrmax, nutrper);
         break;
   }

   if ( status ) {
      fprintf (stderr, "Failed to populate problem.\n");
      goto TERMINATE;
   }


   /* Optimize the problem and obtain solution. */

   status = CPXlpopt (env, lp);
   if ( status ) {
      fprintf (stderr, "Failed to optimize LP.\n");
      goto TERMINATE;
   }

   x = (double *) malloc (nfoods * sizeof(double));
   if ( x == NULL ) {
      status = CPXERR_NO_MEMORY;
      fprintf (stderr, "Could not allocate memory for solution.\n");
      goto TERMINATE;
   }

   status = CPXsolution (env, lp, &solstat, &objval, x, NULL, NULL, NULL);
   if ( status ) {
      fprintf (stderr, "Failed to obtain solution.\n");
      goto TERMINATE;
   }

   /* Write the output to the screen. */

   printf ("\nSolution status = %d\n", solstat);
   printf ("Solution value  = %f\n\n", objval);

   for (j = 0; j < nfoods; j++) {
      printf ("Food %d:  Buy = %10f\n", j, x[j]);
   }

   /* Finally, write a copy of the problem to a file. */

   status = CPXwriteprob (env, lp, "diet.lp", NULL);
   if ( status ) {
      fprintf (stderr, "Failed to write LP to disk.\n");
      goto TERMINATE;
   }

   
TERMINATE:

   /* Free up the problem as allocated by CPXcreateprob, if necessary */

   if ( lp != NULL ) {
      status = CPXfreeprob (env, &lp);
      if ( status ) {
         fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
      }
   }

   /* Free up the CPLEX environment, if necessary */

   if ( env != NULL ) {
      status = CPXcloseCPLEX (&env);

      /* Note that CPXcloseCPLEX produces no output,
         so the only way to see the cause of the error is to use
         CPXgeterrorstring.  For other CPLEX routines, the errors will
         be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON. */

      if ( status > 0 ) {
         char  errmsg[CPXMESSAGEBUFSIZE];
         fprintf (stderr, "Could not close CPLEX environment.\n");
         CPXgeterrorstring (env, status, errmsg);
         fprintf (stderr, "%s", errmsg);
      }
   }

   if ( nutrper != NULL ) {
      for (i = 0; i < nnutr; ++i) {
         free_and_null ((char **) &(nutrper[i]));
      }
   }
   free_and_null ((char **) &nutrper);
   free_and_null ((char **) &cost);
   free_and_null ((char **) &cost);
   free_and_null ((char **) &lb);
   free_and_null ((char **) &ub);
   free_and_null ((char **) &nutrmin);
   free_and_null ((char **) &nutrmax);
   free_and_null ((char **) &x);

   return (status);

}  /* END main */


static int
populatebyrow (CPXENVptr env, CPXLPptr lp,
               int nfoods, double *cost, double *lb, double *ub, 
               int nnutr, double *nutrmin, double *nutrmax,
               double **nutrper)
{
   int status = 0;

   int zero = 0;
   int *ind = NULL;
   int i, j;

   ind = (int*) malloc(nfoods * sizeof(int));
   if ( ind == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }
   for (j = 0; j < nfoods; j++) {
      ind[j] = j;
   }

   status = CPXnewcols (env, lp, nfoods, cost, lb, ub, NULL, NULL);
   if ( status )  goto TERMINATE;

   for (i = 0; i < nnutr; i++) {
      double rng  = nutrmax[i] - nutrmin[i];

      status = CPXaddrows (env, lp, 0, 1, nfoods, nutrmin+i, "R",
                           &zero, ind, nutrper[i], NULL, NULL);
      if ( status )  goto TERMINATE;

      status = CPXchgrngval (env, lp, 1, &i, &rng);
      if ( status )  goto TERMINATE;
   }

TERMINATE:

   free_and_null ((char **)&ind);

   return (status);

}  /* END populatebyrow */



/* To populate by column, we first create the rows, and then add the
   columns.  */

static int
populatebycolumn (CPXENVptr env, CPXLPptr lp,
                  int nfoods, double *cost, double *lb, double *ub, 
                  int nnutr, double *nutrmin, double *nutrmax,
                  double **nutrper)
{
   int status = 0;

   int i, j;

   int    zero    = 0;
   int    *ind    = NULL;
   double *val    = NULL;
   char   *sense  = NULL;
   double *rngval = NULL;

   sense = (char*)malloc(nnutr * sizeof(char));
   if ( sense == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }
   for (i = 0; i < nnutr; i++) {
      sense[i] = 'R';
   }

   val = (double*)malloc(nnutr * sizeof(double));
   if ( val == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }

   rngval = (double*)malloc(nnutr * sizeof(double));
   if ( rngval == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }
   for (i = 0; i < nnutr; i++) {
      rngval[i] = nutrmax[i] - nutrmin[i];
   }

   ind = (int*) malloc(nfoods * sizeof(int));
   if ( ind == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }
   for (i = 0; i < nnutr; i++) {
      ind[i] = i;
   }

   status = CPXnewrows (env, lp, nnutr, nutrmin, sense, rngval, NULL);
   if ( status )  goto TERMINATE;

   for (j = 0; j < nfoods; ++j) {
      for (i = 0; i < nnutr; i++) {
         val[i] = nutrper[i][j];
      }

      status = CPXaddcols (env, lp, 1, nnutr, cost+j, &zero,
                           ind, val, lb+j, ub+j, NULL);
      if ( status )  goto TERMINATE;
   }

TERMINATE:

   free_and_null ((char **)&sense);
   free_and_null ((char **)&rngval);
   free_and_null ((char **)&ind);
   free_and_null ((char **)&val);

   return (status);

}  /* END populatebycolumn */


/* This simple routine frees up the pointer *ptr, and sets *ptr to NULL */

static void
free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */  



static void
usage (char *progname)
{
   fprintf (stderr,"Usage: %s -X <datafile>\n", progname);
   fprintf (stderr,"   where X is one of the following options: \n");
   fprintf (stderr,"      r          generate problem by row\n");
   fprintf (stderr,"      c          generate problem by column\n");
   fprintf (stderr," Exiting...\n");
} /* END usage */


static int
readarray (FILE *in, int *num_p, double **data_p)
{
   int  status = 0;
   int  max, num;
   char ch;

   num = 0;
   max = 10;

   *data_p = (double*)malloc(max * sizeof(double));
   if ( *data_p == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }

   for (;;) {
      fscanf (in, "%c", &ch);
      if ( ch == '\t' ||
           ch == '\r' ||
           ch == ' '  ||
           ch == '\n'   ) continue;
      if ( ch == '[' ) break;
      status = -1;
      goto TERMINATE;
   }

   for(;;) {
      int read;
      read = fscanf (in, "%lg", (*data_p)+num);
      if ( read == 0 ) {
         status = -1;
         goto TERMINATE;
      }
      num++;
      if ( num >= max ) {
         max *= 2;
         *data_p = (double*)realloc(*data_p, max * sizeof(double));
         if ( *data_p == NULL ) {
            status = CPXERR_NO_MEMORY;
            goto TERMINATE;
         }
      }
      do {
         fscanf (in, "%c", &ch);
      } while (ch == ' ' || ch == '\n' || ch == '\t'  || ch == '\r');
      if ( ch == ']' ) break;
      else if ( ch != ',' ) {
         status = -1;
         goto TERMINATE;
      }
   }

   *num_p = num;

TERMINATE:

   return (status);

} /* END readarray */


static int
readdata (char* file,
          int *nfoods_p, double **cost_p, double **lb_p, double **ub_p, 
          int *nnutr_p, double **nutrmin_p, double **nutrmax_p,
          double ***nutrper_p)
{
   int status = 0;

   int ncost, nlb, nub;
   int nmin, nmax;

   int  i, n;
   char ch;
   FILE *in = NULL;

   in = fopen(file, "r");
   if ( in == NULL ) {
      status = -1;
      goto TERMINATE;
   }

   if ( (status = readarray(in, &ncost, cost_p)) ) goto TERMINATE;
   if ( (status = readarray(in, &nlb,   lb_p))   ) goto TERMINATE;
   if ( (status = readarray(in, &nub,   ub_p))   ) goto TERMINATE;
   if ( ncost != nlb  ||  ncost != nub ) {
      status = -1;
      goto TERMINATE;
   }
   *nfoods_p = ncost;

   if ( (status = readarray(in, &nmin, nutrmin_p)) ) goto TERMINATE;
   if ( (status = readarray(in, &nmax, nutrmax_p)) ) goto TERMINATE;
   if ( nmax != nmin ) {
      status = -1;
      goto TERMINATE;
   }
   *nnutr_p = nmin;

   *nutrper_p = (double**)malloc(nmin * sizeof(double*));
   if ( *nutrper_p == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }

   for (;;) {
      fscanf (in, "%c", &ch);
      if ( ch == '\t' ||
           ch == '\r' ||
           ch == ' '  ||
           ch == '\n'   ) continue;
      if ( ch == '[' ) break;
      status = -1;
      goto TERMINATE;
   }
   for ( i = 0; i < nmin; i++ ) {
      if ( (status = readarray(in, &n, (*nutrper_p)+i)) ) goto TERMINATE;
      if ( n != ncost ) {
         status = -1;
         goto TERMINATE;
      }
      fscanf (in, "%c", &ch);
      if ( i < nmin-1  &&  ch != ',' ) {
         status = -1;
         goto TERMINATE;
      }
   }
   if ( ch != ']' ) {
      status = -1;
      goto TERMINATE;
   }


TERMINATE:
   if ( in != NULL )
      fclose (in);

   return (status);

} /* END readdata */

