/* --------------------------------------------------------------------------
 * File: xdistmipex2.c
 * Version 12.6.1
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
 * Copyright IBM Corporation 2013, 2014. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 * --------------------------------------------------------------------------
 */

/* xdistmipex2.c - Reading a MIP problem from a file and solving
 *                 it via distributed parallel optimization and
 *                 using the informational callback.
 * See the usage() function for details about how to run this.
 */
#ifdef USE_MPI
#   include <mpi.h>
#endif

#include <ilcplex/cplexx.h>
#include <ilcplex/cplexdistmipx.h>

#include <string.h>
#include <stdlib.h>
#include <math.h>

struct loginfo {
   double timestart;
   double dettimestart;
   double lastincumbent;
   double lastdettime;
   CPXDIM numcols;
};
typedef struct loginfo LOGINFO, *LOGINFOptr;

static int CPXPUBLIC
   logcallback (CPXCENVptr env, void *cbdata, int wherefrom,
                void *cbhandle);

static void
usage (char const *program)
{
   fprintf (stderr, "Usage: %s <vmc> <model>\n"
            "   Solves a model specified by a model file using\n"
            "   distributed parallel MIP.\n"
            "   Arguments:\n"
            "    <vmc>    The virtual machine configuration file\n"
            "             that describes the machine that can be\n"
            "             used for parallel distributed MIP.\n"
            "    <model>  Model file with the model to be solved.\n"
            "   Example:\n"
            "     %s  process.vmc model.lp\n",
            program, program);
}

int
main (int argc, char **argv)
{
   CPXENVptr env = NULL;
   CPXLPptr lp = NULL;
   char errbuf[CPXMESSAGEBUFSIZE];
   int status = 0;
   char const *vmconfig = NULL;
   char const *model = NULL;
   double objval;
   LOGINFO myloginfo;

#ifdef USE_MPI
   MPI_Init (&argc, &argv);
#endif

   if ( argc != 3 ) {
      usage (argv[0]);
      return -1;
   }
   
   vmconfig = argv[1];
   model = argv[2];

   /* Create a new CPLEX environment for each problem to solve.
    */
   env = CPXXopenCPLEX (&status);
   if ( env == NULL || status != 0 ) {
      fprintf (stderr, "Failed to open CPLEX: %s\n",
               CPXXgeterrorstring (NULL, status, errbuf));
      goto TERMINATE;
   }
   /* Load a virtual machine configuration.
    */
   status = CPXXreadcopyvmconfig (env, vmconfig);
   if ( status != 0 ) {
      fprintf (stderr, "Failed to load VMC %s: %s\n",
               vmconfig, CPXXgeterrorstring (env, status, errbuf));
      goto TERMINATE;
   }
   /* Create and input a problem.
    */
   lp = CPXXcreateprob (env, &status, model);
   if ( lp == NULL || status != 0 ) {
      fprintf (stderr, "Failed to create problem: %s\n",
               CPXXgeterrorstring (env, status, errbuf));
      goto TERMINATE;
   }
   status = CPXXreadcopyprob (env, lp, model, NULL);
   if ( status != 0 ) {
      fprintf (stderr, "Failed to read problem %s: %s\n",
               model, CPXXgeterrorstring (env, status, errbuf));
      goto TERMINATE;
   }
   /* Turn off CPLEX logging.
    */
   status = CPXXsetintparam (env, CPXPARAM_MIP_Display, 0);
   if ( status )  goto TERMINATE;
   /* Install an incumbent callback for logging.
    */
   status = CPXXgettime (env, &myloginfo.timestart);
   if ( status ) {
      fprintf (stderr, "Failed to query time.\n");
      goto TERMINATE;
   }
   status = CPXXgetdettime (env, &myloginfo.dettimestart);
   if ( status ) {
      fprintf (stderr, "Failed to query deterministic time.\n");
      goto TERMINATE;
   }
   myloginfo.numcols       = CPXXgetnumcols (env, lp);
   myloginfo.lastincumbent = CPXXgetobjsen (env, lp) * 1e+35;
   myloginfo.lastdettime   = -10000.0;
   status = CPXXsetinfocallbackfunc (env, logcallback, &myloginfo);
   if ( status ) {
      fprintf (stderr, "Failed to set logging callback function.\n");
      goto TERMINATE;
   }
   /* Solve the problem using parallel distributed MIP.
    */
   status = CPXXdistmipopt (env, lp);
   if ( status != 0 ) {
      fprintf (stderr, "Failed to optimize: %s\n",
               CPXXgeterrorstring (env, status, errbuf));
      goto TERMINATE;
   }
   /* Print some solution information.
    */
   status = CPXXgetobjval (env, lp, &objval);
   if ( status == CPXERR_NO_SOLN ) {
      printf ("No solution available.\n");
   }
   else if ( status == 0 ) {
      printf ("Solution value %f\n", objval);
   }
   else {
      fprintf (stderr, "Error %d: %s\n", status,
               CPXXgeterrorstring (env, status, errbuf));
   }
   printf ("Solution status %d\n", CPXXgetstat (env, lp));


 TERMINATE:
   CPXXfreeprob (env, &lp);
   CPXXcloseCPLEX (&env);

#ifdef USE_MPI
   CPXXfinalizeMPIworkers (-1, NULL, 0, NULL, 1);
   MPI_Finalize ();
#endif

   return status;
}


/* This simple routine frees up the pointer *ptr, and sets *ptr to NULL */
static void
free_and_null (char **ptr)
{
   if ( *ptr != NULL ) {
      free (*ptr);
      *ptr = NULL;
   }
} /* END free_and_null */

/* Log new incumbents if they are at better than the old by a
 * relative tolerance of 1e-5; also log progress info every
 * 100 nodes.
 */
static int CPXPUBLIC
logcallback (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle)
{
   int status = 0;

   LOGINFOptr info = (LOGINFOptr) cbhandle;
   int        hasincumbent = 0;
   int        newincumbent = 0;
   double     dettime;
   double     objval;
   double     bound;
   double     *x = NULL;

   status = CPXXgetcallbackinfo (env, cbdata, wherefrom,
                                CPX_CALLBACK_INFO_MIP_FEAS, &hasincumbent);
   if ( status )  goto TERMINATE;

   if ( hasincumbent ) {
      status = CPXXgetcallbackinfo (env, cbdata, wherefrom,
                                   CPX_CALLBACK_INFO_BEST_INTEGER, &objval);
      if ( status )  goto TERMINATE;

      if ( fabs(info->lastincumbent - objval) > 1e-5*(1.0 + fabs(objval)) ) {
         newincumbent = 1;
         info->lastincumbent = objval;
      }
   }

   status = CPXXgetdettime (env, &dettime);
   if ( status )  goto TERMINATE;

   if ( dettime >= info->lastdettime + 1000.0  ||  newincumbent ) {
      double walltime;

      status = CPXXgetcallbackinfo (env, cbdata, wherefrom,
                                   CPX_CALLBACK_INFO_BEST_REMAINING, &bound);
      if ( status )  goto TERMINATE;

      if ( !newincumbent )  info->lastdettime = dettime;

      status = CPXXgettime (env, &walltime);
      if ( status )  goto TERMINATE;

      printf ("Time = %.2f  Dettime = %.2f  Best objective = %g",
              walltime - info->timestart, dettime - info->dettimestart, bound);
      if ( hasincumbent )  printf ("  Incumbent objective = %g\n", objval);
      else                 printf ("\n");
   }

   if ( newincumbent ) {
      int j;
      CPXDIM numcols = info->numcols;

      x = malloc (numcols*sizeof(*x));
      if ( x == NULL ) {
         status = CPXERR_NO_MEMORY;
         goto TERMINATE;
      }
      status = CPXXgetcallbackincumbent (env, cbdata, wherefrom,
                                        x, 0, numcols-1);
      if ( status )  goto TERMINATE;

      printf ("New incumbent values:\n");
      for (j = 0; j < numcols; j++) {
         if ( fabs(x[j]) > 1e-6 ) {
            printf ("  Column %d:  %g\n", j, x[j]);
         }
      }
   }

TERMINATE:

   free_and_null ((char **) &x);
   return (status);

} /* END logcallback */
