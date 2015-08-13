/* --------------------------------------------------------------------------
 * File: xdistmipex1.c
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

/* xdistmipex1.c - Reading a MIP problem from a file and solving
 *                 it via distributed parallel optimization.
 * See the usage() function for details about how to run this.
 */
#ifdef USE_MPI
#   include <mpi.h>
#endif

#include <ilcplex/cplexx.h>
#include <ilcplex/cplexdistmipx.h>

#include <string.h>
#include <stdlib.h>


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
   /* Enable screen output so that we can see any progress or
    * error messages.
    */
   status = CPXXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
   if ( status != 0 ) {
      fprintf (stderr, "Failed to enable screen output: %s\n",
               CPXXgeterrorstring (env, status, errbuf));
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
