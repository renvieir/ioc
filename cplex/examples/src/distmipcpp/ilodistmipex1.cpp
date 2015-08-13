// -------------------------------------------------------------- -*- C++ -*-
// File: ilomdistmipex1.cpp
// Version 12.6.1
// --------------------------------------------------------------------------
// Licensed Materials - Property of IBM
// 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
// Copyright IBM Corporation 2013, 2014. All Rights Reserved.
//
// US Government Users Restricted Rights - Use, duplication or
// disclosure restricted by GSA ADP Schedule Contract with
// IBM Corp.
// --------------------------------------------------------------------------
//
// ilodistmipex1.cpp - Reading a MIP problem from a file and solving
//                     it via distributed parallel optimization.
// See the usage() function for details about how to run this.

#include <ilcplex/ilocplex.h>
#include <cstring>

ILOSTLBEGIN

static void
usage (char const *program)
{
   cerr << "Usage: " << program << " <vmc> <model>" << endl
        << "   Solves a model specified by a model file using " << endl
        << "   distributed parallel MIP." << endl
        << "   Arguments:" << endl
        << "    <vmc>    The virtual machine configuration file" << endl
        << "             that describes the machine that can be" << endl
        << "             used for parallel distributed MIP." << endl
        << "    <model>  Model file with the model to be solved." << endl
        << "   Example:" << endl
        << "     " << program << " process.vmc model.lp" << endl;
}

int
main (int argc, char **argv)
{
   char const *vmconfig = NULL;

   // Check command line length (exactly two arguments are required).
   if ( argc != 3 ) {
      usage (argv[0]);
      return -1;
   }

   // Pick up VMC from command line.
   vmconfig = argv[1];

   // Solve the model.
   int exitcode = 0;
   IloEnv env;
   try {
      // Create and read the model.
      IloModel model(env);
      IloCplex cplex(model);
      cplex.importModel(model, argv[2]);

      // Load the virtual machine configuration.
      // This will force solve() to use parallel distributed MIP.
      cplex.readVMConfig(vmconfig);

      // Solve the problem and display some results.
      if ( cplex.solve() )
         env.out() << "Solution value  = " << cplex.getObjValue() << endl;
      else
         env.out() << "No solution" << endl;
      env.out() << "Solution status = " << cplex.getStatus() << endl;

      // Cleanup.
      cplex.end();
      model.end();
   }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
      exitcode = -1;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
      exitcode = -1;
   }

   env.end();

   return exitcode;

}  // END main
