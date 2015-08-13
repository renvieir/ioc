// -------------------------------------------------------------- -*- C++ -*-
// File: ilomdistmipex2.cpp
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
// ilodistmipex2.cpp - Reading a MIP problem from a file and solving
//                     it via distributed parallel optimization using
//                     an informational callback.
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

// Log new incumbents if they are at better than the old by a
// relative tolerance of 1e-5; also log progress info every
// 100 nodes.

ILOMIPINFOCALLBACK5(loggingCallback,
                    IloNumVarArray, vars,
                    IloNum,         lastDettime, 
                    IloNum,         lastIncumbent,
                    IloNum,         startTime,
                    IloNum,         startDetTime)
{
   int newIncumbent = 0;
   double dettime = getDetTime();

   if ( hasIncumbent()                                  && 
        fabs(lastIncumbent - getIncumbentObjValue())
              > 1e-5*(1.0 + fabs(getIncumbentObjValue())) ) {
      lastIncumbent = getIncumbentObjValue();
      newIncumbent = 1;
   }
     
   if ( dettime >= lastDettime + 1000.0  ||  newIncumbent ) {  

      if ( !newIncumbent )  lastDettime = dettime;
      getEnv().out() << "Time = " << getCplexTime() - startTime
                     << "  Dettime = " << dettime - startDetTime
                     << "  Best objective = " << getBestObjValue();

      if ( hasIncumbent() ) {
         getEnv().out() << "  Incumbent objective = " << getIncumbentObjValue()
                        << endl;
      }
      else {
         getEnv().out() << endl;
      }

   }
   if ( newIncumbent ) {
      IloNumArray val(vars.getEnv());
      getIncumbentValues(val, vars);
      val[0] = getIncumbentValue(vars[0]);
      getEnv().out() << "New incumbent variable values: " << endl
                     << val << endl;
      val.end();
   }
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

      IloObjective   obj;
      IloNumVarArray var(env);
      IloRangeArray  rng(env);
      IloSOS1Array   sos1(env);
      IloSOS2Array   sos2(env);
      IloRangeArray  lazy(env);
      IloRangeArray  cuts(env);

      cplex.importModel(model, argv[2], obj, var, rng, sos1, sos2,
                        lazy, cuts);

      cplex.extract(model);

      if ( lazy.getSize() > 0 )  cplex.addLazyConstraints (lazy);
      if ( cuts.getSize() > 0 )  cplex.addUserCuts (cuts);

      // Load the virtual machine configuration.
      // This will force solve() to use parallel distributed MIP.
      cplex.readVMConfig(vmconfig);

      // Install logging info callback.
      IloNum lastObjVal = (obj.getSense() == IloObjective::Minimize ) ?
         IloInfinity : -IloInfinity;
      cplex.use(loggingCallback(env, var, -100000, lastObjVal,
                                cplex.getCplexTime(), cplex.getDetTime()));
      // Turn off CPLEX logging
      cplex.setParam(IloCplex::Param::MIP::Display, 0);


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
