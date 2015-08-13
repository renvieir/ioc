/* --------------------------------------------------------------------------
 * File: DistMIPex2.java
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
 *
 * DistMIPex2.java - Reading a MIP problem from a file and solving
 *                   it via distributed parallel optimization using
 *                   an informational callback.
 * See the usage() function for details about how to run this.
 */

import ilog.concert.*;
import ilog.cplex.*;


public class DistMIPEx2 {
   private static void usage() {
      final String program = "java DistMIPEx2";
      System.err.println("Usage: " + program + " <vmc> <model>");
      System.err.println("   Solves a model specified by a model file using");
      System.err.println("   distributed parallel MIP.");
      System.err.println("   Arguments:");
      System.err.println("    <vmc>    The virtual machine configuration file");
      System.err.println("             that describes the machine that can be");
      System.err.println("             used for parallel distributed MIP.");
      System.err.println("    <model>  Model file with the model to be solved.");
      System.err.println("   Example:");
      System.err.println("     " + program + " process.vmc model.lp");
   }

   // Log new incumbents if they are at better than the old by a
   // relative tolerance of 1e-5; also log progress info every
   // 100 nodes.

   static class LogCallback extends IloCplex.MIPInfoCallback {
      IloNumVar[] _var;
      double      _lastDettime;
      double      _lastIncumbent;
      double      _startTime;
      double      _startDetTime;

      LogCallback(IloNumVar[] var, double lastDettime, double lastIncumbent,
		  double startTime, double startDetTime) {
         _var = var;
         _lastDettime = lastDettime;
         _lastIncumbent = lastIncumbent;
	 _startTime = startTime;
	 _startDetTime = startDetTime;
      }

      public void main() throws IloException {
         boolean newIncumbent = false;
         double  dettime      = getDetTime();

         if ( hasIncumbent()  &&
              Math.abs(_lastIncumbent - getIncumbentObjValue())
                 > 1e-5*(1.0 + Math.abs(getIncumbentObjValue())) ) {
            _lastIncumbent = getIncumbentObjValue();
            newIncumbent = true;
         }
         if ( dettime >= _lastDettime + 1000.0  ||  newIncumbent ) {
            if ( !newIncumbent )  _lastDettime = dettime;
            System.out.print("Time = " + (getCplexTime() - _startTime)
                             + "  Dettime = " +  (dettime - _startDetTime)
                             + "  Best objective = " + getBestObjValue());
            if ( hasIncumbent() ) {
               System.out.println ("  Incumbent objective = " +
                                   getIncumbentObjValue());
            }
            else {
               System.out.println("");
            }
         }

         if ( newIncumbent ) {
            System.out.println("New incumbent values: ");

            int n = _var.length;
            double[] x = getIncumbentValues(_var, 0, n);
            for (int j = 0; j < n; j++) {
               System.out.println("Variable " + j + ": Value = " + x[j]);
            }
         }
      }
   }

   public static void main(String[] args) {
      // Check command line length.
      if ( args.length != 2 ) {
         usage();
         System.exit(-1);
      }

      // Pick up VMC from command line.
      String vmconfig = args[0];

      // Solve the model.
      IloCplex cplex = null;
      try {
         // Create CPLEX solver and load model.
         cplex = new IloCplex();
         cplex.importModel(args[1]);

         // Load the virtual machine configuration.
         // This will force solve() to use distributed parallel MIP.
         cplex.readVMConfig(vmconfig);

	 // Install logging info callback.
         IloLPMatrix lp = (IloLPMatrix)cplex.LPMatrixIterator().next();
         IloObjective obj = cplex.getObjective();
	 double lastObjVal =
	     (obj.getSense() == IloObjectiveSense.Minimize ) ?
	     Double.MAX_VALUE : -Double.MAX_VALUE;
	 cplex.use(new LogCallback(lp.getNumVars(), -100000, lastObjVal,
				   cplex.getCplexTime(), cplex.getDetTime()));

	 // Turn off CPLEX logging
	 cplex.setParam(IloCplex.Param.MIP.Display, 0);

         // Solve the model and print some solution information.
         if ( cplex.solve() )
            System.out.println("Solution value  = " + cplex.getObjValue());
         else
            System.out.println("No solution available");
         System.out.println("Solution status = " + cplex.getStatus());
      }
      catch (IloException e) {
         System.err.println("Concert exception caught '" + e + "' caught");
      } finally {
         if ( cplex != null )
            cplex.end();
      }
   }
}
