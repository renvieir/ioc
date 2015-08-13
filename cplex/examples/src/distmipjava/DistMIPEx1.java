/* --------------------------------------------------------------------------
 * File: DistMIPex1.java
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
 * DistMIPex1.java - Reading a MIP problem from a file and solving
 *                   it via distributed parallel optimization.
 * See the usage() function for details about how to run this.
 */

import ilog.concert.*;
import ilog.cplex.*;


public class DistMIPEx1 {
   private static void usage() {
      final String program = "java DistMIPEx1";
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
