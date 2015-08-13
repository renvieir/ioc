/* --------------------------------------------------------------------------
 * File: Diet.java   
 * Version 12.6.1  
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
 * Copyright IBM Corporation 2001, 2014. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 * --------------------------------------------------------------------------
 *
 * A dietary model.
 *
 * Input data:
 * foodMin[j]          minimum amount of food j to use
 * foodMax[j]          maximum amount of food j to use 
 * foodCost[j]         cost for one unit of food j
 * nutrMin[i]          minimum amount of nutrient i
 * nutrMax[i]          maximum amount of nutrient i
 * nutrPerFood[i][j]   nutrition amount of nutrient i in food j
 *
 * Modeling variables:
 * Buy[j]          amount of food j to purchase
 *
 * Objective:
 * minimize sum(j) Buy[j] * foodCost[j]
 *
 * Constraints:
 * forall foods i: nutrMin[i] <= sum(j) Buy[j] * nutrPer[i][j] <= nutrMax[j]
 */

import ilog.concert.*;
import ilog.cplex.*;


public class Diet {
   static class Data {
      int        nFoods;
      int        nNutrs;
      double[]   foodCost;
      double[]   foodMin;
      double[]   foodMax;
      double[]   nutrMin;
      double[]   nutrMax;
      double[][] nutrPerFood; 
    
      Data(String filename) throws IloException, java.io.IOException,
                                   InputDataReader.InputDataReaderException {
         InputDataReader reader = new InputDataReader(filename);
         
         foodCost = reader.readDoubleArray();
         foodMin  = reader.readDoubleArray();
         foodMax  = reader.readDoubleArray();
         nutrMin  = reader.readDoubleArray();
         nutrMax  = reader.readDoubleArray();
         nutrPerFood = reader.readDoubleArrayArray();
       
         nFoods = foodMax.length;
         nNutrs = nutrMax.length;
       
         if ( nFoods != foodMin.length  ||
              nFoods != foodMax.length    )
            throw new IloException("inconsistent data in file " + filename);
         if ( nNutrs != nutrMin.length    ||
              nNutrs != nutrPerFood.length  )
            throw new IloException("inconsistent data in file " + filename);
         for (int i = 0; i < nNutrs; ++i) {
            if ( nutrPerFood[i].length != nFoods )
               throw new IloException("inconsistent data in file " + filename);
         }
      }
   }

   static void buildModelByRow(IloModeler    model,
                               Data          data,
                               IloNumVar[]   Buy,
                               IloNumVarType type) throws IloException {
      int nFoods = data.nFoods;
      int nNutrs = data.nNutrs;

      for (int j = 0; j < nFoods; j++) {
         Buy[j] = model.numVar(data.foodMin[j], data.foodMax[j], type);
      }
      model.addMinimize(model.scalProd(data.foodCost, Buy));

      for (int i = 0; i < nNutrs; i++) {
         model.addRange(data.nutrMin[i],
                        model.scalProd(data.nutrPerFood[i], Buy),
                        data.nutrMax[i]);
      }
   }

   static void buildModelByColumn(IloMPModeler  model,
                                  Data          data,
                                  IloNumVar[]   Buy,
                                  IloNumVarType type) throws IloException {
      int nFoods = data.nFoods;
      int nNutrs = data.nNutrs;

      IloObjective cost       = model.addMinimize();
      IloRange[]   constraint = new IloRange[nNutrs];
    
      for (int i = 0; i < nNutrs; i++) {
         constraint[i] = model.addRange(data.nutrMin[i], data.nutrMax[i]);
      }
   
      for (int j = 0; j < nFoods; j++) {
         IloColumn col = model.column(cost, data.foodCost[j]);
         for (int i = 0; i < nNutrs; i++) {
            col = col.and(model.column(constraint[i], data.nutrPerFood[i][j]));
         }
         Buy[j] = model.numVar(col, data.foodMin[j], data.foodMax[j], type);
      }
   }


   public static void main(String[] args) {

      try {
         String          filename  = "../../../examples/data/diet.dat";
         boolean         byColumn  = false;
         IloNumVarType   varType   = IloNumVarType.Float;
        
         for (int i = 0; i < args.length; i++) {
            if ( args[i].charAt(0) == '-') {
               switch (args[i].charAt(1)) {
               case 'c':
                  byColumn = true;
                  break;
               case 'i':
                  varType = IloNumVarType.Int;
                  break;
               default:
                  usage();
                  return;
               }
            }
            else {
               filename = args[i];
               break;
            }
         }
        
         Data data = new Data(filename);
         int nFoods = data.nFoods;
       
         // Build model
         IloCplex     cplex = new IloCplex();
         IloNumVar[]  Buy   = new IloNumVar[nFoods];
       
         if ( byColumn ) buildModelByColumn(cplex, data, Buy, varType);
         else            buildModelByRow   (cplex, data, Buy, varType);

         // Solve model
       
         if ( cplex.solve() ) { 
            System.out.println();
            System.out.println("Solution status = " + cplex.getStatus());
            System.out.println();
            System.out.println(" cost = " + cplex.getObjValue());
            for (int i = 0; i < nFoods; i++) {
               System.out.println(" Buy" + i + " = " + cplex.getValue(Buy[i]));
            }
            System.out.println();
         }
         cplex.end();
      }
      catch (IloException ex) {
         System.out.println("Concert Error: " + ex);
      }
      catch (InputDataReader.InputDataReaderException ex) {
         System.out.println("Data Error: " + ex);
      }
      catch (java.io.IOException ex) {
         System.out.println("IO Error: " + ex);
      }
   }

   static void usage() {
      System.out.println(" ");
      System.out.println("usage: java Diet [options] <data file>");
      System.out.println("options: -c  build model by column");
      System.out.println("         -i  use integer variables");
      System.out.println(" ");
   }
}

/*  Sample output

Solution status = Optimal

cost   = 14.8557
  Buy0 = 4.38525
  Buy1 = 0
  Buy2 = 0
  Buy3 = 0
  Buy4 = 0
  Buy5 = 6.14754
  Buy6 = 0
  Buy7 = 3.42213
  Buy8 = 0
*/
