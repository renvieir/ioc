// -------------------------------------------------------------- -*- C++ -*-
// File: ilodiet.cpp
// Version 12.6.1  
// --------------------------------------------------------------------------
// Licensed Materials - Property of IBM
// 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
// Copyright IBM Corporation 2000, 2014. All Rights Reserved.
//
// US Government Users Restricted Rights - Use, duplication or
// disclosure restricted by GSA ADP Schedule Contract with
// IBM Corp.
// --------------------------------------------------------------------------
//
// A dietary model.
//
// Input data:
// foodMin[j]      minimum amount of food j to use
// foodMax[j]      maximum amount of food j to use 
// foodCost[j]     cost for one unit of food j
// nutrMin[i]      minimum amount of nutrient i
// nutrMax[i]      maximum amount of nutrient i
// nutrPer[i][j]   nutrition amount of nutrient i in food j
//
// Modeling variables:
// Buy[j]          amount of food j to purchase
//
// Objective:
// minimize sum(j) Buy[j] * foodCost[j]
//
// Constraints:
// forall foods i: nutrMin[i] <= sum(j) Buy[j] * nutrPer[i][j] <= nutrMax[j]
//

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN



void usage(const char* name) {
   cerr << endl;
   cerr << "usage:   " << name << " [options] <file>" << endl;
   cerr << "options: -c  build model by column" << endl;
   cerr << "         -i  use integer variables" << endl;
   cerr << endl;
}


void
buildModelByRow(IloModel           mod,
                IloNumVarArray     Buy,
                const IloNumArray  foodMin,
                const IloNumArray  foodMax,
                const IloNumArray  foodCost, 
                const IloNumArray  nutrMin,
                const IloNumArray  nutrMax,
                const IloNumArray2 nutrPer,
                IloNumVar::Type type) {
   IloEnv env = mod.getEnv();

   Buy.clear();
   IloNumVarArray tmp(env, foodMin, foodMax, type);
   Buy.add(tmp);
   tmp.end();

   IloInt i, j;
   IloInt n = foodCost.getSize();
   IloInt m = nutrMin.getSize();

   mod.add(IloMinimize(env, IloScalProd(Buy,foodCost)));
   for (i = 0; i < m; i++) {
      IloExpr expr(env);
      for (j = 0; j < n; j++) {
         expr += Buy[j] * nutrPer[i][j];
      }
      mod.add(nutrMin[i] <= expr <= nutrMax[i]);
      expr.end();
   }
}


void
buildModelByColumn(IloModel           mod,
                   IloNumVarArray     Buy,
                   const IloNumArray  foodMin,
                   const IloNumArray  foodMax,
                   const IloNumArray  foodCost, 
                   const IloNumArray  nutrMin,
                   const IloNumArray  nutrMax,
                   const IloNumArray2 nutrPer,
                   IloNumVar::Type    type) {
   IloEnv env = mod.getEnv();

   IloInt i, j;
   IloInt n = foodCost.getSize();
   IloInt m = nutrMin.getSize();

   IloRangeArray range (env, nutrMin, nutrMax);
   mod.add(range);
   IloObjective cost = IloAdd(mod, IloMinimize(env));

   for (j = 0; j < n; j++) {
      IloNumColumn col = cost(foodCost[j]);
      for (i = 0; i < m; i++) {
         col += range[i](nutrPer[i][j]);
      }
      Buy.add(IloNumVar(col, foodMin[j], foodMax[j], type));
      col.end();
   }
   range.end();
}


int
main(int argc, char **argv)
{
   IloEnv env;

   try {
      const char*     filename  = "../../../examples/data/diet.dat";
      IloBool         byColumn  = IloFalse;
      IloNumVar::Type varType   = ILOFLOAT;

      IloInt i;

      for (i = 1; i < argc; i++) {
         if (argv[i][0] == '-') {
            switch (argv[i][1]) {
            case 'c':
               byColumn = IloTrue;
               break;
            case 'i':
               varType = ILOINT;
               break;
            default:
               usage(argv[0]);
               throw (-1);
            }
         }
         else {
            filename = argv[i];
            break;
         }
      }

      ifstream file(filename);
      if ( !file ) {
         cerr << "ERROR: could not open file '" << filename
              << "' for reading" << endl;
         usage(argv[0]);
         throw (-1);
      }

      // model data

      IloNumArray  foodCost(env), foodMin(env), foodMax(env);
      IloNumArray  nutrMin(env), nutrMax(env);
      IloNumArray2 nutrPer(env);

      file >> foodCost >> foodMin >> foodMax;
      file >> nutrMin >> nutrMax;
      file >> nutrPer;

      IloInt nFoods = foodCost.getSize();
      IloInt nNutr  = nutrMin.getSize();

      if ( foodMin.getSize() != nFoods ||
           foodMax.getSize() != nFoods ||
           nutrPer.getSize() != nNutr  ||
           nutrMax.getSize() != nNutr    ) {
         cerr << "ERROR: Data file '" << filename
              << "' contains inconsistent data" << endl;
         throw (-1);
      }

      for (i = 0; i < nNutr; i++) {
         if (nutrPer[i].getSize() != nFoods) {
            cerr << "ERROR: Data file '" << argv[0]
                 << "' contains inconsistent data" << endl;
            throw (-1);
         }
      }

      // Build model

      IloModel       mod(env);
      IloNumVarArray Buy(env);
      if ( byColumn ) {
         buildModelByColumn(mod, Buy, foodMin, foodMax, foodCost,
                            nutrMin, nutrMax, nutrPer, varType);
      }
      else {
         buildModelByRow(mod, Buy, foodMin, foodMax, foodCost,
                         nutrMin, nutrMax, nutrPer, varType);
      }

      // Solve model

      IloCplex cplex(mod);
      cplex.exportModel("diet.lp");

      cplex.solve();
      cplex.out() << "solution status = " << cplex.getStatus() << endl;

      cplex.out() << endl;
      cplex.out() << "cost   = " << cplex.getObjValue() << endl;
      for (i = 0; i < foodCost.getSize(); i++) {
         cplex.out() << "  Buy" << i << " = " << cplex.getValue(Buy[i]) << endl;
      }
   }
   catch (IloException& ex) {
      cerr << "Error: " << ex << endl;
   }
   catch (...) {
      cerr << "Error" << endl;
   }

   env.end();

   return 0;
}

/*
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
