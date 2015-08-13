// -------------------------------------------------------------- -*- C++ -*-
// File: facility.cpp
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

#include <ilcplex/ilocplex.h>

ILOSTLBEGIN

typedef IloArray<IloNumArray>    FloatMatrix;
typedef IloArray<IloNumVarArray> NumVarMatrix;

int
main(int argc, char **argv)
{
   IloEnv env;
   try {
      IloInt i, j;

      IloNumArray capacity(env), fixedCost(env);
      FloatMatrix cost(env);
      IloInt      nbLocations;
      IloInt      nbClients;

      const char* filename  = "../../../examples/data/facility.dat";
      if (argc > 1) 
         filename = argv[1];
      ifstream file(filename);
      if (!file) {
         cerr << "ERROR: could not open file '" << filename
              << "' for reading" << endl;
         cerr << "usage:   " << argv[0] << " <file>" << endl;
         throw(-1);
      }

      file >> capacity >> fixedCost >> cost;
      nbLocations = capacity.getSize();
      nbClients   = cost.getSize(); 

      IloBool consistentData = (fixedCost.getSize() == nbLocations);
      for(i = 0; consistentData && (i < nbClients); i++)
         consistentData = (cost[i].getSize() == nbLocations);
      if (!consistentData) {
         cerr << "ERROR: data file '" 
              << filename << "' contains inconsistent data" << endl;
         throw(-1);
      }

      IloNumVarArray open(env, nbLocations, 0, 1, ILOINT);
      NumVarMatrix supply(env, nbClients);
      for(i = 0; i < nbClients; i++)
         supply[i] = IloNumVarArray(env, nbLocations, 0, 1, ILOINT);

      IloModel model(env);
      for(i = 0; i < nbClients; i++)
         model.add(IloSum(supply[i]) == 1);
      for(j = 0; j < nbLocations; j++) {
         IloExpr v(env);
         for(i = 0; i < nbClients; i++)
            v += supply[i][j];
         model.add(v <= capacity[j] * open[j]);
         v.end();
      }

      IloExpr obj = IloScalProd(fixedCost, open);
      for(i = 0; i < nbClients; i++) {
         obj += IloScalProd(cost[i], supply[i]);
      }
      model.add(IloMinimize(env, obj));
      obj.end();

      IloCplex cplex(env);
      cplex.extract(model);
      cplex.solve();
	
      cplex.out() << "Solution status: " << cplex.getStatus() << endl;

      IloNum tolerance = cplex.getParam(
         IloCplex::Param::MIP::Tolerances::Integrality);
      cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
      for(j = 0; j < nbLocations; j++) {
         if (cplex.getValue(open[j]) >= 1 - tolerance) {
            cplex.out() << "Facility " << j << " is open, it serves clients ";
            for(i = 0; i < nbClients; i++) {
               if (cplex.getValue(supply[i][j]) >= 1 - tolerance)
                  cplex.out() << i << " ";
            }
            cplex.out() << endl; 
         }
      }
   }
   catch(IloException& e) {
      cerr  << " ERROR: " << e << endl;   
   }
   catch(...) {
      cerr  << " ERROR" << endl;   
   }
   env.end();
   return 0;
}
