/* --------------------------------------------------------------------------
 * File: iloparbenders.cpp
 * Version 12.6.1
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
 * Copyright IBM Corporation, 2012, 2014. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 * --------------------------------------------------------------------------
 */
 
/* ********************************************************************** *
 *                                                                        *
 *    Parallel distributed Benders decomposition                          *
 *                                                                        *
 *    This file illustrates how a distributed Benders' decomposition      *
 *    can be implemented by means of the CPLEX remote object.             *
 *                                                                        *
 * ********************************************************************** */

/* The example is focused on a MILP model for the facility location problem,
 * but can be easily extended to any MILP model.
 *
 * The MILP model used in this example is
 *   minimize sum(f in F) sum(j in J) C[f,j]*x[f,j] + sum(f in F) F[f]*y[f]
 *     such that    
 *                  sum(f in F) x[f,j] >= 1         forall j in J
 *                   sum (f in F) y[f] <= 2
 *                                y[f] >= x[f,j]    forall f in F, j in J
 *                              x[f,j] >= 0         forall f in F, j in J
 *                                y[f] in { 0, 1 }  forall f in F
 * where y[f] = 1 if factory f is opened and 0 otherwise and
 * x[f,j] is the fraction of demand serviced from factory f to customer j.
 * F[f] is the fixed charge cost for opening factory f and C[f,j] is the
 * cost for servicing customer j from f.
 *
 * The model is decomposed in a master block involving the y binary variables
 * only, and in |J| blocks B[j], one for each customer j in J, that are used
 * to separate Bender's cuts.
 *
 * The master block is solved on a local machine, by adding Benders' cuts
 * through a lazyconstraint callback, while the blocks B[j] (j in J) are
 * solved in parallel on different remote machines, and are used to separate
 * violated Benders' cuts to be added to the current master.
 *
 * The MILP master block is:
 * minimize sum(j in J) eta[j] + sum(f in F) F[f]*y[f]
 *     such that    
 *                   sum (f in F) y[f] <= 2
 *                                y[f] in { 0, 1 }  forall f in F
 *                              eta[j] >= 0         for all j in J
 * where, for each j in J, the variable eta[j] = sum(f in F) C[f,j]*x[f,j]
 * takes into account the objective function associated with variables
 * x[f,j] (f in F)
 *
 * The LP block B[j] associated with each fixed customer j in J is
 * minimize sum(f in F) C[f,j]*x[f,j]
 *     such that    
 *                  sum(f in F) x[f,j] >= 1         
 *                              x[f,j] <= y[f]      forall f in F
 *                              x[f,j] >= 0         forall f in F
 * where the binary variables y are fixed and their value is provided
 * by the current solution of the master block.
 *
 * Given the current solution (y[F], eta[j]) of the master, for each customer
 * j in J a Bender's cut is separated by solving the dual of B[j], that is,
 * the following LP D[j]:
 * maximize beta - sum(f in F) y[f]*alpha[f]
 *     such that    
 *                     beta - alpha[f] <= C[f,j]    forall f in F
 *                            alpha[f] >= 0         forall f in F
 *                                beta >= 0         
 *
 * An unbounded ray (beta, alpha[f]) of D[j] gives raise to the ray cut
 * (feasibility cut)
 *    beta - sum(f in F) alpha[f]*y[f] <= 0
 * in the master problem.
 * A vertex (beta, alpha[f]) of D[j] give raises to the point cut
 * (optimality cut)
 *    eta[j] >= beta - sum(f in F) alpha[f]*y[f]
 * in the master problem.
 */
#include <cmath>
#include <iostream>
#include <sstream>

#if defined(USE_MPI)
#   define OMPI_SKIP_MPICXX 1 /* We don't use the C++ bindings of OpenMPI. */
#   include <mpi.h>
#   define TRANSPORT "mpitransport"
#elif defined(USE_PROCESS)
#   define TRANSPORT "processtransport"
#elif defined(USE_TCPIP)
#   define TRANSPORT "tcpiptransport"
#else
#   error "No transport type selected"
#endif


// The epsilon we use in comparion of double precision values
#define EPSILON 1e-6

/* ********************************************************************** *
 *                                                                        *
 *    M a s t e r   s i d e   i m p l e m e n t a t i o n                 *
 *                                                                        *
 * ********************************************************************** */

#if defined(COMPILE_MASTER)

#include <set>
#include <map>
#include <vector>
#include <ilcplex/ilocplex.h>

// ----------------------------------------------------------------------

/** A comparator for extractables.
 * All sub-classes of IloExtractable have a getId() member function that
 * returns a unique id for the extractable. We can use this id to compare
 * extractables that are used as keys in maps or in sets.
 */
template<typename T>
struct ExtractableLess {
   bool operator() (T const &e1, T const &e2) const {
      return e1.getId() < e2.getId();
   }
};
/** A set of rows. */
typedef std::set<IloRange,ExtractableLess<IloRange> > RowSet;

// ----------------------------------------------------------------------

/** A class that solves problems by means of parallel distributed Benders
 * decomposition. The public API of this class consists of only two things:
 * - An abstract (interface) class Problem that describes problems that
 *   can be solved by this class.
 * - A function solve() that solves instances of Problem.
 */
struct BendersOpt {
   /** Interface to problems that can be solved by the BendersOpt class. */
   struct Problem {
      Problem() {}
      virtual ~Problem() {}
      /** Get the number of blocks in this problem. */
      virtual IloInt getNBlocks() const = 0;
      /** Get the model for this problem. */
      virtual IloModel getModel() const = 0;
      /** Get the variables in this problem. */
      virtual IloNumVarArray getVariables() const = 0;
      /** Get the linear constraints in this problem. */
      virtual IloRangeArray getRows() const = 0;
      /** Get the block number for variable <code>x</code>.
       * @return The block number (in [0,getNBlocks()-1]) for <code>x</code>
       *         or -1 if <code>x</code> is in the master.
       */
      virtual IloInt getBlock(IloNumVar x) const = 0;
      /** Get the set of rows that intersect <code>x</code>. */
      virtual RowSet const &getIntersectedRows(IloNumVar x) const = 0;
      /** Get the objective function coefficient for <code>x</code>. */
      virtual double getObjCoef(IloNumVar x) const = 0;
      /** Get the objective function sense in this model. */
      virtual IloObjective::Sense getObjSense() const = 0;
   };

   static bool solve (Problem const *problem,
                      int argc, char const *const *argv,
                      std::vector<char const *> const &machines);

private:
   struct Block;
   struct LazyConstraintCallback;
   typedef std::vector<Block *> BlockVector;

   /** A block (either master or Benders) in a Benders decomposition.
    */
   struct Block {
      typedef std::map<IloNumVar,double,ExtractableLess<IloNumVar> > ObjMap;
      typedef std::map<IloNumVar,IloInt,ExtractableLess<IloNumVar> > IdxMap;
      typedef std::map<IloNumVar,IloNumVar,ExtractableLess<IloNumVar> >VarMap;


      IloEnv env;
      IloInt const number;  /**< Serial number of this block. */
      IloNumVarArray vars;  /**< Variables in this block's model. */
      IloRangeArray  rows;  /**< Rows in this block's model. */
      IloCplex       cplex; /**< The solver that (re)solves this block's model. */
      LazyConstraintCallback *cb;
      /** Description of variables that are fixed by master solves. */
      struct FixData {
         IloInt row;
         IloInt col;
         double val;
         FixData(IloInt r, IloInt c, double v) : row(r), col(c), val(v) {}
      };
      std::vector<Block::FixData> fixed;
      IloObjective obj;     /**< Objective function in this block. */
      ObjMap objMap;        /**< Objective coefficient for each variable. */
      VarMap varMap;        /**< Map variables in this block to variables
                             *   in original model. */

      // Extract a block from a problem.
      Block(Problem const *problem, IloInt n,
            int argc, char const *const *argv,
            std::vector<char const *> const &machines);
      // Extract master block from a problem.
      Block(Problem const *problem, BlockVector const &blocks);
      ~Block();
   };

   /** The callback that is used to separate Benders cuts at integer
    * feasible solutions.
    */
   struct LazyConstraintCallback : public IloCplex::LazyConstraintCallbackI {
      Block *const master;       /**< Master block. */
      BlockVector const &blocks; /**< Array of sub-blocks. */
      IloInt const etaind;       /**< Index of first eta variablein master. */
      IloCplex::AsyncHandle *handles;
      LazyConstraintCallback(IloEnv env, Block *m, BlockVector const &b,
                             IloInt e);
      ~LazyConstraintCallback();
      IloCplex::CallbackI *duplicateCallback() const;
      void main();
   };

   // Create dual of a linear program.
   static void makeDual(IloObjective const &primalObj,
                        IloNumVarArray const &primalVars,
                        IloRangeArray const &primalRows,
                        IloObjective *dualObj,
                        IloNumVarArray *dualVars,
                        IloRangeArray *dualRows);
};

/** Extract sub block number <code>n</code> from <code>problem</code>.
 * The constructor creates a representation of block number <code>n</code>
 * as described in <code>problem</code>.
 * The constructor will also connect the newly created block to a remote
 * object solver instance.
 * @param problem  The problem from which the block is to be extracted.
 * @param n        Index of the block to be extracted.
 * @param argc     Argument for IloCplex constructor.
 * @param argv     Argument for IloCplex constructor.
 * @param machines List of machines to which to connect. If the code is
 *                 compiled for the TCP/IP transport then the block will
 *                 be connected to <code>machines[n]</code>.
 */
BendersOpt::Block::Block(Problem const *problem, IloInt n, int argc, char const *const *argv,
                         std::vector<char const *> const &machines)
   : env(), number(n), vars(0), rows(0), cplex(0), cb(0)
{
   IloNumVarArray problemVars = problem->getVariables();
   IloRangeArray problemRanges = problem->getRows();

   // Create a map that maps variables in the original model to their
   // respective index in problemVars.
   std::map<IloNumVar,IloInt,ExtractableLess<IloNumVar> > origIdxMap;
   for (IloInt j = 0; j < problemVars.getSize(); ++j)
      origIdxMap.insert(std::map<IloNumVar,IloInt,ExtractableLess<IloNumVar> >::value_type(problemVars[j], j));

   // Copy non-fixed variables from original problem into primal problem.
   IloExpr primalObj(env);
   IloNumVarArray primalVars(env);
   IloRangeArray primalRows(env);
   IdxMap idxMap; // Index of original variable in block's primal model
   RowSet rowSet;
   for (IloInt j = 0; j < problemVars.getSize(); ++j) {
      IloNumVar x = problemVars[j];
      if ( problem->getBlock(x) == number ) {
         // Create column in block LP with exactly the same data.
         if ( x.getType() != IloNumVar::Float ) {
            std::stringstream s;
            s << "Cannot create non-continuous block variable " << x;
            std::cerr << s.str() << std::endl;
            throw s.str();
         }
         IloNumVar v(env, x.getLB(), x.getUB(), x.getType(), x.getName());
         // Normalize objective function to 'minimize'
         double coef = problem->getObjCoef(x);
         if ( problem->getObjSense() != IloObjective::Minimize )
            coef *= -1.0;
         primalObj += coef * v;
            
         // Record the index that the copied variable has in the
         // block model.
         idxMap.insert(IdxMap::value_type(x, primalVars.getSize()));
         primalVars.add(v);
            
         // Mark the rows that are intersected by this column
         // so that we can collect them later.
         RowSet const &intersected = problem->getIntersectedRows(x);
         for (RowSet::const_iterator it = intersected.begin();
              it != intersected.end(); ++it)
            rowSet.insert(*it);
      }
      else
         idxMap.insert(IdxMap::value_type(x, -1));
   }

   // Now copy all rows that intersect block variables.
   for (IloInt i = 0; i < problemRanges.getSize(); ++i) {
      IloRange r = problemRanges[i];
      if ( rowSet.find(r) == rowSet.end() )
         continue;

      // Create a copy of the row, normalizing it to '<='
      double factor = 1.0;
      if ( r.getLB() > -IloInfinity )
         factor = -1.0;
      IloRange primalR(env,
                       factor < 0 ? -r.getUB() : r.getLB(),
                       factor < 0 ? -r.getLB() : r.getUB(), r.getName());
      IloExpr lhs(env);
      for (IloExpr::LinearIterator it = r.getLinearIterator(); it.ok(); ++it)
      {
         IloNumVar v = it.getVar();
         double const val = factor * it.getCoef();
         if ( problem->getBlock(v) != number ) {
            // This column is not explicitly in this block. This means
            // that it is a column that will be fixed by the master.
            // We collect all such columns so that we can adjust the
            // dual objective function according to concrete fixings.
            // Store information about variables in this block that
            // will be fixed by master solves.
            fixed.push_back(FixData(primalRows.getSize(), origIdxMap[v], -val));
         }
         else {
            // The column is an ordinary in this block. Just copy it.
            lhs += primalVars[idxMap[v]] * val;
         }
      }
      primalR.setExpr(lhs);
      primalRows.add(primalR);
      lhs.end();
   }

   // Create the dual of the primal model we just created.
   // Note that makeDual _always_ returns a 'maximize' objective.
   IloObjective objective(env, primalObj, IloObjective::Minimize);
      
   makeDual(objective, primalVars, primalRows,
            &obj, &vars, &rows);
   objective.end();
   primalRows.endElements();
   primalRows.end();
   primalVars.endElements();
   primalVars.end();
   primalObj.end();
   // Create a model.
   IloModel model(env);
   model.add(obj);
   model.add(vars);
   model.add(rows);
   for (IloExpr::LinearIterator it = obj.getLinearIterator(); it.ok(); ++it)
      objMap.insert(ObjMap::value_type(it.getVar(), it.getCoef()));

   // Finally create the IloCplex instance that will solve
   // the problems associated with this block.
   char const **transargv = new char const *[argc + 3];
   for (int i = 0; i < argc; ++i)
      transargv[i] = argv[i];
#if defined(USE_MPI)
   char extra[128];
   sprintf (extra, "-remoterank=%d", static_cast<int>(number + 1));
   transargv[argc++] = extra;
   (void)machines;
#elif defined(USE_PROCESS)
   char extra[128];
   sprintf (extra, "-logfile=block%04d.log", static_cast<int>(number));
   transargv[argc++] = extra;
   (void)machines;
#elif defined(USE_TCPIP)
   transargv[argc++] = machines[number];
#endif
   cplex = IloCplex(model, TRANSPORT, argc, transargv);
   delete[] transargv;

   // Suppress output from this block's solver.
   cplex.setOut(env.getNullStream());
   cplex.setWarning(env.getNullStream());
}

/** Extract the master block from <code>problem</code>.
 * The constructor also sets up the solver for the newly created master
 * block. The master block can only be extracted if all sub-blocks have
 * already been extracted.
 * @param problem The problem from which to extract the master.
 * @param blocks  The sub blocks that have already been extracted.
 */
BendersOpt::Block::Block(Problem const *problem, BlockVector const &blocks)
   : env(), number(-1), vars(0), rows(0), cplex(0), cb(0)
{
   IloNumVarArray problemVars = problem->getVariables();
   IloRangeArray problemRanges = problem->getRows();

   IloExpr masterObj(env);
   IloNumVarArray masterVars(env);
   IloRangeArray masterRows(env);

   // Find columns that do not intersect block variables and
   // copy them to the master block.
   IdxMap idxMap;
   RowSet rowSet;
   for (IloInt j = 0; j < problemVars.getSize(); ++j) {
      IloNumVar x = problemVars[j];
      if ( problem->getBlock(x) < 0 ) {
         // Column is not in a block. Copy it to the master.
         IloNumVar v(env, x.getLB(), x.getUB(), x.getType(), x.getName());
         varMap.insert(VarMap::value_type(v, x));
         masterObj += problem->getObjCoef(x) * v;

         idxMap[x] = masterVars.getSize();
         masterVars.add(v);
      }
      else {
         // Column is in a block. Collect all rows that intersect
         // this column.
         RowSet const &intersected = problem->getIntersectedRows(x);
         for (RowSet::const_iterator it = intersected.begin();
              it != intersected.end(); ++it)
            rowSet.insert(*it);
         idxMap[x] = -1;
      }
   }

   // Pick up the rows that we need to copy.
   // These are the rows that are only intersected by master variables,
   // that is, the rows that are not in any block's rowset.
   for (IloInt i = 0; i < problemRanges.getSize(); ++i) {
      IloRange r = problemRanges[i];
      if ( rowSet.find(r) == rowSet.end() ) {
         IloRange masterRow(env, r.getLB(), r.getUB(), r.getName());
         IloExpr lhs(env);
         for (IloExpr::LinearIterator it = r.getLinearIterator(); it.ok(); ++it)
         {
            lhs += it.getCoef() * masterVars[idxMap[it.getVar()]];
         }
         masterRow.setExpr(lhs);
         masterRows.add(masterRow);
      }
   }

   // Adjust variable indices in blocks so that reference to variables
   // in the original problem become references to variables in the master.
   for (BlockVector::const_iterator b = blocks.begin(); b != blocks.end(); ++b) {
      for (std::vector<FixData>::iterator it = (*b)->fixed.begin(); it != (*b)->fixed.end(); ++it)
         it->col = idxMap[problemVars[it->col]];
   }

   // Create the eta variables, one for each block.
   // See the comments at the top of this file for details about the
   // eta variables.
   IloInt const firsteta = masterVars.getSize();
   for (BlockVector::size_type i = 0; i < blocks.size(); ++i) {
      std::stringstream s;
      s << "_eta" << i;
      IloNumVar eta(env, 0.0, IloInfinity, s.str().c_str());
      masterObj += eta;
      masterVars.add(eta);
   }

   // Create model and solver instance
   vars = masterVars;
   rows = masterRows;
   IloModel model(env);
   model.add(obj = IloObjective(env, masterObj, problem->getObjSense()));
   model.add(vars);
   model.add(rows);
   cplex = IloCplex(model);

   cplex.use(cb = new (env) LazyConstraintCallback(env, this, blocks,
                                              firsteta));

   for (IloExpr::LinearIterator it = obj.getLinearIterator(); it.ok(); ++it)
      objMap.insert(ObjMap::value_type(it.getVar(), it.getCoef()));
}

/** Destructor. */
BendersOpt::Block::~Block() {
   if ( cb ) cb->~LazyConstraintCallback();
   cplex.end();
   env.end();
}

/** Create the dual of a linear program.
 * The function can only dualize programs of the form
 * <code>Ax <= b, x >= 0</code>. The data in <code>primalVars</code> and
 * <code>dualRows</code> as well as in <code>primalRows</code> and
 * <code>dualVars</code> is in 1-to-1-correspondence.
 * @param primalObj  Objective function of primal problem.
 * @param primalVars Variables in primal problem.
 * @param primalRows Rows in primal problem.
 * @param dualObj    Objective function of dual will be stored here.
 * @param dualVars   All dual variables will be stored here.
 * @param dualRows   All dual rows will be stored here.
 */
void BendersOpt::makeDual(IloObjective const &primalObj,
                          IloNumVarArray const &primalVars,
                          IloRangeArray const &primalRows,
                          IloObjective *dualObj,
                          IloNumVarArray *dualVars,
                          IloRangeArray *dualRows)
{
   // To keep the code simple we only support problems
   // of the form Ax <= b, b >= 0 here. We leave it as a reader's
   // exercise to extend the function to something that can handle
   // any kind of linear model.
   for (IloInt j = 0; j < primalVars.getSize(); ++j)
      if ( primalVars[j].getLB() != 0 ||
           primalVars[j].getUB() < IloInfinity )
      {
         std::stringstream s;
         s << "Cannot dualize variable " << primalVars[j];
         throw s.str();
      }
   for (IloInt i = 0; i < primalRows.getSize(); ++i)
      if ( primalRows[i].getLB() > -IloInfinity ||
           primalRows[i].getUB() >= IloInfinity )
      {
         std::stringstream s;
         s << "Cannot dualize constraint " << primalRows[i];
         std::cerr << s.str() << std::endl;
         throw s.str();
      }

   // The dual of
   //   min/max c^T x
   //       Ax <= b
   //        x >= 0
   // is
   //   max/min y^T b
   //       y^T A <= c
   //           y <= 0
   // We scale y by -1 to get >= 0

   IloEnv env = primalVars.getEnv();
   IloObjective obj(env, 0.0,
                    primalObj.getSense() == IloObjective::Minimize ?
                    IloObjective::Maximize : IloObjective::Minimize);
   IloRangeArray rows(env);
   IloNumVarArray y(env);
   std::map<IloNumVar,IloInt,ExtractableLess<IloNumVar> > v2i;
   for (IloInt j = 0; j < primalVars.getSize(); ++j) {
      IloNumVar x = primalVars[j];
      v2i.insert(std::map<IloNumVar,IloInt,ExtractableLess<IloNumVar> >::value_type(x, j));
      rows.add(IloRange(env, -IloInfinity, 0, x.getName()));
   }
   for (IloExpr::LinearIterator it = primalObj.getLinearIterator(); it.ok(); ++it)
      rows[v2i[it.getVar()]].setUB(it.getCoef());

   for (IloInt i = 0; i < primalRows.getSize(); ++i) {
      IloRange r = primalRows[i];
      IloNumColumn col(env);
      col += obj(-r.getUB());
      for (IloExpr::LinearIterator it = r.getLinearIterator(); it.ok(); ++it)
         col += rows[v2i[it.getVar()]](-it.getCoef());
      y.add(IloNumVar(col, 0, IloInfinity, IloNumVar::Float, r.getName()));
   }

   *dualObj = obj;
   *dualVars = y;
   *dualRows = rows;
}

/** Create callback.
 * The newly created callback separates Benders cuts for the problem
 * described by <code>m</code> (the master) and <code>b</code> (the
 * sub blocks).
 * @param env Environment in which the callback is created.
 * @param m   Master block for Benders decomposition.
 * @param b   Sub blocks for Benders decomposition.
 * @param e   Index of first eta variable in master problem.
 */
BendersOpt::LazyConstraintCallback::LazyConstraintCallback(IloEnv env,
                                                           Block *m,
                                                           BlockVector const &b,
                                                           IloInt e)
   : IloCplex::LazyConstraintCallbackI(env),
     master(m), blocks(b), etaind(e),
     handles(new IloCplex::AsyncHandle [blocks.size()])
{
}

/** Destructor. */
BendersOpt::LazyConstraintCallback::~LazyConstraintCallback() {
   delete[] handles;
}

/** Clone function as required by IloCplex::CallbackI. */
IloCplex::CallbackI *BendersOpt::LazyConstraintCallback::duplicateCallback() const {
   return new (getEnv()) LazyConstraintCallback(getEnv(), master,
                                                blocks, etaind);
}

/** Separation function.
 * This function is invoked whenever CPLEX finds an integer feasible
 * solution. It then separates either feasibility or optimality cuts
 * on this solution.
 */
void BendersOpt::LazyConstraintCallback::main() {
   std::cout << "Callback invoked. Separate Benders cuts." << std::endl;

   IloNumArray x(getEnv());
   IloNumArray rayVals(getEnv());
   IloNumVarArray rayVars(getEnv());
   IloNumArray cutVal(getEnv());
   IloNumVarArray cutVar(getEnv());

   getValues(x, master->vars);

   bool error = false;

   // Iterate over blocks and trigger a separation on each of them.
   // The separation is triggered asynchronously so that it can happen
   // on different remote objects simultaneously.
   for (BlockVector::size_type b = 0; b < blocks.size(); ++b) {
      Block *const block = blocks[b];

      // Remove current objective from the block's model.
      IloModel model = block->cplex.getModel();
      IloObjective obj = block->obj;
      model.remove(obj);
      IloExpr newObj = obj.getExpr();
         
      // Iterate over the fixed master variables in this block to update
      // the block's objective function.
      // Each fixed variable goes to the right-hand side and therefore
      // into the objective function.
      for (std::vector<Block::FixData>::const_iterator it = block->fixed.begin(); it != block->fixed.end(); ++it)
         newObj -= block->vars[it->row] * (it->val * x[it->col]);
      obj.setExpr(newObj);
      model.add(obj);
      newObj.end();

      // If the problem is unbounded we need to get an infinite ray in
      // order to be able to generate the respective Benders cut. If
      // CPLEX proves unboundedness in presolve then it will return
      // CPX_STAT_INForUNBD and no ray will be available. So we need to
      // disable presolve.
      block->cplex.setParam(IloCplex::Param::Preprocessing::Presolve, false);
      block->cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);

      // Solve the updated problem to optimality.
      block->cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
      try {
         handles[b] = block->cplex.solve(true);
      } catch (...) {
         // If there is an exception then we need to kill and join
         // all remaining solves. Otherwise we may leak handles.
         while (--b > 0) {
            handles[b].kill();
            handles[b].join();
         }
         throw;
      }
   }

   // Wait for the various LP solves to complete.
   for (BlockVector::size_type b = 0; b < blocks.size(); ++b)
      handles[b].join();

   // See if we need to generate cuts.
   for (BlockVector::size_type b = 0; b < blocks.size(); ++b) {
      Block *const block = blocks[b];
      cutVal.clear();
      cutVar.clear();
      double cutlb = -IloInfinity;
      double cutub = IloInfinity;

      // We ust STL types here since they are exception safe.
      std::vector<double> tmp(master->vars.getSize(), 0);
      std::map<IloNumVar,double,ExtractableLess<IloNumVar> > rayMap;

      // Depending on the status either seperate a feasibility or an
      // optimality cut.
      switch (block->cplex.getStatus()) {
      case IloAlgorithm::Unbounded:
         {
            // The subproblem is unbounded. We need to extract a feasibility
            // cut from an unbounded ray of the problem (see also the comments
            // at the top of this file).
            std::cout << "Block " << b << " unbounded ";
            block->cplex.getRay(rayVals, rayVars);
            cutub = 0.0;
            for (IloInt j = 0; j < rayVars.getSize(); ++j) {
               cutub -= rayVals[j] * block->objMap[rayVars[j]];
               rayMap[rayVars[j]] = rayVals[j];
            }
            for (std::vector<Block::FixData>::const_iterator it = block->fixed.begin(); it != block->fixed.end(); ++it)
               tmp[it->col] -= it->val * rayMap[block->vars[it->row]];
            for (IloInt j = 0; j < master->vars.getSize(); ++j) {
               if ( fabs (tmp[j]) > EPSILON ) {
                  cutVar.add(master->vars[j]);
                  cutVal.add(tmp[j]);
               }
            }
         }
         break;
      case IloAlgorithm::Optimal:
         {
            // The subproblem has a finite optimal solution.
            // We need to check if this gives rise to an optimality cut (see
            // also the comments at the top of this file).            
            std::cout << "Block " << b << " optimal ";
            double const objval = block->cplex.getObjValue();
            double const eta = x[etaind + b];
            block->cplex.getValues(block->vars, rayVals);
            
            if ( objval > eta + EPSILON ) {
               cutub = 0.0;
               for (IloInt j = 0; j < block->vars.getSize(); ++j)
                  cutub -= rayVals[j] * block->objMap[block->vars[j]];
               for (std::vector<Block::FixData>::const_iterator it = block->fixed.begin(); it != block->fixed.end(); ++it)
                  tmp[it->col] -= it->val * rayVals[it->row];
               for (IloInt j = 0; j < master->vars.getSize(); ++j) {
                  if ( fabs (tmp[j]) > EPSILON ) {
                     cutVal.add(tmp[j]);
                     cutVar.add(master->vars[j]);
                  }
               }
               cutVal.add(-1.0);
               cutVar.add(master->vars[etaind + b]);
            }
         }
         break;
      default:
         std::cerr << "Block " << b << " Unexpected status "
                   << block->cplex.getStatus() << std::endl;
         error = true;
         break;
      }
      
      // If a cut was found then add that.
      if ( cutVar.getSize() > 0 ) {
         IloExpr expr(master->env);
         for (IloInt i = 0; i < cutVar.getSize(); ++i)
            expr += cutVar[i] * cutVal[i];
         IloRange cut(getEnv(), cutlb, expr, cutub);
         expr.end();
         std::cout << "cut found: " << cut << std::endl;
         add(cut).end();
      }
      else
         std::cout << "no cuts." << std::endl;
   }
   cutVar.end();
   cutVal.end();
   rayVars.end();
   rayVals.end();
   x.end();
   if ( error )
      throw -1;
}

/** Solve the <code>problem</code> using a distributed implementation of
 * Benders' decomposition.
 */
bool
BendersOpt::solve (Problem const *problem,
                   int argc, char const *const *argv,
                   std::vector<char const *> const &machines)
{
   IloEnv env = problem->getModel().getEnv();

   std::vector<Block *> blocks;
   Block *master = 0;
   bool result = false;
   try {
      // Extract blocks and master problem.
      std::cout << "Extracting " << problem->getNBlocks() << " blocks."
                << std::endl;
      for (IloInt b = 0; b < problem->getNBlocks(); ++b)
         blocks.push_back(new Block(problem, b, argc, argv, machines));
      master = new Block(problem, blocks);

      // Write out the master and all blocks (for debugging).
      master->cplex.exportModel("master.lp");
      for (BlockVector::size_type b = 0; b < blocks.size(); ++b) {
         std::stringstream s;
         s << "block" << b << ".lp";
         blocks[b]->cplex.exportModel(s.str().c_str());
      }

      // Solve the master.
      // If we find a feasible solution then perform a final solve
      // with the original problem so that we get solution values for
      // all variables in the original problem.
      if ( master->cplex.solve() ) {
         // Perform a last solve to get the solution values.
         IloNumArray vals(env), startVals(env);
         IloNumVarArray startVars(env);
         master->cplex.getValues(vals, master->vars);
         IloCplex cplex(problem->getModel());
         // Fix integral variables to their value in the master solution
         // and add them as MIP start.
         for (IloInt i = 0; i < vals.getSize(); ++i)
            if ( master->vars[i].getType() != IloNumVar::Float ) {
               double const v = IloRound(vals[i]);
               IloNumVar x = master->varMap[master->vars[i]];
               startVals.add(v);
               startVars.add(x);
               // We add lazy constraints so as to make sure that
               // - we don't modify the original model
               // - the problem has only the unique optimal solution
               //   we are interested in
               cplex.addLazyConstraint(x == v);
            }
         cplex.addMIPStart(startVars, startVals);
         cplex.solve();

         // Report the results.
         std::cout << "#### Problem solved (" << cplex.getObjValue() << ", "
                   << cplex.getCplexStatus() << ")." << std::endl;
         IloNumVarArray vars = problem->getVariables();
         for (IloInt i = 0; i < vars.getSize(); ++i)
            std::cout << "#### \tx[" << i << "] = " << cplex.getValue(vars[i])
                      << std::endl;
         cplex.end();
         result = true;
      }
   } catch (...) {
      if ( master )
         delete master;
      while ( blocks.size() > 0 ) {
         delete blocks.back();
         blocks.pop_back();
      }
      throw;
   }

   delete master;
   while ( blocks.size() > 0 ) {
      delete blocks.back();
      blocks.pop_back();
   }

   return result;
}

// ----------------------------------------------------------------------

/** Create an example problem.
 * The function creates the following facility location problem:
 *
 *   minimize sum(f in F) sum(j in J) C[f,j]*x[f,j] + sum(f in F) F[f]*y[f]
 *     such that    sum(f in F) x[f,j] >= 1         forall j in J
 *                   sum (f in F) y[f] <= 2
 *                                y[f] >= x[f,j]    forall f in F, j in J
 *                              x[f,j] >= 0         forall f in F, j in J
 *                                y[f] in { 0, 1 }  forall f in F
 * where y[f] = 1 if factory f is opened and 0 otherwise and
 * x[f,j] is the fraction of demand serviced from factory f to customer j.
 * F[f] is the fixed charge cost for opening factory f and C[f,j] is the
 * cost for servicing customer j from f.
 *
 * The output of the function is an CPXLPptr instance and a description
 * of the blocks in this problem. The description of blocks is as follows:
 *   *NBLOCKS_P specifies the number of blocks (not accounting for the
 *              master block).
 *   *BLOCK_P   is an array that has one entry for each column in the
 *              returned CPXLPptr object. A negative entry specifies that
 *              the column appears only in the master. A non-negative entry
 *              specifies that the variable is in a block and gives the
 *              block index.
 */
class Example : public BendersOpt::Problem {
   typedef std::map<IloNumVar,IloInt,ExtractableLess<IloNumVar> > BlockMap;
   typedef std::map<IloNumVar,RowSet,ExtractableLess<IloNumVar> > IntersectMap;
   typedef std::map<IloNumVar,double,ExtractableLess<IloNumVar> > ObjMap;

   IloInt nblocks;
   IloModel model;
   IloNumVarArray vars;
   IloRangeArray ranges;
   BlockMap blockMap;
   IntersectMap intersectMap;
   ObjMap objMap;
   IloObjective::Sense objSense;

public:
   Example(IloEnv env)
      : nblocks(0), model(env), vars(env), ranges(env)
   {
      // Model data.
      // fixed[] is the fixed cost for opening a facility,
      // cost[i,j] is the cost for serving customer i from facility j.
      static double const fixed[] = { 2.0, 3.0, 3.0 };
      static double const cost[] = { 2.0, 3.0, 4.0, 5.0, 7.0,
                                     4.0, 3.0, 1.0, 2.0, 6.0,
                                     5.0, 4.0, 2.0, 1.0, 3.0 };
#define NFACTORY ((CPXDIM)(sizeof(fixed) / sizeof(fixed[0])))
#define NCUSTOMER ((CPXDIM)((sizeof(cost) / sizeof(cost[0])) / NFACTORY))
      nblocks = NCUSTOMER;

      IloExpr obj(env);
      // Create integer y  variables.
      IloNumVarArray y(env);
      for (IloInt f = 0; f < NFACTORY; ++f) {
         std::stringstream s;
         s << "y" << f;
         IloIntVar v(env, 0, 1, s.str().c_str());
         obj += fixed[f] * v;
         objMap[v] = fixed[f];
         y.add(v);
         blockMap.insert(BlockMap::value_type(v, -1));
         intersectMap.insert(IntersectMap::value_type(v, RowSet()));
      }

      // Create continuous x variables.
      IloNumVarArray x(env);
      for (IloInt f = 0; f < NFACTORY; ++f) {
         for (IloInt c = 0; c < NCUSTOMER; ++c) {
            std::stringstream s;
            s << "x" << f << "#" << c;
            IloNumVar v(env, 0.0, IloInfinity, s.str().c_str());
            obj += v * cost[f * NCUSTOMER + c];
            objMap[v] = cost[f * NCUSTOMER + c];
            x.add(v);
            blockMap.insert(BlockMap::value_type(v, c));
            intersectMap.insert(IntersectMap::value_type(v, RowSet()));
         }
      }
      vars.add(y);
      vars.add(x);
      model.add(vars);

      // Add objective function.
      model.add(IloMinimize(env, obj, "obj"));
      objSense = IloObjective::Minimize;
      obj.end();

      // Satisfy each customer's demand.
      for (IloInt c = 0; c < NCUSTOMER; ++c) {
         std::stringstream s;
         s << "c1_" << c;
         IloRange r(env, 1.0, IloInfinity, s.str().c_str());
         IloExpr lhs(env);
         for (IloInt f = 0; f < NFACTORY; ++f) {
            lhs += x[f * NCUSTOMER + c];
            intersectMap[x[f * NCUSTOMER + c]].insert(r);
         }
         r.setExpr(lhs);
         ranges.add(r);
         lhs.end();
      }

      // A factory must be open if we service from it.
      for (IloInt c = 0; c < NCUSTOMER; ++c) {
         for (IloInt f = 0; f < NFACTORY; ++f) {
            std::stringstream s;
            s << "c2_" << c << "#" << f;
            IloRange r(env, 0.0, IloInfinity, s.str().c_str());
            intersectMap[x[f * NCUSTOMER + c]].insert(r);
            intersectMap[y[f]].insert(r);
            r.setExpr(-x[f * NCUSTOMER + c] + y[f]);
            ranges.add(r);
         }
      }

      // Capacity constraint.
      IloRange r(env, -IloInfinity, NFACTORY - 1, "c3");
      IloExpr lhs(env);
      for (IloInt f = 0; f < NFACTORY; ++f) {
         lhs += y[f];
         intersectMap[y[f]].insert(r);
      }
      r.setExpr(lhs);
      ranges.add(r);
      lhs.end();

      model.add(ranges);

#undef NFACTORY
#undef NCUSTOMER
   }

   // Implementation of functions required by the BendersOpt::Problem
   // interface.
   IloInt getNBlocks() const { return nblocks; }
   IloModel getModel() const { return model; }
   IloNumVarArray getVariables() const { return vars; }
   IloRangeArray getRows() const { return ranges; }
   IloInt getBlock(IloNumVar x) const {
      BlockMap::const_iterator const it = blockMap.find(x);
      return (it == blockMap.end()) ? -1 : it->second;
   }
   RowSet const &getIntersectedRows(IloNumVar x) const {
      static RowSet const empty;
      IntersectMap::const_iterator const it = intersectMap.find(x);
      return (it == intersectMap.end()) ? empty : it->second;
   }
   double getObjCoef(IloNumVar x) const {
      ObjMap::const_iterator const it = objMap.find(x);
      return (it == objMap.end()) ? 0.0 : it->second;
   }
   IloObjective::Sense getObjSense() const { return objSense; }
};

   // ----------------------------------------------------------------------

int
main (int argc, char **argv)
{
   int myargc;
   char const **myargv = NULL;
   std::vector<char const *> machines;
   int nmachines = 0;

#if defined(USE_MPI)
   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &nmachines);
   --nmachines;
   myargc = 0;
#elif defined(USE_PROCESS)
   nmachines = INT_MAX;
   myargc = 2;
#elif defined(USE_TCPIP)
   nmachines = 0;
   myargc = 0;
#endif


   // Process the command line.
   myargv = new char const *[argc];

   for (int i = 1; i < argc; ++i) {
#if defined(USE_MPI)
      /* Nothing to do. */
      if ( 0 )
         ;
#elif defined(USE_PROCESS)
      if ( strncmp (argv[i], "-bin=", 5) == 0 ) {
         myargv[0] = argv[i] + 5;
         myargv[1] = "-worker=process";
      }
#elif defined(USE_TCPIP)
      if ( strncmp (argv[i], "-address=", 9) == 0 ) {
         machines.push_back (argv[i]);
         ++nmachines;
      }
#endif
      else
         myargv[myargc++] = argv[i];
   }

   IloEnv env;

   try {
      Example problem(env);

      if ( nmachines < problem.getNBlocks() ) {
         std::cerr << "We have " << problem.getNBlocks()
                   << " blocks blocks but only "
                   << nmachines << " machines!" << std::endl;
         throw -1;
      }

      BendersOpt::solve (&problem, myargc, myargv, machines);
   } catch (...) {
      env.end();
      delete[] myargv;
      throw;
   }
   env.end();

   delete[] myargv;

#if defined(USE_MPI)
   MPI_Finalize ();
#endif


   return 0;
}

#endif /* COMPILE_MASTER */
