/* --------------------------------------------------------------------------
 * File: iloparbenders.java
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
import ilog.cplex.*;
import ilog.concert.*;

import java.util.HashMap;
import java.util.Vector;

/* ********************************************************************** *
 *                                                                        *
 *    M a s t e r   s i d e   i m p l e m e n t a t i o n                 *
 *                                                                        *
 * ********************************************************************** */

/** A class that solves problems by means of parallel distributed Benders
 * decomposition. The public API of this class consists of only two things:
 * - An interface Problem that describes problems that
 *   can be solved by this class.
 * - A function solve() that solves instances of Problem.
 */
public final class RemoteParBenders {

   // ----------------------------------------------------------------------

   private static final class RowSet extends java.util.HashSet<IloRange> {}

   /** Interface to problems that can be solved by the {@link RemoteParBenders} class. */
   public interface Problem {
      /** Get the number of blocks in this problem. */
      abstract int getNBlocks();
      /** Get the model for this problem. */
      abstract IloCplex getModel();
      /** Get the variables in this problem. */
      abstract IloNumVar[] getVariables();
      /** Get the linear constraints in this problem. */
      abstract IloRange[] getRows();
      /** Get the block number for variable <code>x</code>.
       * @return The block number (in [0,getNBlocks()-1]) for <code>x</code>
       *         or -1 if <code>x</code> is in the master.
       */
      abstract int getBlock(IloNumVar x);
      /** Get the set of rows that intersect <code>x</code>. */
      abstract RowSet getIntersectedRows(IloNumVar x);
      /** Get the objective function coefficient for <code>x</code>. */
      abstract double getObjCoef(IloNumVar x);
      /** Get the objective function sense in this model. */
      abstract IloObjectiveSense getObjSense();
   }

   /** A block (either master or Benders) in a Benders decomposition.
    */
   private static final class Block {
      final int         number; /**< Serial number of this block. */
      final IloNumVar[] vars;   /**< Variables in this block's model. */
      final IloRange[]  rows;   /**< Rows in this block's model. */
      final IloCplex    cplex;  /**< The solver that (re)solves this block's model. */
      /** Description of variables that are fixed by master solves. */
      public static final class FixData {
         public int row;
         public int col;
         public double val;
         public FixData(int row, int col, double val) {
            this.row = row;
            this.col = col;
            this.val = val;
         }
      }
      final Vector<FixData> fixed;
      final IloObjective obj;                    /**< Objective function in this block. */
      final HashMap<IloNumVar,Double> objMap;    /**< Objective coefficient for each variable. */
      final HashMap<IloNumVar,IloNumVar> varMap; /**< Map variables in this block to variables
                                                   *   in original model. */


      /** Extract sub block number <code>n</code> from <code>problem</code>.
       * The constructor creates a representation of block number <code>n</code>
       * as described in <code>problem</code>.
       * The constructor will also connect the newly created block to a remote
       * object solver instance.
       * @param problem   The problem from which the block is to be extracted.
       * @param n         Index of the block to be extracted.
       * @param transport Argument for IloCplex constructor.
       * @param args      Argument for IloCplex constructor.
       * @param argv      Argument for IloCplex constructor.
       * @param machines  List of machines to which to connect. If the code is
       *                  compiled for the TCP/IP transport then the block will
       *                  be connected to <code>machines[n]</code>.
       */
      public Block(Problem problem, int n,
                   String transport, String[] args,
                   String[] machines) throws IloException
      {
         this.number = n;
         this.fixed = new Vector<FixData>();
         this.objMap = new HashMap<IloNumVar,Double>();
         this.varMap = new HashMap<IloNumVar,IloNumVar>();
         final IloNumVar[] problemVars = problem.getVariables();
         final IloRange[] problemRanges = problem.getRows();

         // Create a map that maps variables in the original model to their
         // respective index in problemVars.
         HashMap<IloNumVar,Integer> origIdxMap = new HashMap<IloNumVar,Integer>();
         for (int j = 0; j < problemVars.length; ++j)
            origIdxMap.put(problemVars[j], j);

         // Copy non-fixed variables from original problem into primal problem.
         final IloCplexModeler primal = new IloCplexModeler();
         final IloLinearNumExpr primalObj = primal.linearNumExpr();
         final Vector<IloNumVar> primalVars = new Vector<IloNumVar>();
         final Vector<IloRange> primalRows = new Vector<IloRange>();
         final HashMap<IloNumVar,Integer> idxMap = new HashMap<IloNumVar,Integer>(); // Index of original variable in block's primal model
         final RowSet rowSet = new RowSet();
         for (int j = 0; j < problemVars.length; ++j) {
            IloNumVar x = problemVars[j];
            if ( problem.getBlock(x) == number ) {
               // Create column in block LP with exactly the same data.
               if ( x.getType() != IloNumVarType.Float )
                  throw new RuntimeException("Cannot create non-continuous block variable " + x);
               IloNumVar v = primal.numVar(x.getLB(), x.getUB(), x.getType(), x.getName());
               // Normalize objective function to 'minimize'
               double coef = problem.getObjCoef(x);
               if ( problem.getObjSense() != IloObjectiveSense.Minimize )
                  coef *= -1.0;
               primalObj.addTerm(coef, v);
            
               // Record the index that the copied variable has in the
               // block model.
               idxMap.put(x, primalVars.size());
               primalVars.add(v);
            
               // Mark the rows that are intersected by this column
               // so that we can collect them later.
               RowSet intersected = problem.getIntersectedRows(x);
               for (IloRange r : problem.getIntersectedRows(x))
                  rowSet.add(r);
            }
            else
               idxMap.put(x, -1);
         }

         // Now copy all rows that intersect block variables.
         for (int i = 0; i < problemRanges.length; ++i) {
            IloRange r = problemRanges[i];
            if ( !rowSet.contains(r) )
               continue;

            // Create a copy of the row, normalizing it to '<='
            double factor = 1.0;
            if ( r.getLB() > Double.NEGATIVE_INFINITY )
               factor = -1.0;
            IloRange primalR = primal.range(factor < 0 ? -r.getUB() : r.getLB(),
                                            factor < 0 ? -r.getLB() : r.getUB(),
                                            r.getName());
            IloLinearNumExpr lhs = primal.linearNumExpr();
            for (IloLinearNumExprIterator it = ((IloLinearNumExpr)r.getExpr()).linearIterator(); it.hasNext(); /* nothing */) {
               IloNumVar v = it.nextNumVar();
               double val = factor * it.getValue();
               if ( problem.getBlock(v) != number ) {
                  // This column is not explicitly in this block. This means
                  // that it is a column that will be fixed by the master.
                  // We collect all such columns so that we can adjust the
                  // dual objective function according to concrete fixings.
                  // Store information about variables in this block that
                  // will be fixed by master solves.
                  fixed.add(new FixData(primalRows.size(), origIdxMap.get(v), -val));
               }
               else {
                  // The column is an ordinary in this block. Just copy it.
                  lhs.addTerm(val, primalVars.elementAt(idxMap.get(v)));
               }
            }
            primalR.setExpr(lhs);
            primalRows.add(primalR);
         }

         // Create the IloCplex instance that will solve
         // the problems associated with this block.
         final Vector<String> transargv = new Vector<String>();
         transargv.addAll(java.util.Arrays.asList(args));
         if ( transport.equals("processtransport") ) {
            transargv.add("-logfile=block" + number + ".log");
         }
         else if ( transport.equals("tcpiptransport") ) {
            transargv.add(machines[number]);
         }
         cplex = new IloCplex(transport, transargv.toArray(new String[0]));
         try {
            // Suppress output from this block's solver.
            cplex.setOut(null);
            cplex.setWarning(null);

            // Create the dual of the primal model we just created.
            // Note that makeDual _always_ returns a 'maximize' objective.
            IloObjective objective = primal.minimize(primalObj);
            Vector<IloNumVar> dualVars = new Vector<IloNumVar>();
            Vector<IloRange> dualRows = new Vector<IloRange>();
      
            makeDual(objective,
                     primalVars.toArray(new IloNumVar[0]),
                     primalRows.toArray(new IloRange[0]),
                     cplex, dualVars, dualRows);
            vars = dualVars.toArray(new IloNumVar[0]);
            rows = dualRows.toArray(new IloRange[0]);
            obj = cplex.getObjective();

            for (IloLinearNumExprIterator it = ((IloLinearNumExpr)obj.getExpr()).linearIterator(); it.hasNext(); /* nothing */) {
               IloNumVar v = it.nextNumVar();
               objMap.put(v, it.getValue());
            }
         } catch (IloException e) {
            cplex.end();
            throw e;
         }
      }

      /** Extract the master block from <code>problem</code>.
       * The constructor also sets up the solver for the newly created master
       * block. The master block can only be extracted if all sub-blocks have
       * already been extracted.
       * @param problem The problem from which to extract the master.
       * @param blocks  The sub blocks that have already been extracted.
       */
      Block(Problem problem, Vector<Block> blocks) throws IloException {
         this.number = -1;
         this.fixed = null;
         this.objMap = new HashMap<IloNumVar,Double>();
         this.varMap = new HashMap<IloNumVar,IloNumVar>();
         this.cplex = new IloCplex();
         try {
            final IloNumVar[] problemVars = problem.getVariables();
            final IloRange[] problemRanges = problem.getRows();

            IloLinearNumExpr masterObj = cplex.linearNumExpr();
            Vector<IloNumVar> masterVars = new Vector<IloNumVar>();
            Vector<IloRange> masterRows = new Vector<IloRange>();

            // Find columns that do not intersect block variables and
            // copy them to the master block.
            HashMap<IloNumVar,Integer> idxMap = new HashMap<IloNumVar,Integer>();
            RowSet rowSet = new RowSet();
            for (int j = 0; j < problemVars.length; ++j) {
               IloNumVar x = problemVars[j];
               if ( problem.getBlock(x) < 0 ) {
                  // Column is not in a block. Copy it to the master.
                  IloNumVar v = cplex.numVar(x.getLB(), x.getUB(), x.getType(), x.getName());
                  varMap.put(v, x);
                  masterObj.addTerm(problem.getObjCoef(x), v);

                  idxMap.put(x, masterVars.size());
                  masterVars.add(v);
               }
               else {
                  // Column is in a block. Collect all rows that intersect
                  // this column.
                  for (IloRange r : problem.getIntersectedRows(x))
                     rowSet.add(r);
                  idxMap.put(x, -1);
               }
            }

            // Pick up the rows that we need to copy.
            // These are the rows that are only intersected by master variables,
            // that is, the rows that are not in any block's rowset.
            for (int i = 0; i < problemRanges.length; ++i) {
               IloRange r = problemRanges[i];
               if ( !rowSet.contains(r) ) {
                  IloRange masterRow = cplex.range(r.getLB(), r.getUB(), r.getName());
                  IloLinearNumExpr lhs = cplex.linearNumExpr();
                  for (IloLinearNumExprIterator it = ((IloLinearNumExpr)r.getExpr()).linearIterator(); it.hasNext(); /* nothing */) {
                     final IloNumVar var = it.nextNumVar();
                     lhs.addTerm(it.getValue(), masterVars.elementAt(idxMap.get(var)));
                  }
                  masterRow.setExpr(lhs);
                  masterRows.add(masterRow);
               }
            }

            // Adjust variable indices in blocks so that reference to variables
            // in the original problem become references to variables in the master.
            for (Block b : blocks) {
               for (FixData f : b.fixed)
                  f.col = idxMap.get(problemVars[f.col]);
            }

            // Create the eta variables, one for each block.
            // See the comments at the top of this file for details about the
            // eta variables.
            final int firsteta = masterVars.size();
            for (int i = 0; i < blocks.size(); ++i) {
               IloNumVar eta = cplex.numVar(0.0, Double.POSITIVE_INFINITY, "_eta" + i);
               masterObj.addTerm(1.0, eta);
               masterVars.add(eta);
            }

            // Create model and solver instance
            vars = masterVars.toArray(new IloNumVar[0]);
            rows = masterRows.toArray(new IloRange[0]);
            obj = cplex.addObjective(problem.getObjSense(), masterObj);
            cplex.add(vars);
            cplex.add(rows);

            cplex.use(new LazyConstraintCallback(this, blocks, firsteta));

            for (IloLinearNumExprIterator it = masterObj.linearIterator(); it.hasNext(); /* nothing */) {
               IloNumVar v = it.nextNumVar();
               objMap.put(v, it.getValue());
            }
         } catch (IloException e) {
            cplex.end();
            throw e;
         }
      }

      public void end() { cplex.end(); }
      public void finalize() { end(); }

   }

   /** The callback that is used to separate Benders cuts at integer
    * feasible solutions.
    */
   private static final class LazyConstraintCallback extends IloCplex.LazyConstraintCallback {
      final Block master;                   /**< Master block. */
      final Vector<Block> blocks; /**< Array of sub-blocks. */
      final int etaind;                     /**< Index of first eta variablein master. */
      final IloCplex.SolveHandle[] handles;
      /** Create callback.
       * The newly created callback separates Benders cuts for the problem
       * described by <code>m</code> (the master) and <code>b</code> (the
       * sub blocks).
       * @param env Environment in which the callback is created.
       * @param m   Master block for Benders decomposition.
       * @param b   Sub blocks for Benders decomposition.
       * @param e   Index of first eta variable in master problem.
       */
      public LazyConstraintCallback(Block m, Vector<Block> b, int e) {
         this.master = m;
         this.blocks = b;
         this.etaind = e;
         this.handles = new IloCplex.SolveHandle[b.size()];
      }

      /** Separation function.
       * This function is invoked whenever CPLEX finds an integer feasible
       * solution. It then separates either feasibility or optimality cuts
       * on this solution.
       */
      protected void main() throws IloException {
         System.out.println("Callback invoked. Separate Benders cuts.");

         double[] cutVal = new double[master.vars.length];
         IloNumVar[] cutVar = new IloNumVar[master.vars.length];
         double[] x = getValues(master.vars);

         boolean error = false;

         // Iterate over blocks and trigger a separation on each of them.
         // The separation is triggered asynchronously so that it can happen
         // on different remote objects simultaneously.
         for (int b = 0; b < blocks.size(); ++b) {
            Block block = blocks.elementAt(b);

            // Remove current objective from the block's model.
            IloObjective obj = block.obj;
            block.cplex.remove(obj);
            IloLinearNumExpr newObj = (IloLinearNumExpr)obj.getExpr();
         
            // Iterate over the fixed master variables in this block to update
            // the block's objective function.
            // Each fixed variable goes to the right-hand side and therefore
            // into the objective function.
            for (Block.FixData it : block.fixed)
               newObj.addTerm(-(it.val * x[it.col]), block.vars[it.row]);
            obj.setExpr(newObj);
            block.cplex.add(obj);

            // If the problem is unbounded we need to get an infinite ray in
            // order to be able to generate the respective Benders cut. If
            // CPLEX proves unboundedness in presolve then it will return
            // CPX_STAT_INForUNBD and no ray will be available. So we need to
            // disable presolve.
            block.cplex.setParam(IloCplex.Param.Preprocessing.Presolve, false);
            block.cplex.setParam(IloCplex.Param.Preprocessing.Reduce, 0);

            // Solve the updated problem to optimality.
            block.cplex.setParam(IloCplex.Param.RootAlgorithm, IloCplex.Algorithm.Primal);
            try {
               handles[b] = block.cplex.solve(true);
            } catch (IloException e) {
               // If there is an exception then we need to kill and join
               // all remaining solves. Otherwise we may leak handles.
               while (--b > 0) {
                  handles[b].kill();
                  handles[b].join();
               }
               throw e;
            }
         }

         // Wait for the various LP solves to complete.
         for (int b = 0; b < blocks.size(); ++b)
            handles[b].join();

         // See if we need to generate cuts.
         for (int b = 0; b < blocks.size(); ++b) {
            Block block = blocks.elementAt(b);
            double cutlb = Double.NEGATIVE_INFINITY;
            double cutub = Double.POSITIVE_INFINITY;
            int cutNz = 0;

            // We ust STL types here since they are exception safe.
            double[] tmp = new double[master.vars.length];
            HashMap<IloNumVar,Double> rayMap = new HashMap<IloNumVar,Double>();

            // Depending on the status either seperate a feasibility or an
            // optimality cut.
            if (block.cplex.getStatus().equals(IloCplex.Status.Unbounded)) {
               // The subproblem is unbounded. We need to extract a feasibility
               // cut from an unbounded ray of the problem (see also the comments
               // at the top of this file).
               System.out.print("unbounded ");
               cutub = 0.0;
               for (IloLinearNumExprIterator it = block.cplex.getRay().linearIterator(); it.hasNext(); /* nothing */) {
                  IloNumVar rayVar = it.nextNumVar();
                  double rayVal = it.getValue();
                  cutub -= rayVal * block.objMap.get(rayVar);
                  rayMap.put(rayVar, rayVal);
               }
               for (Block.FixData f : block.fixed)
                  tmp[f.col] -= f.val * rayMap.get(block.vars[f.row]);
               for (int j = 0; j < master.vars.length; ++j) {
                  if ( Math.abs(tmp[j]) > 1e-6 ) {
                     cutVar[cutNz] = master.vars[j];
                     cutVal[cutNz] = tmp[j];
                     ++cutNz;
                  }
               }
            }
            else if (block.cplex.getStatus().equals(IloCplex.Status.Optimal)) {
               // The subproblem has a finite optimal solution.
               // We need to check if this gives rise to an optimality cut (see
               // also the comments at the top of this file).            
               System.out.print("optimal ");
               double objval = block.cplex.getObjValue();
               double eta = x[etaind + b];
               double[] rayVals = block.cplex.getValues(block.vars);
               
               if ( objval > eta + 1e-6 ) {
                  cutub = 0.0;
                  for (int j = 0; j < block.vars.length; ++j)
                     cutub -= rayVals[j] * block.objMap.get(block.vars[j]);
                  for (Block.FixData f : block.fixed)
                     tmp[f.col] -= f.val * rayVals[f.row];
                  for (int j = 0; j < master.vars.length; ++j) {
                     if ( Math.abs(tmp[j]) > 1e-6 ) {
                        cutVal[cutNz] = tmp[j];
                        cutVar[cutNz] = master.vars[j];
                        ++cutNz;
                     }
                  }
                  cutVal[cutNz] = -1.0;
                  cutVar[cutNz] = master.vars[etaind + b];
                  ++cutNz;
               }
            }
            else {
               System.err.println("Unexpected status " +block.cplex.getStatus());
               error = true;
            }
      
            // If a cut was found then add that.
            if ( cutNz > 0 ) {
               IloLinearNumExpr expr = master.cplex.linearNumExpr();
               for (int i = 0; i < cutNz; ++i)
                  expr.addTerm(cutVal[i], cutVar[i]);
               IloRange cut = master.cplex.range(cutlb, expr, cutub);
               System.out.println("cut found: " + cut);
               add(cut);
            }
            else
               System.out.println("no cuts.");
         }
         if ( error )
            throw new RuntimeException();
      }

   }

   /** Create the dual of a linear program.
    * The function can only dualize programs of the form
    * <code>Ax <= b, x >= 0</code>. The data in <code>primalVars</code> and
    * <code>dualRows</code> as well as in <code>primalRows</code> and
    * <code>dualVars</code> is in 1-to-1-correspondence.
    * @param primalObj  Objective function of primal problem.
    * @param primalVars Variables in primal problem.
    * @param primalRows Rows in primal problem.
    * @param dual       The dual model will be stored in this modeler.
    * @param dualVars   All dual variables will be stored here.
    * @param dualRows   All dual rows will be stored here.
    */
   private static void makeDual(IloObjective primalObj,
                                IloNumVar[] primalVars,
                                IloRange[] primalRows,
                                IloCplexModeler dual,
                                Vector<IloNumVar> dualVars,
                                Vector<IloRange> dualRows) throws IloException
   {
      // To keep the code simple we only support problems
      // of the form Ax <= b, b >= 0 here. We leave it as a reader's
      // exercise to extend the function to something that can handle
      // any kind of linear model.
      for (int j = 0; j < primalVars.length; ++j)
         if ( primalVars[j].getLB() != 0 ||
              primalVars[j].getUB() < 1e20 )
            throw new RuntimeException("Cannot dualize variable " + primalVars[j]);
      for (int i = 0; i < primalRows.length; ++i)
         if ( primalRows[i].getLB() > -1e20 ||
              primalRows[i].getUB() >= 1e20 )
            throw new RuntimeException("Cannot dualize constraint " + primalRows[i]);

      // The dual of
      //   min/max c^T x
      //       Ax <= b
      //        x >= 0
      // is
      //   max/min y^T b
      //       y^T A <= c
      //           y <= 0
      // We scale y by -1 to get >= 0

      IloObjective obj = dual.objective(primalObj.getSense() == IloObjectiveSense.Minimize ?
                                        IloObjectiveSense.Maximize : IloObjectiveSense.Minimize);
      Vector<IloRange> rows = new Vector<IloRange>();
      Vector<IloNumVar> y = new Vector<IloNumVar>();
      HashMap<IloNumVar,Integer> v2i = new HashMap<IloNumVar,Integer>();
      for (int j = 0; j < primalVars.length; ++j) {
         IloNumVar x = primalVars[j];
         v2i.put(x, j);
         rows.add(dual.range(Double.NEGATIVE_INFINITY, 0, x.getName()));
      }
      for (IloLinearNumExprIterator it = ((IloLinearNumExpr)primalObj.getExpr()).linearIterator(); it.hasNext(); /* nothing */) {
         final IloNumVar v = it.nextNumVar();
         rows.elementAt(v2i.get(v)).setUB(it.getValue());
      }
      for (int i = 0; i < primalRows.length; ++i) {
         IloRange r = primalRows[i];
         IloColumn col = dual.column(obj, -r.getUB());
         for (IloLinearNumExprIterator it = ((IloLinearNumExpr)r.getExpr()).linearIterator(); it.hasNext(); /* nothing */) {
            IloNumVar v = it.nextNumVar();
            col = col.and(dual.column(rows.elementAt(v2i.get(v)), -it.getValue()));
         }
         y.add(dual.numVar(col, 0, Double.POSITIVE_INFINITY, IloNumVarType.Float, r.getName()));
      }
      dual.add(obj);
      dual.add(y.toArray(new IloNumVar[0]));
      dual.add(rows.toArray(new IloRange[0]));
      if ( dualVars != null )
         dualVars.addAll(y);
      if ( dualRows != null )
         dualRows.addAll(rows);
   }

   /** Solve the <code>problem</code> using a distributed implementation of
    * Benders' decomposition.
    */
   static boolean solve (Problem problem, String transport, String[] args, String[] machines) throws IloException {
      Vector<Block> blocks = new Vector<Block>();
      Block master = null;
      boolean result = false;
      try {
         // Extract blocks and master problem.
         System.out.println("Extracting " + problem.getNBlocks() + " blocks.");
         for (int b = 0; b < problem.getNBlocks(); ++b)
            blocks.add(new Block(problem, b, transport, args, machines));
         master = new Block(problem, blocks);
         
         // Write out the master and all blocks (for debugging).
         master.cplex.exportModel("master.lp");
         for (int b = 0; b < blocks.size(); ++b)
            blocks.elementAt(b).cplex.exportModel("block" + b + ".lp");

         // Solve the master.
         // If we find a feasible solution then perform a final solve
         // with the original problem so that we get solution values for
         // all variables in the original problem.
         if ( master.cplex.solve() ) {
            // Perform a last solve to get the solution values.
            double[] vals = master.cplex.getValues(master.vars);
            IloNumVar[] ox = new IloNumVar[master.vars.length];
            double[] olb = new double[master.vars.length];
            double[] oub = new double[master.vars.length];
            int next = 0;
            for (int i = 0; i < master.vars.length; ++i) {
               // Fix integral variables to their value in the master.
               final IloNumVar x = master.varMap.get(master.vars[i]);
               if ( x != null && x.getType() != IloNumVarType.Float ) {
                  ox[next] = x;
                  olb[next] = x.getLB();
                  oub[next] = x.getUB();
                  ++next;
                  final double bnd = Math.round(vals[i]);
                  x.setLB(bnd);
                  x.setUB(bnd);
               }
            }
            try {
               IloCplex cplex = problem.getModel();
               if ( cplex.solve() ) {
                  vals = cplex.getValues(problem.getVariables());
                  System.out.println("#### Problem solved (" + cplex.getObjValue() + ", " + cplex.getCplexStatus() + ")");
                  for (int i = 0; i < vals.length; ++i)
                     System.out.println("\tx[" + i + "] = " + vals[i]);
                  result = true;
               }
            } finally {
               // Restore original variable bounds.
               for (int i = 0; i < next; ++i) {
                  ox[i].setLB(olb[i]);
                  ox[i].setLB(oub[i]);
               }
            }
         }
      } catch (IloException e) {
         if ( master != null )
            master.end();
         while ( blocks.size() > 0 ) {
            blocks.remove(blocks.size() - 1).end();
         }
         throw e;
      }

      master.end();
      while ( blocks.size() > 0 ) {
         blocks.remove(blocks.size() - 1).end();
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
   private static final class Example implements Problem {

      private final int nblocks;
      private final IloCplex cplex;
      private final IloNumVar[] vars;
      private final IloRange[] ranges;
      private final HashMap<IloNumVar,Integer> blockMap = new HashMap<IloNumVar,Integer>();
      private final HashMap<IloNumVar,RowSet> intersectMap = new HashMap<IloNumVar,RowSet>();
      private final HashMap<IloNumVar,Double> objMap = new HashMap<IloNumVar,Double>();
      private final IloObjectiveSense objSense;

      public Example() throws IloException {
         // Model data.
         // fixed[] is the fixed cost for opening a facility,
         // cost[i,j] is the cost for serving customer i from facility j.
         final double[] fixed = new double[]{ 2.0, 3.0, 3.0 };
         final double[] cost = new double[]{ 2.0, 3.0, 4.0, 5.0, 7.0,
                                             4.0, 3.0, 1.0, 2.0, 6.0,
                                             5.0, 4.0, 2.0, 1.0, 3.0 };
         nblocks = (cost.length/fixed.length);
         cplex = new IloCplex();

         IloLinearNumExpr obj = cplex.linearNumExpr();

         // Create integer y  variables.
         Vector<IloNumVar> y = new Vector<IloNumVar>();
         for (int f = 0; f < fixed.length; ++f) {
            IloIntVar v = cplex.intVar(0, 1, "y" + f);
            obj.addTerm(fixed[f], v);
            objMap.put(v, fixed[f]);
            y.add(v);
            blockMap.put(v, -1);
            intersectMap.put(v, new RowSet());
         }

         // Create continuous x variables.
         Vector<IloNumVar> x = new Vector<IloNumVar>();
         for (int f = 0; f < fixed.length; ++f) {
            for (int c = 0; c < (cost.length/fixed.length); ++c) {
               IloNumVar v = cplex.numVar(0.0, Double.POSITIVE_INFINITY, "x" + f + "#" + c);
               obj.addTerm(cost[f * (cost.length/fixed.length) + c], v);;
               objMap.put(v, cost[f * (cost.length/fixed.length) + c]);
               x.add(v);
               blockMap.put(v, c);
               intersectMap.put(v, new RowSet());
            }
         }
         y.addAll(x);
         vars = y.toArray(new IloNumVar[0]);
         cplex.add(vars);

         // Add objective function.
         cplex.addMinimize(obj, "obj");
         objSense = IloObjectiveSense.Minimize;

         // Satisfy each customer's demand.
         final Vector<IloRange> rngs = new Vector<IloRange>();
         for (int c = 0; c < (cost.length/fixed.length); ++c) {
            IloRange r = cplex.range(1.0, Double.POSITIVE_INFINITY, "c1_" + c);
            IloLinearNumExpr lhs = cplex.linearNumExpr();
            for (int f = 0; f < fixed.length; ++f) {
               lhs.addTerm(1.0, x.elementAt(f * (cost.length/fixed.length) + c));
               intersectMap.get(x.elementAt(f * (cost.length/fixed.length) + c)).add(r);
            }
            r.setExpr(lhs);
            rngs.add(r);
         }

         // A factory must be open if we service from it.
         for (int c = 0; c < (cost.length/fixed.length); ++c) {
            for (int f = 0; f < fixed.length; ++f) {
               IloRange r = cplex.range(0.0, Double.POSITIVE_INFINITY, "c2_" + c + "#" + f);
               intersectMap.get(x.elementAt(f * (cost.length/fixed.length) + c)).add(r);
               intersectMap.get(y.elementAt(f)).add(r);
               r.setExpr(cplex.sum(cplex.prod(-1.0, x.elementAt(f * (cost.length/fixed.length) + c)),
                                   y.elementAt(f)));
               rngs.add(r);
            }
         }

         // Capacity constraint.
         IloRange r = cplex.range(Double.NEGATIVE_INFINITY, fixed.length - 1, "c3");
         IloLinearNumExpr lhs = cplex.linearNumExpr();
         for (int f = 0; f < fixed.length; ++f) {
            lhs.addTerm(1.0, y.elementAt(f));
            intersectMap.get(y.elementAt(f)).add(r);
         }
         r.setExpr(lhs);
         rngs.add(r);

         ranges = rngs.toArray(new IloRange[0]);
         cplex.add(ranges);

      }

      // Implementation of functions required by the BendersOpt::Problem
      // interface.
      public int getNBlocks() { return nblocks; }
      public IloCplex getModel() { return cplex; }
      public IloNumVar[] getVariables() { return vars; }
      public IloRange[] getRows() { return ranges; }
      public int getBlock(IloNumVar x) {
         final Integer block = blockMap.get(x);
         return block != null ? block.intValue() : -1;
      }
      private static final RowSet EMPTY_ROWSET = new RowSet();
      public RowSet getIntersectedRows(IloNumVar x) {
         final RowSet set = intersectMap.get(x);
         return set != null ? set : EMPTY_ROWSET;
      }
      public double getObjCoef(IloNumVar x) {
         final Double d = objMap.get(x);
         return d != null ? d.doubleValue() : 0.0;
      }
      public IloObjectiveSense getObjSense() { return objSense; }
   }

   // ----------------------------------------------------------------------

   public static void main(String[] args) {
      final Vector<String> transportArgs = new Vector<String>();
      final Vector<String> machines = new Vector<String>();
      String transportName = System.getProperty("ilog.cplex.transport");
      int nmachines = 0;

      // Check that a supported transport was selected.
      if ( transportName == null )
         throw new RuntimeException("Property ilog.cplex.transport not defined");

      if ( transportName.equals("mpitransport") )
         throw new RuntimeException("MPI transport not supported for Java");
      else if ( transportName.equals("processtransport") ) {
         nmachines = Integer.MAX_VALUE;
      }
      else if ( transportName.equals("tcpiptransport") ) {
         nmachines = 0;
      }
      else
         throw new RuntimeException("Unknown transport " + transportName);

      // Process the command line.
      for (String arg : args) {
         if ( transportName.equals("processtransport") ) {
            if ( arg.startsWith("-bin=") ) {
               transportArgs.add(0, arg.substring(5));
               transportArgs.add(1, "-worker=process");
               continue;
            }
         }
         else if ( transportName.equals("tcpiptransport") ) {
            if ( arg.startsWith("-address=") ) {
               machines.add(arg);
               ++nmachines;
               continue;
            }
         }
         else
            transportArgs.add(arg);
      }

      try {
         Example problem = new Example();

         if ( nmachines < problem.getNBlocks() )
            throw new RuntimeException("We have " + problem.getNBlocks()
                                       + " blocks blocks but only "
                                       + nmachines + " machines!");

         solve (problem, transportName, transportArgs.toArray(new String[0]),
                machines.toArray(new String[0]));
      } catch (IloException e) {
         System.err.println(e.getMessage());
         e.printStackTrace();
         System.exit(-1);
      }
   }
}
