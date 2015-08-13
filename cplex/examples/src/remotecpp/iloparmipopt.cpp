/* --------------------------------------------------------------------------
 * File: parmipopt.cpp
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
 *    Parallel concurrent mipopt                                          *
 *                                                                        *
 *   The code in this file runs parallel distributed concurrent mipopt.   *
 *   This means it solves the same problem on different machines, using   *
 *   a different parameter setting on each machine. In the example        *
 *   implemented here we use two different parameter settings:            *
 *   "primal only": This parameter setup configures the solver to         *
 *                  mainly work on the primal part of the problem, i.e.   *
 *                  find good feasible solutions. To this end heuristic   *
 *                  parameters are cranked up and cuts are disabled.      *
 *   "dual only":   This parameter setup configures the solver to         *
 *                  mainly work on the dual part of the problem, i.e.,    *
 *                  improve the dual bound.                               *
 *   Additional jobs will alternate "primal only" and "dual only"         *
 *   settings, yet starting with a different random seed each time.       *
 *   We run each parameter setting on a different machine and keep track  *
 *   of the best known primal and dual bounds on each machine. Once the   *
 *   best known primal and best known dual bound meet the gap stopping    *
 *   criterion we kill all solves since the optimal solution is then      *
 *   found. Note that the best primal and best dual bound do not          *
 *   necessarily have to come from the same machine. Instead the best     *
 *   primal bound is likely to come from the "primal only" solver while   *
 *   the best dual bound is likely to come from the "dual only" solver.   *
 *                                                                        *
 * ********************************************************************** */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* ********************************************************************** *
 *                                                                        *
 *  C o m m o n   s t u f f   f o r   m a s t e r   a n d   w o r k e r   *
 *                                                                        *
 * ********************************************************************** */

/** The different actions our userfunction can perform.
 * We use the userfunction to either install callbacks on the remote
 * worker or to change the minimum difference between consecutive objective
 * functions that the worker reports to the master.
 */
enum {
   USERACTION_ADDCALLBACK,    /**< Install the info callback in the worker. */
   USERACTION_REMOVECALLBACK, /**< Remove the info callback from the worker. */
   USERACTION_CHANGEOBJDIFF   /**< Change the minimal difference for objective
                               *   changes. Consecutive objective values must
                               *   differ by at least this amount to be
                               *   considered different. */
};

/** The different types of info messages that workers can send.
 * Workers may inform the master about new primal or dual bounds they found.
 * They may also send their deterministic time so that the master can check
 * in a deterministic way how much work each worker did so far.
 */
enum {
   INFO_NEWDUAL,   /**< Reports a new dual bound. */
   INFO_NEWPRIMAL, /**< Reports a new primal bound. */
   INFO_DETTIME    /**< Reports the current deterministic time stamp. */
};

/* ********************************************************************** *
 *                                                                        *
 *    W o r k e r   s i d e   u s e r f u n c t i o n                     *
 *                                                                        *
 *    On the worker we register a user function that installs an info     *
 *    callback. That callback reports new bounds to the master so that    *
 *    the master can monitor progress.                                    *
 *    Note that the userfunction and the callback must interact with      *
 *    the callable library's C API, not with the C++ API. This is         *
 *    because the remote worker does not have a Concert model.            *
 *                                                                        *
 * ********************************************************************** */

#if defined(COMPILE_USERFUNCTION)
#if defined(__hpux)
#define ILOUSEMT
#endif
extern "C" {
#include <ilcplex/cplexremoteworkerx.h>
}

/** Best known primal and dual bounds in this remote worker.
 *  We need to keep track of these values so that we only report bound
 *  improvements to the master process.
 */
static struct {
   bool            haveDual;   /**< Is the 'dual' field valid? */
   bool            havePrimal; /**< Is the 'primal' field valid? */
   double          dual;       /**< Best known dual bound in this solver. */
   double          primal;     /**< Best known primal bound in this solver. */
   double          objdiff;    /**< Minimal delta between consecutive bounds.
                                *   If two consecutive bounds differ by less
                                *   than this value they are not considered
                                *   to have changed. */
} best;

extern "C" {
/** MIP info callback that is registered with CPLEX.
 *  This callback picks up primal and dual bounds as well as the current
 *  deterministic time. If bounds changed then updated bounds are send
 *  to the master. Deterministic time stamps are always send to the master.
 */
static int CPXPUBLIC
infocallback (CPXCENVptr cbenv, void *cbdata, int wherefrom, void *cbhandle)
{
   CPXCENVptr env = static_cast<CPXCENVptr>(cbhandle);
   double dual, primal, ts;

   // Test if we have improved the primal bound and report if so.
   if ( CPXXgetcallbackinfo (cbenv, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &primal) == 0 ) {
      if ( !best.havePrimal || fabs (best.primal - primal) > best.objdiff ) {
         best.havePrimal = true;
         best.primal = primal;
         (void)CPXXsendinfodouble (env, INFO_NEWPRIMAL, 1, &primal);
      }
   }

   // Test if we have improved the dual bound and report if so.
   if ( CPXXgetcallbackinfo (cbenv, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &dual) == 0 ) {
      if ( !best.haveDual || fabs (best.dual - dual) > best.objdiff ) {
         best.haveDual = true;
         best.dual = dual;
         (void)CPXXsendinfodouble (env, INFO_NEWDUAL, 1, &dual);
      }
   }

   // Always report the current deterministic time.
   if ( CPXXgetdettime (cbenv, &ts) == 0 ) {
      (void)CPXXsendinfodouble (env, INFO_DETTIME, 1, &ts);
   }

   return 0;
}

/** User function implementation.
 *  This function is executed when the master invokes a user function
 *  on this remote solver.
 */
static int CPXPUBLIC
userfunction (CPXENVptr env, int id, CPXLONG inlen,
              void const *indata, CPXLONG maxout,
              CPXLONG *outlen_p, void *outdata, void *handle)
{
   CPXDESERIALIZERptr d = NULL;
   int status = 0;

   (void)maxout;
   (void)outdata;
   (void)handle;

   *outlen_p = 0;
   CPXXdeserializercreate (&d, inlen, indata);

   switch (id) {
   case USERACTION_ADDCALLBACK:
      best.havePrimal = false;
      best.haveDual = false;
      (void)CPXXsetinfocallbackfunc (env, infocallback, env);
      break;
   case USERACTION_REMOVECALLBACK:
      (void)CPXXsetinfocallbackfunc (env, NULL, NULL);
      break;
   case USERACTION_CHANGEOBJDIFF:
      {
         double dbl;

         dbl = best.objdiff;
         d->getdouble (d, &dbl);
         best.objdiff = dbl;
      }
      break;
   }

   CPXXdeserializerdestroy (d);

   return status;
}

/** Register user function handler.
 *  This function is invoked by the remote worker at startup to allow us
 * to register a user function handler.
 */
CPXEXPORT void CPXPUBLIC REGISTER_USERFUNCTION (struct messagehandler *handler);
CPXEXPORT void CPXPUBLIC
REGISTER_USERFUNCTION (struct messagehandler *handler)
{
   fprintf (stderr, "Registering user function.\n");
   CPXXsetuserfunction (handler, userfunction, NULL);
}
} // extern "C"

#endif // COMPILE_USERFUNCTION

/* ********************************************************************** *
 *                                                                        *
 *    M a s t e r   s i d e   i m p l e m e n t a t i o n                 *
 *                                                                        *
 *    Implement of the master consists of the following major parts:      *
 *    - class ParamValue. This class is used to setup (sets of)           *
 *                        predefined parameter settings.                  *
 *    - class SolveState. There is only one global instance of this       *
 *                        class that represents the progress and current  *
 *                        state of the solve. It keeps track of the       *
 *                        best known primal and dual bounds over all      *
 *                        workers.                                        *
 *    - class CustumOutput. This class is used in the special case in     *
 *                        which we want to see the log from the workers   *
 *                        but want to prefix each line by the             *
 *                        respective worker's index.                      *
 *    - class Worker.     This class represents a solve that is running   *
 *                        on a remote object worker.                      *
 *    - main().           The main function does three things:            *
 *                        1. Parse the command line.                      *
 *                        2. Create an instance of Worker for each        *
 *                           configured remote object worker (the         *
 *                           constructor implicitly starts an             *
 *                           asynchronous solve).                         *
 *                        3. Wait until optimality was proven.            *
 *                                                                        *
 * ********************************************************************** */

#if defined(COMPILE_MASTER)

#include <ilcplex/ilocplex.h>


#include <limits.h>
#ifdef _WIN32
#   include <windows.h>
#   include <direct.h>
#   define MAX_PATH_LEN MAX_PATH
#   define getcwd _getcwd
#   define millisleep(x) Sleep(x)
#else
#   include <unistd.h>
#   define MAX_PATH_LEN PATH_MAX
#   define millisleep(x) usleep((x) * 1000)
#endif
#ifdef USE_MPI
#   define OMPI_SKIP_MPICXX 1 // We don't use the C++ bindings of OpenMPI.
#   include <mpi.h>
#endif

// -------------------- class ParamValue ------------------------------

/** Class to pre-define parameter settings.
 * An instance of this class specifies a pre-defined parameter setting
 * for one particular parameter.
 */
class ParamValue {
   char type; /**< Parameter type. */
   int num;   /**< Parameter number. */
   union {
      CPXINT  i;
      CPXLONG l;
      IloNum  d;
   } value;   /**< The value for the parameter. */
public:
   /** Special constructor to define an "invalid" setting that can be
    * used as end marker in arrays of ParamValue.
    */
   ParamValue() : type(' '), num(-1) {}
   /** Create a pre-defined integer parameter. */
   ParamValue(IloCplex::IntParam n, CPXINT i) : type('i'), num(n) {
      value.i = i;
   }
   /** Create a pre-defined long parameter. */
   ParamValue(IloCplex::LongParam n, CPXLONG l) : type('l'), num(n) {
      value.l = l;
   }
   /** Create a pre-defined double parameter. */
   ParamValue(IloCplex::NumParam n, IloNum d) : type('d'), num(n) {
      value.d = d;
   }
   /** Create a pre-defined boolean parameter. */
   ParamValue(IloCplex::BoolParam n, bool b) : type('b'), num(n) {
      value.i = b ? 1 : 0;
   }
   /** Test if this is a valid parameter setting or one that was created
    * by the special "invalid" constructor.
    */
   bool isValid() const { return num >= 0; }
   /** Apply the setting defined by this instance to <code>cplex</code>. */
   void apply(IloCplex cplex) const {
      switch (type) {
      case 'b': cplex.setParam(IloCplex::BoolParam(num), value.i != 0); break;
      case 'i': cplex.setParam(IloCplex::IntParam(num), value.i); break;
      case 'l': cplex.setParam(IloCplex::LongParam(num), value.l); break;
      case 'd': cplex.setParam(IloCplex::NumParam(num), value.d); break;
      }
   }
};

/** Parameter settings to have CPLEX work mainly on the primal bound. */
static ParamValue const primalbound[] = {
   ParamValue(IloCplex::Param::MIP::Strategy::HeuristicFreq, 1), // Run heuristics at each node.
   ParamValue(IloCplex::Param::MIP::Strategy::RINSHeur, 1), // Run RINS every two nodes
   ParamValue(IloCplex::Param::MIP::Limits::CutPasses, -1), // Disable cuts
   ParamValue()
};
/** Parameter settings to have CPLEX work exclusively on the dual bound. */
static ParamValue const dualbound[] = {
   ParamValue(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1), // Heuristics off.
   // All cuts aggressive.
   ParamValue(IloCplex::Param::MIP::Cuts::Cliques, 3),
   ParamValue(IloCplex::Param::MIP::Cuts::Covers, 3),
   ParamValue(IloCplex::Param::MIP::Cuts::Disjunctive, 3),
   ParamValue(IloCplex::Param::MIP::Cuts::FlowCovers, 2),
   ParamValue(IloCplex::Param::MIP::Cuts::Gomory, 2),
   ParamValue(IloCplex::Param::MIP::Cuts::GUBCovers, 2),
   ParamValue(IloCplex::Param::MIP::Cuts::Implied, 2),
   ParamValue(IloCplex::Param::MIP::Cuts::MIRCut, 2),
   ParamValue(IloCplex::Param::MIP::Cuts::ZeroHalfCut, 2),
   ParamValue(IloCplex::Param::MIP::Cuts::MCFCut, 2),
   ParamValue()
};

/** Predefined parameter settings. */
static struct {
   char const        *name;
   class ParamValue const *values;
} const settings[] = {
   { "primal only", primalbound },
   { "dual bound",  dualbound }
};
#define NUMSETTINGS ((int)(sizeof (settings) / sizeof (settings[0])))

// -------------------- class SolveState ------------------------------

/** Current state of a solve.
 * This structure keeps track of the best known primal and dual bounds.
 */
struct SolveState {
   struct Bound {
      bool volatile    valid; /**< True if the bound/env fields are valid. */
      double volatile  bound; /**< Best known primal bound. */
      int volatile     idx;   /**< The environment that reported bound. */
      Bound(bool primal) : valid(false), bound(0.0), idx(-1), isPrimal(primal)
      {}
      void update(double b, int i, bool isMax) {
         if ( !valid
              || (isPrimal && ((isMax && b > bound) || (!isMax && b < bound)))
              || (!isPrimal && ((isMax && b < bound) || (!isMax && b > bound)))
              )
         {
            idx = i;
            bound = b;
            valid = true;
         }
      }
   private:
      bool const isPrimal;
   };
   Bound primal;
   Bound dual;
   SolveState() : primal(true), dual(false) {}
};

// -------------------- class CustomOutput ----------------------------

/** Stream buffer used for custom output.
 * We use this stream buffer if we want to print the CPLEX log messages
 * from each worker and want to prepend them with the workers index.
 * The class derives from std::streambuf and overrides the bare minimum
 * of functions to get the desired behavior.
 */
struct CustomOutput : public std::streambuf {
   /** Create a new instance.
    * @param p The index of the worker. This value is written as first
    *          thing on each line printed.
    */
   CustomOutput(int p) : prefix(p), startLine(true) {}
protected:
   /** Write <code>n</code> characters from <code>s</code> and return
    * the number of characters written.
    */
   std::streamsize xsputn (const char* s, std::streamsize n) {
      for (std::streamsize i = 0; i < n; ++i)
         writechar(s[i]);
      return n;
   }
   /** Write a single character <code>c</code> and return the character
    * just written (or <code>EOF</code> on error.
    */
   int overflow (int c = EOF) {
      if ( c == EOF )
         return EOF;
      else {
         writechar(static_cast<char>(c));
         return traits_type::to_int_type(c);
      }
   }
private:
   int const prefix; /**< Output prefix. */
   bool startLine;   /**< Are we at the start of a line. */
   /* Write a single character <code>c</code> to standard output.
    */
   void writechar(char c) {
      // If we are at the start of a line then first write the prefix.
      if ( startLine ) {
         std::cout << "[" << prefix << "] ";
         startLine = false;
      }
      // Write character.
      std::cout << c;
      // If we just wrote a line break then we are now at the start
      // of the next line.
      if ( c == '\n' )
         startLine = true;
   }
};

// -------------------- class Worker ----------------------------------

/** A remote worker.
 */
class Worker {
   /** Handler for info messages from the remote worker.
    * Instances of this class handle info messages that are sent from
    * the worker (from function infocallback() defined above). Whenever
    * such a message arrives function main() is invoked.
    * The information reported by remote workers are:
    * - new dual bounds,
    * - new primal bounds,
    * - current deterministic timestamps.
    * The main() function updates the respective fields in the correspoding
    * Worker instance and also updates that global best known primal/dual
    * bounds if necessary.
    */
   struct InfoHandler : public IloCplex::RemoteInfoHandler {
      Worker *const worker;
      /** Create a new info message handler for Worker <code>w</code>. */
      InfoHandler(Worker *w) : worker(w) {}
      void main(CPXINFOTYPE type, int tag, CPXLONG length,
                void const *data)
      {
         SolveState *const s = worker->state;
         double d;

         switch (tag) {
         case INFO_NEWDUAL:
            assert (type == CPXINFO_DOUBLE);
            assert (length == 1);
            d = *static_cast<double const *>(data);
            std::cout << "[" << worker->idx << "] New dual bound: " << d
                      << std::endl;
            worker->dual = d;
            s->dual.update(d, worker->idx,
                           worker->getObjectiveSense() == IloObjective::Maximize);
            break;
         case INFO_NEWPRIMAL:
            assert (type == CPXINFO_DOUBLE);
            assert (length == 1);
            d = *static_cast<double const *>(data);
            std::cout << "[" << worker->idx << "] New primal bound: " << d
                      << std::endl;
            worker->primal = d;
            s->primal.update(d, worker->idx,
                             worker->getObjectiveSense() == IloObjective::Maximize);
            break;
         case INFO_DETTIME:
            assert (type == CPXINFO_DOUBLE);
            assert (length == 1);
            d = *static_cast<double const *>(data);
            std::cout << "[" << worker->idx << "] Current dettime: " << d
                      << " (" << worker->dual << ", " << worker->primal << ")"
                      << std::endl;
            worker->dettime = d;
            break;
         default:
            std::cout << "[" << worker->idx << "] Unknown info " << tag
                      << std::endl;
         }
      }
   };
   friend class InfoHandler;



   int const idx;                 /**< Index of worker. */
   SolveState *const state;       /**< Pointer to global solve state. */
   IloModel model;                /**< Model solve on this worker. */
   IloCplex cplex;                /**< Solver used by this worker. */
   IloCplex::SolveHandle handle;  /**< Handle for asynchronous solve. */
   double dettime;                /**< Last dettime stamp. */
   double primal;                 /**< Best known primal bound. */
   double dual;                   /**< Best known dual bound. */
   IloObjective obj;              /**< Objective function from model. */
   IloNumVarArray x;              /**< All variables in model. */
   IloRangeArray rng;             /**< All linear constraints in model. */
   InfoHandler infoHandler;       /**< Handler for info messages. */
   CustomOutput outb;             /**< Stream buffer for custom output. */
   std::ostream outs;             /**< Output stream for custom output. */
public:
   /** The different log message output modes. */
   typedef enum {
      OUTPUT_SILENT,   /**< No log messages are printed. */
      OUTPUT_PREFIXED, /**< Log message are printed and prefixed by worker id. */
      OUTPUT_LOG       /**< Log message are printed without prefix. */
   } OUTPUT;

   /** Create a new worker.
    * The constructor mainly does the following:
    * - create an IloCplex instance that refers to a remote worker,
    * - load the model in <code>modelfile</code>,
    * - setup parameters depending on this worker's index,
    * - start an asynchronous solve.
    * If anything fails then an exception will be thrown.
    * @param env The environment used for instantiating Ilo* objects.
    * @param i The index of the worker to be created. This also
    *          determines the parameter settings to use in this worker.
    * @param s  A pointer to the global solve state.
    * @param transport  The transport name for the IloCplex constructor.
    * @param argc       The argument count for the IloCplex constructor.
    * @param argv       The array of transport arguments for the IloCplex
    *                   constructor.
    * @param modelfile  Name of the model to be loaded into the worker.
    * @param output     The output mode.
    * @param objdiff    The minimal difference between so that two
    *                   consecutive objective function values are considered
    *                   different.
    */
   Worker(IloEnv env, int i, SolveState *s, char const *transport,
          int argc, char const **argv, char const *modelfile,
          OUTPUT output, double objdiff)
      : idx(i), state(s), model(env), cplex(0), handle(0),
        primal(IloInfinity), dual(-IloInfinity),
        obj(env), x(env), rng(env), infoHandler(this), outb(idx), outs(&outb)
   {
      try {
         // Create remote object, setup output and load the model.
         cplex = IloCplex(model, transport, argc, argv);
         switch (output) {
         case OUTPUT_SILENT:
            // Disable output on the output and warning stream.
            cplex.setOut(env.getNullStream());
            cplex.setWarning(env.getNullStream());
            break;
         case OUTPUT_PREFIXED:
            // Redirect output to our custom stream.
            cplex.setOut(outs);
            cplex.setWarning(outs);
            break;
         case OUTPUT_LOG:
            // Nothing to do here. By default output is enabled.
            break;
         }
         cplex.importModel(model, modelfile, obj, x, rng);
         if ( obj.getSense() == IloObjective::Minimize ) {
            primal = -IloInfinity;
            dual = IloInfinity;
         }

         // We set the thread count for each solver to 1 so that we do not
         // run into problems if multiple solves are performed on the same
         // machine.
         cplex.setParam(IloCplex::Param::Threads, 1);
         // Each worker runs with a different random seed. This way we
         // get different paths through the tree even if the other
         // parameter settings are the same.
         cplex.setParam(IloCplex::Param::RandomSeed, idx);
         // Apply parameter settings.
         for (class ParamValue const *vals = settings[idx % NUMSETTINGS].values;
              vals->isValid(); ++vals)
            vals->apply(cplex);

         // Install callback and set objective change.
         int status = cplex.userfunction (USERACTION_ADDCALLBACK,
                                          0, NULL, 0, 0, NULL);
         if ( status )
            throw status;
         IloCplex::Serializer s;
         s.add(objdiff);
         status = cplex.userfunction (USERACTION_CHANGEOBJDIFF,
                                      s.getRawLength(), s.getRawData(),
                                      0, 0, NULL);
         if ( status )
            throw status;

         // Register the handler that will process info messages sent
         // from the worker.
         cplex.setRemoteInfoHandler(&infoHandler);

         // Everything is setup. Launch the asynchronous solve.
         handle = cplex.solve(true);
      } catch (...) {
         // In case of an exception we need to take some special
         // cleanup actions. Note that if we get here then the
         // solve cannot have been started and we don't need to
         // kill or join the asynchronous solve.
         if ( cplex.getImpl() )
            cplex.end();
         rng.end();
         x.end();
         obj.end();
         model.end();
         throw;
      }
   }

   /** Destructor. */
   ~Worker() {
      if ( handle ) {
         handle.kill();
         handle.join();
      }
      if ( cplex.getImpl() ) {
         cplex.userfunction(USERACTION_REMOVECALLBACK, 0, NULL, 0, 0, NULL);
         cplex.end();
      }
      rng.end();
      x.end();
      obj.end();
      model.end();
   }

   /** Get the best incumbent in this worker.
    * Calling this function before calling join() is an error.
    */
   IloNumArray getX() const {
      IloNumArray values(cplex.getEnv());
      cplex.getValues(x, values);
      return values;
   }

   /** Get the best integer feasible objective value in this worker.
    * Calling this function before calling join() is an error.
    */
   double getObjective() const {  return cplex.getObjValue(); }

   /** Get the solution status for this worker.
    * Calling this function before calling join() is an error.
    */
   IloAlgorithm::Status getStatus() const { return cplex.getStatus(); }

   /** Get the current best dual bound in this worker. */
   double getDual() const { return dual; }

   /** Get the current best primal bound in this worker. */
   double getPrimal() const { return primal; }

   /** Get the current deterministic time stamp for this worker. */
   double getDetTime() const { return dettime; }

   /** Get the objective sense for this worker's objective. */
   IloObjective::Sense getObjectiveSense() const { return obj.getSense(); }

   /** Test if this worker is still running. */
   bool isRunning() const {
      if ( handle.test() )
         return true;
      return false;
   }

   /** Kill the solve carried out by this worker. */
   void kill() { handle.kill(); }

   /** Join this worker.
    * Waits until the asynchronous solve performed by this worker is complete.
    * @return <code>true</code> if the worker found a feasible solution,
    *         <code>false</code> otherwise.
    */
   bool join() {
      bool result = false;
      result = handle.joinSolve();
      handle = 0;
      return result;
   }
};

// -------------------- main() ----------------------------------------

int
main(int argc, char **argv)
{
   char const *modelfile = NULL;
   int jobs = 0;
   char const **machine = 0;

#if defined(USE_MPI)
   int numprocs, rank;

   // Initialize MPI.
   // We must have at least three processors (master and 2 workers) and
   // the master must have rank 0.
   MPI_Init(&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
   if ( numprocs < 3 ) {
      fprintf (stderr, "Invalid number of processors (%d)\n", numprocs);
      abort ();
   }
   MPI_Comm_rank (MPI_COMM_WORLD, &rank);
   if ( rank != 0 ) {
      fprintf (stderr, "Master must have rank 0!\n");
      MPI_Finalize ();
      abort ();
   }
   machine = new char const*[numprocs];
   for (int i = 0; i < numprocs; ++i)
      machine[i] = "mpimachine";
   jobs = numprocs - 1;
#elif defined(USE_PROCESS)
   // The default binary is CPLEX but that can be overwritten
   // by command line arguments.
   char const *bin = "cplex";

   machine = new char const *[argc];

#elif defined(USE_TCPIP)
   machine = new char const *[argc];
#else
#   error "No transport type selected"
#endif


   // Parse the command line.
   Worker::OUTPUT output = Worker::OUTPUT_SILENT;
   double absgap = 1e-6;
   for (int i = 1; i < argc; ++i) {
      if ( strncmp (argv[i], "-model=", 7) == 0 )
         modelfile = argv[i] + 7;
#if defined(USE_MPI)
      // no MPI specific commands supported
#elif defined(USE_PROCESS)
      // For process transport
      // - machine=<command> specifies the name of a machine to which
      //   we connect (via ssh)
      // -bin=<binary> specifies the binary to execute (default is "cplex")
      else if ( strncmp (argv[i], "-machine=", 9) == 0 )
         machine[jobs++] = argv[i] + 9;
      else if ( strncmp (argv[i], "-bin=", 5) == 0 )
         bin = argv[i] + 5;
#elif defined(USE_TCPIP)
      // For TCP/IP process
      // -address=<host>:<port> specifies a worker address to which to
      //  connect
      else if ( strncmp (argv[i], "-address=", 9) == 0 )
         machine[jobs++] = argv[i];
#endif
      // Further arguments
      // -absgap=<gap>     stop if that absolute gap is reached
      // -output-prefixed  prefix all worker output by the worker number
      // -output-log       output worker log messages
      else if ( strncmp (argv[i], "-absgap=", 8) == 0 )
         absgap = strtod (argv[i] + 8, NULL);
      else if ( strcmp (argv[i], "-output-prefixed") == 0 )
         output = Worker::OUTPUT_PREFIXED;
      else if ( strcmp (argv[i], "-output-log") == 0 )
         output = Worker::OUTPUT_LOG;
   }

   // Validate arguments.
   if ( !modelfile )
      throw "No model file specified with -model=<modelfile>";
   if ( jobs < 1 )
      throw "Invalid job count";
   
   // Find the place at which to find the shared object that implements
   // the user function. On Windows the path to the current directory is
   // likely to contain blanks, so better quote it.
   char cwd[MAX_PATH_LEN];
   char usrfunc[MAX_PATH_LEN];
   getcwd (cwd, sizeof (cwd));
   usrfunc[0] = 0;
#ifdef _WIN32
   strcat (usrfunc, "-libpath=\"");
   strcat (usrfunc, cwd);
   strcat (usrfunc, "\"");
#else
   strcat (usrfunc, "-libpath=");
   strcat (usrfunc, cwd);
#endif

   IloEnv env;
   SolveState state;

   // Initialize the workers.
   // The main thing to do here is to set up the connection arguments
   // for the IloCplex constructor. Once we have them we just instantiate
   // the Worker class. The constructor for this class will also start
   // an asynchronous solve immediately.
   Worker **workers = new Worker*[jobs];
   for (int i = 0; i < jobs; ++i) {
      char const *transport = 0;
      char const *args[16];
      int nextarg = 0;

#if defined(USE_MPI)
      char rankbuf[256];

      sprintf (rankbuf, "-remoterank=%d", i + 1);
      transport = "mpitransport";
      args[nextarg++] = rankbuf;
#elif defined(USE_PROCESS)
      char logbuf[1024];

      transport = "processtransport";
      // If the machine is not "localhost" then use ssh to connect to
      // this machine. Otherwise just fork a process on the local
      // machine.
      if ( machine[i] != NULL && strcmp (machine[i], "localhost") != 0 ) {
         args[nextarg++] = "/usr/bin/ssh";
         args[nextarg++] = machine[i];
      }
      args[nextarg++] = bin;
      args[nextarg++] = "-worker=process";
      if ( machine[i] != NULL )
         args[nextarg++] = "-stdio";
      else
         args[nextarg++] = "-namedpipes=.";
      args[nextarg++] = usrfunc;
      args[nextarg++] = "-userfunction=iloparmipopt_userfunction=REGISTER_USERFUNCTION";
      sprintf (logbuf, "-logfile=server%d.log", i);
      args[nextarg++] = logbuf;
#elif defined(USE_TCPIP)
      transport = "tcpiptransport";
      args[nextarg++] = machine[i];
#endif


      std::cout << "Initializing worker for " << machine[i] << std::endl;
      try {
         workers[i] = new Worker(env, i, &state,
                                 transport, nextarg, args,
                                 modelfile, output, 1e-5);
      } catch (...) {
         while (--i >= 0)
            delete workers[i];
         delete[] workers;
         throw;
      }
   }
   delete[] machine;

   try {
      // At this point all workers have been started and are
      // solving the problem. We just wait until either the first
      // worker has finished or if the best known global primal and dual
      // bounds are close enough.
      IloObjective::Sense const objsen = workers[0]->getObjectiveSense();
      int active = jobs;
      int frequency = 10000; // Print current bounds every two seconds.
      while (active > 0) {
         // Check if we should stop all solves.
         // We stop them if the absolute mipgap is reached.
         if ( state.primal.valid && state.dual.valid &&
              ((objsen == IloObjective::Minimize &&
                state.dual.bound + absgap >= state.primal.bound) ||
               (objsen == IloObjective::Maximize &&
                state.dual.bound - absgap <= state.primal.bound)) )
         {
            std::cout << "Stopping criterion reached. Stopping all pending solves."
                      << std::endl;
            for (int i = 0; i < jobs; ++i)
               workers[i]->kill();
            break;
         }
         if ( --frequency == 0 ) {
            std::cout << "dual=" << state.dual.bound << ", "
                      << "primal=" << state.primal.bound
                      << std::endl;
            frequency = 10000;
         }

         // Loop over all solvers and test if they are still running.
         for (int i = 0; i < jobs; ++i) {
            if ( !workers[i]->isRunning() ) {
               // The job is finished. We have a solution, so kill all
               // others.
               --active;
               std::cout << "First job (" << i << ") is finished, killing the rest"
                         << std::endl;
               for (int j = 0; j < jobs; ++j) {
                  if ( j != i ) {
                     workers[j]->kill();
                     --active;
                  }
               }
               break;
            }
         }
         // Sleep a little so that we do not poll the workers like crazy.
         millisleep (10);
      }

      // All workers have finished or have been killed. Join them.
      // For each worker we print its status, its dettime and its bounds.
      for (int i = 0; i < jobs; ++i) {
         double obj = -IloInfinity;
         bool const result = workers[i]->join();
         if ( !result ) {
            // No feasible solution found (yet) on this machine.
            // Just set objective function to a very big value
            obj = (objsen == IloObjective::Minimize) ? IloInfinity : -IloInfinity;
         }
         else {
            obj = workers[i]->getObjective();
         }

         std::cout << "Job " << i << ": " << obj << ", stat " << workers[i]->getStatus()
                   << std::endl
                   << "\t" << workers[i]->getDetTime() << ", "
                   << workers[i]->getDual() << ", "
                   << workers[i]->getPrimal() << std::endl;
      }

      // Fetch the x vector from the solver that produced the best
      // primal bound.
      int const bestidx = state.primal.idx;
      if ( bestidx < 0 ) {
         std::cout << "No solution (model infeasible)" << std::endl;
      }
      else {
         IloNumArray x = workers[bestidx]->getX();
         std::cout << "Optimal solution:" << std::endl;
         for (IloInt c = 0; c < x.getSize(); ++c)
            std::cout << "x[" << c << "]: " << x[c] << std::endl;
         x.end();
      }

      // Release the workers.
      for (int i = jobs - 1; i >= 0; --i)
         delete workers[i];
      delete[] workers;

      env.end();
   } catch (...) {
      // In case of any error we still need to delete the workers.
      // This is to make sure that we properly disconnect from the
      // remote workers. The destructor of a worker will automatically
      // kill and join the worker if that was not already done.
      for (int i = jobs - 1; i >= 0; --i)
         delete workers[i];
      delete[] workers;
      throw;
   }

#ifdef USE_MPI
   MPI_Finalize ();
#endif


   return 0;
}

#endif // COMPILE_MASTER
