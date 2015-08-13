/* --------------------------------------------------------------------------
 * File: parmipopt.c
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
 *   criterion we kill both solves since the optimal solution is then     *
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
 *    S t u f f   t h a t   i s   c o m m o n   t o   b o t h             *
 *                                                                        *
 * ********************************************************************** */

enum {
   USERACTION_ADDCALLBACK,    /**< Install the info callback in the worker. */
   USERACTION_REMOVECALLBACK, /**< Remove the info callback from the worker. */
   USERACTION_CHANGEOBJDIFF   /**< Change the minimal difference for objective
                               *   changes. Consecutive objective values must
                               *   differ by at least this amount to be
                               *   considered different. */
};

enum {
   INFO_NEWDUAL,   /**< Reports a new dual bound. */
   INFO_NEWPRIMAL, /**< Reports a new primal bound. */
   INFO_DETTIME    /**< Reports the current deterministic time stamp. */
};

/* ********************************************************************** *
 *                                                                        *
 *    S e r v e r   s i d e   u s e r f u n c t i o n                     *
 *                                                                        *
 *    On the server we register a user function that installs an info     *
 *    callback. That callback reports new bounds to the client so that    *
 *    the client can monitor progress.                                    *
 *                                                                        *
 * ********************************************************************** */

#if defined(COMPILE_USERFUNCTION)

#include <ilcplex/cplexremoteworkerx.h>

#ifdef _WIN32
#   include <windows.h>
#   define MUTEX CRITICAL_SECTION
#   define MUTEX_INIT(mtx)    (InitializeCriticalSection (mtx), 0)
#   define MUTEX_DESTROY(mtx) (DeleteCriticalSection (mtx), 0)
#   define MUTEX_LOCK(mtx)    (EnterCriticalSection (mtx), 0)
#   define MUTEX_UNLOCK(mtx)  (LeaveCriticalSection (mtx), 0)
#else
#   include <pthread.h>
#   define MUTEX pthread_mutex_t
#   define MUTEX_INIT(mtx)    pthread_mutex_init ((mtx), NULL)
#   define MUTEX_DESTROY(mtx) pthread_mutex_destroy (mtx)
#   define MUTEX_LOCK(mtx)    pthread_mutex_lock (mtx)
#   define MUTEX_UNLOCK(mtx)  pthread_mutex_unlock (mtx)
#endif

/** Best known primal and dual bounds in this remote worker.
 *  We need to keep track of these values so that we only report bound
 *  improvements to the master process.
 */
static struct {
   MUTEX           objmutex;   /**< Mutex for synchronizing access to
                                *   this structure. */
   int             haveDual;   /**< Is the 'dual' field valid? */
   int             havePrimal; /**< Is the 'primal' field valid? */
   double          dual;       /**< Best known dual bound in this solver. */
   double          primal;     /**< Best known primal bound in this solver. */
   double          objdiff;    /**< Minimal delta between consecutive bounds.
                                *   If two consecutive bounds differ by less
                                *   than this value they are not considered
                                *   to have changed. */
} best;

/** MIP info callback that is registered with CPLEX.
 *  This callback picks up primal and dual bounds as well as the current
 *  deterministic time. If bounds changed then updated bounds are send
 *  to the master. Deterministic time stamps are always send to the master.
 */
static int CPXPUBLIC
infocallback (CPXCENVptr cbenv, void *cbdata, int wherefrom, void *cbhandle)
{
   CPXCENVptr env = cbhandle;
   double dual, primal, ts;

   /* Test if we have improved the primal bound and report if so. */
   if ( CPXXgetcallbackinfo (cbenv, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &primal) == 0 ) {
      MUTEX_LOCK (&best.objmutex);
      if ( !best.havePrimal || fabs (best.primal - primal) > best.objdiff ) {
         best.havePrimal = 1;
         best.primal = primal;
         (void)CPXXsendinfodouble (env, INFO_NEWPRIMAL, 1, &primal);
      }
      MUTEX_UNLOCK (&best.objmutex);
   }

   /* Test if we have improved the dual bound and report if so. */
   if ( CPXXgetcallbackinfo (cbenv, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &dual) == 0 ) {
      MUTEX_LOCK (&best.objmutex);
      if ( !best.haveDual || fabs (best.dual - dual) > best.objdiff ) {
         best.haveDual = 1;
         best.dual = dual;
         (void)CPXXsendinfodouble (env, INFO_NEWDUAL, 1, &dual);
      }
      MUTEX_UNLOCK (&best.objmutex);
   }

   /* Always report the current deterministic time. */
   if ( CPXXgetdettime (cbenv, &ts) == 0 ) {
      MUTEX_LOCK (&best.objmutex);
      (void)CPXXsendinfodouble (env, INFO_DETTIME, 1, &ts);
      MUTEX_UNLOCK (&best.objmutex);
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
      best.havePrimal = 0;
      best.haveDual = 0;
      if ( MUTEX_INIT (&best.objmutex) )
         status = -1;
      else
         (void)CPXXsetinfocallbackfunc (env, infocallback, env);
      break;
   case USERACTION_REMOVECALLBACK:
      (void)CPXXsetinfocallbackfunc (env, NULL, NULL);
      status = MUTEX_DESTROY (&best.objmutex);
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

#endif /* COMPILE_USERFUNCTION */

/* ********************************************************************** *
 *                                                                        *
 *    C l i e n t   s i d e   i m p l e m e n t a t i o n                 *
 *                                                                        *
 * ********************************************************************** */

#if defined(COMPILE_MASTER)

#include <ilcplex/cplexremotemasterx.h>


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
#   include <mpi.h>
#endif

struct paramvalue {
   int num;  /**< Parameter number. */
   int type; /**< Parameter type. */
   struct {
      CPXINT  i;
      CPXLONG l;
      double  d;
   } value;  /**< The value for the parameter. This should better be a
              *   union but some compilers do not support a simple way of
              *   statically initializing a union so we use a struct.
              */
};
#define INTPARAM(num,val)    { num, CPX_PARAMTYPE_INT, { val, 0, 0 } }
#define LONGPARAM(num,val)   { num, CPX_PARAMTYPE_LONG, { 0, val, 0 } }
#define DOUBLEPARAM(num,val) { num, CPX_PARAMTYPE_DOUBLE, { 0, 0, val } }
#define ENDPARAM             { -1,  CPX_PARAMTYPE_NONE, { 0, 0, 0 } }

/** Parameter settings to have CPLEX work mainly on the primal bound. */
static struct paramvalue primalbound[] = {
   LONGPARAM (CPXPARAM_MIP_Strategy_HeuristicFreq, 1), /* Run heuristics at each node. */
   LONGPARAM (CPXPARAM_MIP_Strategy_RINSHeur, 2), /* Run RINS every two nodes. */
   LONGPARAM (CPXPARAM_MIP_Limits_CutPasses, -1), /* Disable cuts. */
   ENDPARAM
};
/** Parameter settings to have CPLEX work exclusively on the dual bound. */
static struct paramvalue dualbound[] = {
   LONGPARAM (CPXPARAM_MIP_Strategy_HeuristicFreq, -1), /* Heuristics off. */
   /* All cuts aggressive. */
   INTPARAM (CPXPARAM_MIP_Cuts_Cliques, 3),
   INTPARAM (CPXPARAM_MIP_Cuts_Covers, 3),
   INTPARAM (CPXPARAM_MIP_Cuts_Disjunctive, 3),
   INTPARAM (CPXPARAM_MIP_Cuts_FlowCovers, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_Gomory, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_GUBCovers, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_Implied, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_MIRCut, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_ZeroHalfCut, 2),
   INTPARAM (CPXPARAM_MIP_Cuts_MCFCut, 2),
   ENDPARAM
};

/** Predefined parameter settings. */
static struct {
   char const        *name;
   struct paramvalue const *values;
} const settings[] = {
   { "primal only", primalbound },
   { "dual bound",  dualbound }
};
#define NUMSETTINGS ((int)(sizeof (settings) / sizeof (settings[0])))

/** Apply a predefined parameter settings.
 * \param env     The environment to which the predefined setting is applied.
 * \param setting Index of the setting to apply.
 */
static void
applySettings (CPXENVptr env, int setting)
{
   int status = 0;
   struct paramvalue const *vals = settings[setting % NUMSETTINGS].values;

   status = CPXXsetintparam (env, CPXPARAM_RandomSeed, setting);
   if ( status ) {
      fprintf (stderr, "CPXXsetintparam(CPXPARAM_RandomSeed): %d\n", status);
      abort ();
   }

   while (vals->type != CPX_PARAMTYPE_NONE) {
      switch (vals->type) {
      case CPX_PARAMTYPE_INT:
         status = CPXXsetintparam (env, vals->num, vals->value.i);
         break;
      case CPX_PARAMTYPE_LONG:
         status = CPXXsetlongparam (env, vals->num, vals->value.l);
         break;
      case CPX_PARAMTYPE_DOUBLE:
         status = CPXXsetdblparam (env, vals->num, vals->value.d);
         break;
      }
      if ( status ) {
         fprintf (stderr, "CPXXset???param(%d): %d\n", vals->num, status);
         abort ();
      }
      ++vals;
   }
}


/** Install/reset the info callback at the remote end.
 */
static int installCallback (CPXENVptr env)
{
   return CPXXuserfunction (env, USERACTION_ADDCALLBACK, 0, NULL, 0, 0, NULL);
}

/** Remove the info callback at the remote end.
 */
static int removeCallback (CPXENVptr env)
{
   return CPXXuserfunction (env, USERACTION_REMOVECALLBACK, 0, NULL, 0, 0, NULL);
}

/** Change the tolerance that controls whether consecutive objective
 * function values are considered different.
 */
static int changeObjdiff (CPXENVptr env, double newdiff)
{
   int status;
   CPXSERIALIZERptr s = NULL;

   CPXXserializercreate (&s);
   s->adddouble (s, newdiff);
   status = CPXXuserfunction (env, USERACTION_CHANGEOBJDIFF,
                              CPXXserializerlength (s),
                              CPXXserializerpayload (s),
                              0, 0, NULL);
   CPXXserializerdestroy (s);

   return status;
}

/** Description of best known primal bound. */
static struct {
   int volatile     valid; /**< True if the bound/env fields are valid. */
   double volatile  bound; /**< Best known primal bound. */
   int volatile     idx;   /**< The environment that reported bound. */
} primal = { 0, 0.0, -1 };

/** Description of best known dual bound. */
static struct {
   int volatile     valid; /**< True if the bound/env fields are valid. */
   double volatile  bound; /**< Best known dual bound. */
   int volatile     idx;   /**< The environment that reported bound. */
} dual = { 0, 0.0, -1 };

static int objsen;

/** Current status of a remote worker. */
static struct remotestat {
   CPXENVptr env;
   double    dettime;
   double    primal;
   double    dual;
   int       idx;
} *remotestats = NULL;

/** This function is invoked whenever a remote worker reports new
 * information.
 * The information reported by remote workers are:
 * - new dual bounds,
 * - new primal bounds,
 * - current deterministic timestamps.
 * The handler functions updates the remotestats[] entry for the respective
 * remote worker and also updates that global best known primal/dual bound
 * if necessary.
 */
static void CPXPUBLIC
infohandler (CPXENVptr xenv, CPXINFOTYPE type, int tag, CPXLONG length,
             void const *data, void *handle)
{
   struct remotestat *const rs = handle;
   double d;
   (void)xenv;
   (void)type;
   (void)length;

   switch (tag) {
   case INFO_NEWDUAL:
      assert (type == CPXINFO_DOUBLE);
      assert (length == 1);
      d = *(double const *)data;
      fprintf (stderr, "[%d] New dual bound: %f\n", rs->idx, d);
      fflush (stderr);
      rs->dual = d;
      if ( !dual.valid ||
           (objsen == CPX_MIN && d > dual.bound) ||
           (objsen == CPX_MAX && d < dual.bound) )
      {
         dual.valid = 1;
         dual.bound = *(double const *)data;
         dual.idx = rs->idx;
      }
      break;
   case INFO_NEWPRIMAL:
      assert (type == CPXINFO_DOUBLE);
      assert (length == 1);
      d = *(double const *)data;
      fprintf (stderr, "[%d] New primal bound: %f\n", rs->idx, d);
      fflush (stderr);
      rs->primal = d;
      if ( !primal.valid ||
           (objsen == CPX_MIN && d < primal.bound) ||
           (objsen == CPX_MAX && d > primal.bound) )
      {
         primal.valid = 1;
         primal.bound = *(double const *)data;
         primal.idx = rs->idx;
      }
      break;
   case INFO_DETTIME:
      assert (type == CPXINFO_DOUBLE);
      assert (length == 1);
      d = *(double const *)data;
      fprintf (stderr, "[%d] Current dettime: %f (%f, %f)\n",
               rs->idx, d, rs->dual, rs->primal);
      fflush (stderr);
      rs->dettime = d;
      break;
   default:
      fprintf (stderr, "[%p] Unknown info %d\n", handle, tag);
      fflush (stderr);
   }
}

/** Function used for prefixed output.
 *  The function just prints \a message prefixed by the value of \a handle.
 */
static void CPXPUBLIC
printer (void *handle, char const *message)
{
   printf ("[%p] %s", handle, message);
}

int
main(int argc, char **argv)
{
   int status;
   CPXENVptr *env;
   CPXLPptr *lp;
   char const *modelfile = NULL;
   CPXASYNCptr *handle;
   CPXENVGROUPptr group;
   int i;
   int jobs = 0;
   int active;
   int *finished;
   char const **machine;
   char cwd[MAX_PATH_LEN];
   char usrfunc[MAX_PATH_LEN];
   int frequency;
   double absgap = 1e-6;
   int bestidx;
   CPXDIM c, cols;
   double *x;
   enum {
      OUTPUT_SILENT, OUTPUT_PREFIXED, OUTPUT_LOG
   } output = OUTPUT_SILENT;

#if defined(USE_MPI)
   int numprocs, rank;

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
   machine = malloc (sizeof (*machine) * numprocs);
   if ( machine == NULL ) {
      fprintf (stderr, "Out of memory!\n");
      abort ();
   }
   for (i = 0; i < numprocs; ++i)
      machine[i] = "mpimachine";
   jobs = numprocs - 1;
#elif defined(USE_PROCESS)
   char const *bin = "./cplex";

   machine = malloc (sizeof (*machine) * argc);
   if ( machine == NULL ) {
      fprintf (stderr, "Out of memory!\n");
      abort ();
   }
#elif defined(USE_TCPIP)
   machine = malloc (sizeof (*machine) * argc);
   if ( machine == NULL ) {
      fprintf (stderr, "Out of memory!\n");
      abort ();
   }
#else
#   error "No transport type selected"
#endif


   /* Parse the command line. */
   for (i = 1; i < argc; ++i) {
      if ( strncmp (argv[i], "-model=", 7) == 0 )
         modelfile = argv[i] + 7;
#if defined(USE_MPI)
#elif defined(USE_PROCESS)
      else if ( strncmp (argv[i], "-machine=", 9) == 0 )
         machine[jobs++] = argv[i] + 9;
      else if ( strncmp (argv[i], "-bin=", 5) == 0 )
         bin = argv[i] + 5;
#elif defined(USE_TCPIP)
      else if ( strncmp (argv[i], "-address=", 9) == 0 )
         machine[jobs++] = argv[i];
#endif
      else if ( strncmp (argv[i], "-absgap=", 8) == 0 )
         absgap = strtod (argv[i] + 8, NULL);
      else if ( strcmp (argv[i], "-output-prefixed") == 0 )
         output = OUTPUT_PREFIXED;
      else if ( strcmp (argv[i], "-output-log") == 0 )
         output = OUTPUT_LOG;
   }

   /* Validate arguments.
    */
   if ( modelfile == NULL ) {
      fprintf (stderr, "No model file specified with -model=<modelfile>\n");
      abort ();
   }
   if ( jobs < 1 ) {
     fprintf (stderr, "Invalid job count %d\n", jobs);
     abort ();
   }

   /* Allocate working arrays. */
   if ( (env = malloc (sizeof (*env) * jobs)) == NULL ||
        (handle = malloc (sizeof (*handle) * jobs)) == NULL ||
        (lp = malloc (sizeof (*lp) * jobs)) == NULL ||
        (finished = calloc (jobs, sizeof (*finished))) == NULL ||
        (remotestats = calloc (jobs, sizeof (*remotestats))) == NULL )
   {
      fprintf (stderr, "Out of memory!\n");
      abort ();
   }

   /* Find the place at which to find the shared object that implements
    * the user function. On Windows the path to the current directory is
    * likely to contain blanks, so better quote it.
    */
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

   /* Create a remote object instances. */
   for (i = 0; i < jobs; ++i) {
      /* These values define how we connect to the remote object. It is
       * important to use a transport configuration that actually supports
       * disconnect/reconnect. For the "processtransport" this means to use
       * named pipes instead of anonymous pipes or stdio.
       */
      char const *transport;
      char const *args[16];
      int nextarg = 0;
      char *logpath = NULL;

#if defined(USE_MPI)
      char rankbuf[256];

      sprintf (rankbuf, "-remoterank=%d", i + 1);
      transport = "mpitransport";
      args[nextarg++] = rankbuf;
#elif defined(USE_PROCESS)
      char logbuf[1024];
      transport = "processtransport";
      /* If the machine is not "localhost" then use ssh to connect to
       * this machine. Otherwise just fork a process on the local
       * machine.
       */
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
      args[nextarg++] = "-userfunction=parmipopt_userfunction=REGISTER_USERFUNCTION";
      sprintf (logbuf, "-logfile=server%d.log", i);
      if ( (args[nextarg] = logpath = strdup (logbuf)) != NULL )
         ++nextarg;
#elif defined(USE_TCPIP)
      transport = "tcpiptransport";
      args[nextarg++] = machine[i];
#endif


      printf ("Creating env on %s\n", machine[i]);
      env[i] = CPXXopenCPLEXremote(transport, nextarg, args, &status);
      if ( status || env[i] == NULL ) {
         fprintf (stderr, "CPXXopenCPLEXremote: %d\n", status);
         abort ();
      }
      free (logpath);

      /* Enable output */
      switch (output) {
      case OUTPUT_SILENT:
         /* nothing */
         break;
      case OUTPUT_PREFIXED:
         {
            CPXCHANNELptr cres, cwar, cerr, clog;

            if ( (status = CPXXgetchannels (env[i], &cres, &cwar, &cerr, &clog)) != 0 ) {
               fprintf (stderr, "CPXXgetchannels: %d\n", status);
               abort ();
            }
            if ( (status = CPXXaddfuncdest (env[i], cres, env[i], printer)) != 0 ||
                 (status = CPXXaddfuncdest (env[i], cwar, env[i], printer)) != 0 ||
                 (status = CPXXaddfuncdest (env[i], clog, env[i], printer)) != 0 )
            {
               fprintf (stderr, "CPXXaddfpdest: %d\n", status);
               abort ();
            }
         }
         break;
      case OUTPUT_LOG:
         {
            if ( (status = CPXXsetintparam (env[i], CPXPARAM_ScreenOutput, CPX_ON)) != 0 ) {
               fprintf (stderr, "CPXXgetchannels: %d\n", status);
               abort ();
            }
         }
         break;
      }

      /* Create empty problem object for this remote solver. */
      printf ("Creating LP %d\n", i);
      lp[i] = CPXXcreateprob(env[i], &status, "problem");
      if ( status || lp[i] == NULL ) {
         fprintf (stderr, "CPXXcreateprob: %d\n", status);
         abort ();
      }

      /* Install and configure callbacks. */
      remotestats[i].env = env[i];
      remotestats[i].idx = i;
      if ( (status = CPXXsetinfohandler (env[i], infohandler, &remotestats[i])) != 0 ) {
         fprintf (stderr, "CPXXsetinfohandler: %d\n", status);
         abort ();
      }

      if ( (status = changeObjdiff (env[i], 1e-5)) != 0 ) {
         fprintf (stderr, "changeObjdiff: %d\n", status);
         abort ();
      }

      if ( (status = installCallback (env[i])) != 0 ) {
         fprintf (stderr, "installCallback: %d\n", status);
         abort ();
      }

      /* Apply predefined perameter settings for this solver. */
      applySettings (env[i], i);
   }

   /* Put all environments into one group so that we can use multicasts
    * on operations that are the same for all solvers and/or imply lots
    * of data exchange.
    */
   status = CPXXcreateenvgroup (&group, jobs, env);
   if ( status != 0 ) {
      fprintf (stderr, "CPXXcreateenvgroup: %d\n", status);
      abort ();
   }

   /* Read the model into all remote solver. */
   status = CPXXreadcopyprob_multicast (group, modelfile, NULL);
   if ( status != 0 ) {
      fprintf (stderr, "CPXXreadcopyprob_multicast: %d\n", status);
      abort ();
   }
   objsen = CPXXgetobjsen (env[0], lp[0]);

   /* We set the thread count for each solver to 1 so that we do not
    * run into problems if multiple solves are performed on the same
    * machine.
    */
   status = CPXXsetintparam_multicast (group, CPXPARAM_Threads, 1);
   if ( status != 0 ) {
      fprintf (stderr, "CPXXsetintparam_multicast: %d\n", status);
      abort ();
   }

   /* Start an asynchronous solve on each remote solver. */
   for (i = 0; i < jobs; ++i) {
      printf ("Solving %d\n", i);
      if ( (status = CPXXmipopt_async (env[i], lp[i], &handle[i])) != 0 ) {
         fprintf (stderr, "CPXXmipopt_async: %d\n", status);
         abort ();
      }
   }

   /* All solves are started. Loop until the stopping criterion is met. */
   active = jobs;
   frequency = 10000; /* Print current bounds every two seconds. */
   while (active > 0) {
      int running = 0;
      /* Check if we shold stop all solves.
       * We stop them if the absolute mipgap is reached.
       */
      if ( primal.valid && dual.valid &&
           ((objsen == CPX_MIN && dual.bound + absgap >= primal.bound) ||
            (objsen == CPX_MAX && dual.bound - absgap <= primal.bound)) )
      {
         printf ("Stopping criterion reached. Stopping all pending solves.\n");
         for (i = 0; i < jobs; ++i) {
            if ( !finished[i] )
               CPXXasynckill (handle[i]);
         }
         break;
      }
      if ( --frequency == 0 ) {
         printf ("dual=%f, primal=%f\n", dual.bound, primal.bound);
         frequency = 10000;
      }

      /* Loop over all solvers and test if they are still running. */
      for (i = 0; i < jobs; ++i) {
         if ( finished[i] )
            continue;
         CPXXasynctest (handle[i], &running);
         if ( !running ) {
            /* The job is finished. We have a solution, so kill all
             * others. */
            int j;
            --active;
            finished[i] = 1;
            printf ("First job (%d) is finished, killing the rest\n", i);
            for (j = 0; j < jobs; ++j) {
               if ( j != i )
                  CPXXasynckill (handle[j]);
            }
            break;
         }
      }
      millisleep (10);
   }

   /* All solves have finished. Join them. */
   for (i = 0; i < jobs; ++i) {
      double obj = -CPX_INFBOUND;
      int stat;
      
      status = CPXXmipopt_join (&handle[i]);
      if ( status ) {
         fprintf (stderr, "CPXXmipopt_join: %d\n", status);
         abort ();
      }

      status = CPXXgetobjval (env[i], lp[i], &obj);
      if ( status == CPXERR_NO_SOLN ) {
         /* No feasible solution found (yet) on this machine.
          * Just set objective function to a very big value
          */
         obj = (objsen == CPX_MIN) ? CPX_INFBOUND : -CPX_INFBOUND;
      }
      else if ( status ) {
         fprintf (stderr, "CPXXgetobjval: %d\n", status);
         abort ();
      }
      stat = CPXXgetstat (env[i], lp[i]);
      printf ("Job %d: %f, stat %d\n", i, obj, stat);
      printf ("\t%f, %f, %f\n", remotestats[i].dettime,
              remotestats[i].dual, remotestats[i].primal);

      if ( (status = removeCallback (env[i])) != 0 ) {
         fprintf (stderr, "removeCallback: %d\n", status);
         abort ();
      }
   }

   /* Fetch the x vector from the solver that produced the best
    * primal bound. */
   bestidx = primal.idx;
   cols = CPXXgetnumcols (env[bestidx], lp[bestidx]);
   if ( (x = malloc (cols * sizeof (*x))) == NULL ) {
      fprintf (stderr, "Out of memory!\n");
      abort ();
   }
   status = CPXXgetx (env[bestidx], lp[bestidx], x, 0, cols - 1);
   if ( status ) {
      fprintf (stderr, "CPXXgetx: %d\n", status);
      abort ();
   }
   printf ("Optimal solution:\n");
   for (c = 0; c < cols; ++c)
      printf ("x[%5d]: %f\n", c, x[c]);
   free (x);

   CPXXfreeenvgroup (&group);

   /* Close the CPLEX objects in _reverse_ order. */
   for (i = jobs - 1; i >= 0; --i)
      CPXXcloseCPLEX (&env[i]);
   free (remotestats);
   free (env);
   free (lp);
   free (handle);
   free (finished);
   free ((char **)machine);

#ifdef USE_MPI
   MPI_Finalize ();
#endif


   return 0;
}

#endif /* COMPILE_MASTER */
