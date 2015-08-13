/* --------------------------------------------------------------------------
 * File: parbenders.c
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#if defined(USE_MPI)
#   include <mpi.h>
#   define TRANSPORT "mpitransport"
#elif defined(USE_PROCESS)
#   define TRANSPORT "processtransport"
#elif defined(USE_TCPIP)
#   define TRANSPORT "tcpiptransport"
#else
#   error "No transport type selected"
#endif


/* ********************************************************************** *
 *                                                                        *
 *    C l i e n t   s i d e   i m p l e m e n t a t i o n                 *
 *                                                                        *
 * ********************************************************************** */

#if defined(COMPILE_MASTER)

#include <ilcplex/cplexremotemasterx.h>

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
static int
init (CPXENVptr env, CPXLPptr *lp_p, CPXDIM** block_p, CPXDIM* nblocks_p)
{
   /* Model data.
    * fixed[] is the fixed cost for opening a facility,
    * cost[i,j] is the cost for serving customer i from facility j.
    */
   static double const fixed[] = { 2.0, 3.0, 3.0 };
   static double const cost[] = { 2.0, 3.0, 4.0, 5.0, 7.0,
                                  4.0, 3.0, 1.0, 2.0, 6.0,
                                  5.0, 4.0, 2.0, 1.0, 3.0 };
#define NFACTORY ((CPXDIM)(sizeof(fixed) / sizeof(fixed[0])))
#define NCUSTOMER ((CPXDIM)((sizeof(cost) / sizeof(cost[0])) / NFACTORY))

   CPXDIM ind[NFACTORY < 2 ? 2 : NFACTORY];
   double val[NFACTORY < 2 ? 2 : NFACTORY];
   double rhs;
   CPXNNZ const rmatbeg = 0;
   int status = 0;
   CPXLPptr lp = NULL;
   CPXDIM *block= NULL;
   char namebuf[32];
   char const *name[1] = { namebuf };
   CPXDIM j, f, c;

   *lp_p = NULL;
   *block_p = NULL;
   *nblocks_p = -1;

   lp = CPXXcreateprob (env, &status, "problem");
   if ( lp == NULL || status != 0 )
      goto TERMINATE;

   block = malloc (sizeof (*cost) * (NFACTORY + NFACTORY * NCUSTOMER));
   if ( block == NULL) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }

   /* Create integer y  variables. */
   j = 0;
   for (f = 0; f < NFACTORY; ++f) {
      double const lb = 0.0;
      double const ub = 1.0;

      sprintf (namebuf, "y%d", f);
      status = CPXXnewcols (env, lp, 1, &fixed[f], &lb, &ub, "B", name);
      if ( status != 0 )
         goto TERMINATE;
      block[j++] = -1;
   }

   /* Create continuous x variables. */
   for (f = 0; f < NFACTORY; ++f) {
      for (c = 0; c < NCUSTOMER; ++c) {
         double const lb = 0.0;
         double const ub = CPX_INFBOUND;

         sprintf (namebuf, "x%d#%d", f, c);
         status = CPXXnewcols (env, lp, 1, &cost[f * NCUSTOMER + c],
                               &lb, &ub, "C", name);
         if ( status != 0 )
            goto TERMINATE;
         block[j++] = c;
      }
   }

   /* Satisfy each customer's demand. */
   for (c = 0; c < NCUSTOMER; ++c) {
      rhs = 1.0;
      for (f = 0; f < NFACTORY; ++f) {
         val[f] = 1.0;
         ind[f] = NFACTORY + f * NCUSTOMER + c;
      }
      sprintf (namebuf, "c1#%d", c);
      status = CPXXaddrows (env, lp, 0, 1, NFACTORY, &rhs, "G",
                            &rmatbeg, ind, val, NULL, name);
      if ( status != 0 )
         goto TERMINATE;
   }

   /* A factory must be open if we service from it. */
   for (c = 0; c < NCUSTOMER; ++c) {
      for (f = 0; f < NFACTORY; ++f) {
         rhs = 0.0;
         ind[0] = NFACTORY + f * NCUSTOMER + c;
         val[0] = -1.0;
         ind[1] = f;
         val[1] = 1.0;
         sprintf (namebuf, "c2#%d#%d", c, f);
         status = CPXXaddrows (env, lp, 0, 1, 2, &rhs, "G", &rmatbeg,
                               ind, val, NULL, name);
         if ( status != 0 )
            goto TERMINATE;
      }
   }

   /* Capacity constraint. */
   rhs = NFACTORY - 1;
   for (f = 0; f < NFACTORY; ++f) {
      ind[f] = f;
      val[f] = 1.0;
   }
   status = CPXXaddrows (env, lp, 0, 1, NFACTORY, &rhs, "L", &rmatbeg, ind,
                         val, NULL, NULL);
   if ( status != 0 )
      goto TERMINATE;

   /* Setup return value and clear arrays so that they are not freed. */
   *lp_p = lp;
   *block_p = block;
   *nblocks_p = NCUSTOMER;
   lp = NULL;
   block = NULL;
 TERMINATE:
   free (block);
   CPXXfreeprob (env, &lp);

   return status;

#undef NFACTORY
#undef NCUSTOMER
}

/** Descriptor for a block into which the problem decomposes.
 * This structure is used to represent the master as well as the sub-blocks.
 */
typedef struct {
   CPXENVptr   env;      /* Environment that stores LP. */
   CPXLPptr    lp;       /* The block described by this instance. */
   CPXDIM      number;   /* The block's serial number. */
   struct {
      CPXNNZ    n;       /* Length of the arrays below. */
      double    *coef;   /* Coefficient of variable. */
      CPXDIM    *col;    /* Column index of fixed variable in master. */
      CPXDIM    *row;    /* Row index of fixed variable in block (column
                          * index in block's dual).
                          */
   }           fixed;    /* Description of variables that appear in this
                          * block but are supposed to be fixed by solves
                          * of the master block. */
   double      objshift; /* Objective shift from CPXXdualwrite(). */

   double      *objval;  /* Worker array for cut separation. */
   CPXDIM      *objind;  /* Worker array for cut separation. */
   double      *actobj;  /* Worker array for cut separation. */

   CPXASYNCptr async;    /* Handle for asynchronous calls. */

   CPXDIM      *map;     /* Map indices in the original problem to
                          * indices in the master.
                          */
} BLOCK;

/** Delete an instance of BLOCK.
 */
static void
free_block (BLOCK **block_p)
{
   if ( block_p && *block_p ) {
      BLOCK *b = *block_p;

      free (b->map);
      free (b->actobj);
      free (b->objind);
      free (b->objval);
      free (b->fixed.row);
      free (b->fixed.col);
      free (b->fixed.coef);
      CPXXfreeprob (b->env, &b->lp);
      CPXXcloseCPLEX (&b->env);
      free (b);
      *block_p = NULL;
   }
}

/* Extract block NUMBER from LP and store it in *BLOCK_P.
 * Note that the functions stores the dual of the BLOCK problem.
 * That dual problem is created on a remote object instance. This instance
 * is created using ARGC and ARGV using the "process" transport.
 */
static int
extract_block (CPXCENVptr env, CPXCLPptr lp, CPXDIM const *block,
               CPXDIM number, BLOCK **block_p,
               int argc, char const *const *argv, char const *const *machines)
{
   int status = 0;
   BLOCK *b = NULL;
   CPXDIM const ocols = CPXXgetnumcols (env, lp);
   CPXDIM const orows = CPXXgetnumrows (env, lp);
   CPXDIM const maxind = (ocols < orows) ? orows : ocols;
   CPXDIM i, j, cols;
   unsigned char *mark = NULL;
   CPXDIM *ind = NULL;
   double *val = NULL;
   CPXDIM *map = NULL;
   CPXLPptr primal = NULL;
   CPXNNZ nfixed = 0;
   CPXNNZ fixcap = 0;
   CPXDIM *fixcol = NULL;
   CPXDIM *fixrow = NULL;
   double *fixcoef = NULL;
   char const **transargv = NULL;
   char extraargs[512];
   char *extra = extraargs;

   printf ("Extracting block %d ...\n", number);
   fflush (stdout);

   *block_p = NULL;

   /* Allocate working arrays. */
   if ( (transargv = malloc ((argc + 3) * sizeof (*transargv))) == NULL ||
        (b = calloc (1, sizeof (*b))) == NULL ||
        (mark = calloc (orows, sizeof (*mark))) == NULL ||
        (ind = malloc (maxind * sizeof (*ind))) == NULL ||
        (val = malloc (maxind * sizeof (*val))) == NULL ||
        (map = malloc (ocols * sizeof (*map))) == NULL )
   {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }

   /* Add a logfile argument to the transport arguments. */
   memcpy (transargv, argv, argc * (sizeof (*argv)));

#if defined(USE_MPI)
   sprintf (extra, "-remoterank=%d", number + 1);
   transargv[argc++] = extra;
   (void)machines;
#elif defined(USE_PROCESS)
   sprintf (extra, "-logfile=block%04d.log", number);
   transargv[argc++] = extra;
   (void)machines;
#elif defined(USE_TCPIP)
   transargv[argc++] = machines[number];
#endif

   /* Create remote object and an empty problem instance on that. */
   b->env = CPXXopenCPLEXremote (TRANSPORT, argc, transargv, &status);
   if ( status != 0 )
      goto TERMINATE;
   b->lp = CPXXcreateprob (b->env, &status, "blockdual");
   if ( status != 0 )
      goto TERMINATE;
   primal = CPXXcreateprob (env, &status, "blockprimal");
   if ( status != 0 )
      goto TERMINATE;

   /* Copy non-fixed variables from original problem into block. */
   for (j = 0, cols = 0; j < ocols; ++j) {
      if ( block[j] == number ) {
         CPXNNZ nnz, surplus, cmatbeg, n;
         double lb, ub, obj;
         char ctype;
         char namestore[4096];
         char *name[1];
         CPXSIZE nsurplus;

         /* Query column information and create new column in block LP
          * with exactly the same information.
          */
         if ( (status = CPXXgetlb (env, lp, &lb, j, j)) != 0 ||
              (status = CPXXgetub (env, lp, &ub, j, j)) != 0 ||
              (status = CPXXgetobj (env, lp, &obj, j, j)) != 0 ||
              (status = CPXXgetctype (env, lp, &ctype, j, j)) != 0 ||
              (status = CPXXgetcolname (env, lp, name, namestore,
                                        sizeof (namestore), &nsurplus,
                                        j, j) != 0) )
            goto TERMINATE;
         if ( ctype != 'C' ) {
            fprintf (stderr,
                     "Cannot handle non-continuous variable %d (%c) in block %d.\n",
                     j, ctype, number);
            status = CPXERR_BAD_ARGUMENT;
            goto TERMINATE;
         }
         /* Normalize objective function to "minimize". */
         if ( CPXXgetobjsen (env, lp) != CPX_MIN )
            obj *= -1.0;
         status = CPXXnewcols (env, primal, 1, &obj, &lb,&ub, NULL,
                               (char const *const *)name);
         if ( status != 0 )
            goto TERMINATE;

         /* Mark the rows that intersect this column so that we can pick
          * them up later.
          */
         status = CPXXgetcols (env, lp, &nnz, &cmatbeg, ind, val,
                               orows, &surplus, j, j);
         if ( status != 0 )
            goto TERMINATE;
         for (n = 0; n < nnz; ++n)
            mark[ind[n]] = 1;

         /* Record the index that the copied variable has in the
          * block model. */
         map[j] = cols++;
      }
      else
         map[j] = -1;
   }

   /* Now copy all rows that intersect block variables. */
   for (i = 0; i < orows; ++i) {
      if ( mark[i] ) {
         CPXNNZ nnz, surplus, n, rmatbeg;
         CPXDIM rowlen;
         double rhs;
         char sense;
         double factor;

         /* Query the row we want to copy. */
         if ( (status = CPXXgetrhs (env, lp, &rhs, i, i)) != 0 ||
              (status = CPXXgetsense (env, lp, &sense, i, i)) != 0 ||
              (status = CPXXgetrows (env, lp, &nnz, &rmatbeg, ind, val,
                                     ocols, &surplus, i, i)) != 0 )
            goto TERMINATE;

         /* Fix up variable indices and normalize constraint to '<='.
          * Normalization is technically not required but makes the models
          * a little easier to read.
          */
         switch (sense) {
         case 'E': factor = 1.0; break;
         case 'L': factor = 1.0; break;
         case 'G': factor = -1.0; sense = 'L'; break;
         default:
            fprintf (stderr, "Unsupported constraint sense '%c'\n", sense);
            status = CPXERR_BAD_ARGUMENT;
            goto TERMINATE;
         }
         rhs *= factor;

         /* Copy the row. */
         for (n = 0, rowlen = 0; n < nnz;++n) {
            val[n] = factor * val[n];
            j = ind[n];
            if ( block[j] != number ) {
               /* This column is not explicitly in this block. This means
                * that it is a column that will be fixed by the master.
                * We collect all such columns so that we can adjust the
                * dual objective function according to concrete fixings.
                */
               if ( nfixed >= fixcap ) {
                  CPXDIM *newcol = NULL, *newrow = NULL;
                  double *newcoef = NULL;

                  /* Make sure the arrays to store information about
                   * fixed variables are large enough.
                   */
                  if ( fixcap == 0 ) {
                     fixcap = 16;
                     if ( (newcol = malloc (sizeof (*newcol) * fixcap)) == NULL ||
                          (newrow = malloc (sizeof (*newrow) * fixcap)) == NULL ||
                          (newcoef = malloc (sizeof (*newcoef) * fixcap)) ==NULL )
                     {
                        free (newrow);
                        free (newcol);
                        status = CPXERR_NO_MEMORY;
                        goto TERMINATE;
                     }
                  }
                  else {
                     fixcap *= 2;
                     if ( (newcol = realloc (fixcol, sizeof (*newcol) * fixcap)) == NULL ||
                          (newrow = realloc (fixrow, sizeof (*newrow) * fixcap)) == NULL ||
                          (newcoef = realloc (fixcoef, sizeof (*newcoef) * fixcap)) == NULL )
                     {
                        free (newrow);
                        free (newcol);
                        status = CPXERR_NO_MEMORY;
                        goto TERMINATE;
                     }
                  }
                  fixrow = newrow;
                  fixcol = newcol;
                  fixcoef = newcoef;
               }

               /* Store information about variables in this block that
                * will be fixed by master solves.
                */
               fixrow[nfixed] = CPXXgetnumrows (env, primal);
               fixcol[nfixed] = j;
               fixcoef[nfixed] = -val[n];
               ++nfixed;
            }
            else {
               /* The column is an ordinary in this block. Just copy it. */
               ind[rowlen] = map[j];
               val[rowlen] = val[n];
               ++rowlen;
            }
         }

         status = CPXXaddrows (env, primal, 0, 1, rowlen, &rhs, &sense,
                               &rmatbeg, ind, val, NULL, NULL);
         if ( status != 0 )
            goto TERMINATE;
      }
   }

   /* Create the dual of the block model we just created and store that
    * in the block descriptor. To keep the code simple we just use a CPLEX
    * function that writes out a dual problem and read that problem back in.
    */
   if ( (status = CPXXdualwrite (env, primal, "dual.mps", &b->objshift)) != 0 ||
        (status = CPXXreadcopyprob (b->env, b->lp, "dual.mps", NULL)) != 0 )
      goto TERMINATE;
   /* dualwrite may write the dual as 'min'. Fix that up. */
   if ( CPXXgetobjsen (b->env, b->lp) != CPX_MAX ) {
      b->objshift *= -1.0;
      cols = CPXXgetnumcols (b->env, b->lp);
      CPXXgetobj (b->env, b->lp, val, 0, cols - 1);
      for (j = 0; j < cols; ++j) {
         val[j] = -val[j];
         ind[j] =j;
      }
      CPXXchgobj (b->env, b->lp, cols, ind, val);
      CPXXchgobjsen (b->env, b->lp, CPX_MAX);
   }

   /* Setup the arrays that we need in order to fiddle with the objective
    * function.
    */
   cols = CPXXgetnumcols (b->env, b->lp);
   if ( (b->objval = malloc (cols * sizeof (*b->objval))) == NULL ||
        (b->objind = malloc (cols * sizeof (*b->objind))) == NULL ||
        (b->actobj = malloc (cols * sizeof (*b->actobj))) == NULL )
   {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }
   if ( (status = CPXXgetobj (b->env, b->lp, b->objval, 0, cols - 1)) != 0 )
      goto TERMINATE;
   for (j = 0; j < cols; ++j)
      b->objind[j] = j;

   /* Setup return value.
    */
   b->fixed.n = nfixed;
   b->fixed.coef = fixcoef;
   b->fixed.col = fixcol;
   b->fixed.row = fixrow;
   *block_p = b;

   /* Set pointers to NULL so that they do not get deleted. */
   fixcoef = NULL;
   fixcol = NULL;
   fixrow = NULL;
   b = NULL;

 TERMINATE:
   free (fixcoef);
   free (fixrow);
   free (fixcol);
   CPXXfreeprob (env, &primal);
   free (map);
   free (ind);
   free (val);
   free (mark);
   free_block (&b);
   free ((char **)transargv);

   printf ("block %d %s (%d).\n", number, status ? "failed" : "ok", status);

   return status;
}

/** Extract the master block from LP.
 * The function also fixes up references to variables in the master
 * for all the NBLOCKS in BLOCKS.
 * The function returns the master block in *BLOCK_P and the index of
 * the first eta variable in *FIRSTETA_P.
 */
static int
extract_master (CPXCENVptr env, CPXCLPptr lp, CPXDIM const *block,
                BLOCK **blocks, CPXDIM nblocks, BLOCK **block_p,
                CPXDIM *firsteta_p)
{
   int status = 0;
   BLOCK *b = NULL;
   unsigned char *mark = NULL;
   CPXDIM const rows = CPXXgetnumrows (env, lp);
   CPXDIM const cols = CPXXgetnumcols (env, lp);
   CPXDIM const maxind = (rows < cols) ? cols : rows;
   CPXDIM *ind = NULL;
   double *val = NULL;
   CPXDIM *map = NULL;
   CPXDIM i, j, newcols;

   *block_p = NULL;

   /* Allocate worker arrays and return value. */
   if ( (b = calloc (1, sizeof (*b))) == NULL ||
        (mark = calloc (rows, sizeof (*mark))) == NULL ||
        (ind = malloc (maxind * sizeof (*ind))) == NULL ||
        (val = malloc (maxind * sizeof (*val))) == NULL ||
        (map = malloc (cols * sizeof (*map))) == NULL )
      goto TERMINATE;

   /* Create an empty master model.
    * This model will be solved on the local machine, so we don't
    * allocate the model on a remote object.
    */
   b->env = CPXXopenCPLEX (&status);
   if ( b->env == NULL || status != 0 )
      goto TERMINATE;
   b->lp = CPXXcreateprob (b->env, &status, "master");
   if ( b->lp == NULL || status != 0 )
      goto TERMINATE;

   /* Find columns that do not intersect block variables and
    * copy them to the master block.
    */
   for (j = 0, newcols = 0; j < cols; ++j) {
      if ( block[j] < 0 ) {
         /* Column is not in a block. Copy it to the master. */
         double lb, ub, obj;
         char ctype;
         char namebuf[4096];
         char *name[1];
         CPXSIZE surplus;

         if ( (status = CPXXgetlb (env, lp, &lb, j, j)) != 0 ||
              (status = CPXXgetub (env, lp, &ub, j, j)) != 0 ||
              (status = CPXXgetobj (env, lp, &obj, j, j)) != 0 ||
              (status = CPXXgetctype (env, lp, &ctype, j, j)) != 0 ||
              (status = CPXXgetcolname (env, lp, name, namebuf,
                                        sizeof (namebuf), &surplus,
                                        j, j)) != 0 )
            goto TERMINATE;

         status = CPXXnewcols (b->env, b->lp, 1, &obj, &lb, &ub, &ctype,
                               (char const *const *)name);
         if ( status != 0 )
            goto TERMINATE;

         /* Record index of newly created variable in the master model. */
         map[j] = newcols++;
      }
      else {
         /* Column is in a block. Mark all rows that intersect
          * this column.
          */
         CPXNNZ nnz, surplus, cmatbeg, n;

         status = CPXXgetcols (env, lp, &nnz, &cmatbeg, ind, val, maxind,
                               &surplus, j, j);
         if ( status )
            goto TERMINATE;
         for (n = 0; n < nnz; ++n)
            mark[ind[n]] = 1;
         map[j] = -1;
      }
   }

   /* Pick up the rows that we need to copy.
    * These are the rows that are only intersected by master variables,
    * that is, the rows that are not marked.
    */
   for (i = 0; i < rows; ++i) {
      if ( !mark[i] ) {
         CPXNNZ nnz, surplus, rmatbeg, n;
         double rhs;
         char sense;
         char namebuf[4096];
         char *name[1];
         CPXSIZE surplus2;

         /* Get the row. */
         if ( (status = CPXXgetrhs (env, lp, &rhs, i, i)) != 0 ||
              (status = CPXXgetsense (env, lp, &sense, i, i)) != 0 ||
              (status = CPXXgetrows (env, lp, &nnz, &rmatbeg, ind, val,
                                     maxind, &surplus, i, i)) != 0 ||
              (status = CPXXgetrowname (env, lp, name, namebuf,
                                        sizeof (namebuf), &surplus2,
                                        i, i)) != 0 )
            goto TERMINATE;
         /* Map indices in the original problem to indices in the master.*/
         for (n = 0; n < nnz; ++n)
            ind[n] = map[ind[n]];
         /* Add the row just read to the master model. */
         status = CPXXaddrows (b->env, b->lp, 0, 1, nnz, &rhs, &sense,
                               &rmatbeg, ind, val, NULL,
                               (char const *const *)name);
         if ( status != 0 )
            goto TERMINATE;
      }
   }

   /* Adjust variable indices in blocks so that reference to variables
    * in the original problem become references to variables in the master.
    */
   for (i = 0; i < nblocks; ++i) {
      for (j = 0; j < blocks[i]->fixed.n; ++j)
         blocks[i]->fixed.col[j] = map[blocks[i]->fixed.col[j]];
   }

   /* Create the eta variables, one for each block.
    * See the comments at the top of this file for details about the
    * eta variables.
    */
   *firsteta_p = CPXXgetnumcols (b->env, b->lp);
   for (i = 0; i < nblocks; ++i) {
      double const obj = 1.0, lb = 0.0, ub = CPX_INFBOUND;
      char namebuf[128];
      char const *name[1] = { namebuf };

      sprintf (namebuf, "_eta%d", i);
      if ( (status = CPXXnewcols (b->env, b->lp, 1, &obj, &lb, &ub, "C",
                                  name)) != 0 )
         goto TERMINATE;
   }

   b->map = map;
   map = NULL;
   *block_p = b;
   b = NULL;
 TERMINATE:
   free_block (&b);
   free (map);
   free (val);
   free (ind);
   free (mark);

   return status;
}

/** Data used in callbacks to separate Benders' cuts.
 */
typedef struct {
   BLOCK  *master;     /* master block. */
   BLOCK  **blocks;    /* sub blocks. */
   CPXDIM nblocks;     /* Length of BLOCKS. */
   CPXDIM *cutind;     /* Worker array for cut separation. */
   double *cutval;     /* Worker array for cut separation. */
   double *x;          /* Worker array to query solution from master. */
   CPXDIM etaind;      /* Index of first eta variable in master. */
} CBHANDLE;

/** Lazy constraint callback that separates Benders cuts.
 * The callback is invoked whenever we find an integer feasible solution
 * for the master block.
 */
static int CPXPUBLIC
callback (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle,
          int *useraction_p)
{
   CBHANDLE *handle = cbhandle;
   int status;
   CPXDIM const mcols = CPXXgetnumcols (handle->master->env,
                                        handle->master->lp);
   CPXDIM b;
   double *tmp = NULL, *ray = NULL;

   printf ("Callback invoked. Separate Benders cuts.\n");

   if ( (tmp = calloc (mcols, sizeof (*tmp))) == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }

   /* Get the current relaxation. */
   status = CPXXgetcallbacknodex (env, cbdata, wherefrom, handle->x,
                                  0, mcols - 1);
   if ( status )
      goto TERMINATE;

   /* Iterate over blocks and trigger a separation on each of them.
    * The separation is triggered asynchronously so that it can happen
    * on different remote objects simultaneously.
    */
   for (b = 0; b < handle->nblocks; ++b) {
      BLOCK *block = handle->blocks[b];
      CPXDIM const ncols = CPXXgetnumcols (block->env, block->lp);
      CPXDIM j;

      /* Iterate over the fixed master variables in this block.
       * Each fixed variable goes to the right-hand side and therefore
       * into the objective function.
       */
      memcpy (block->actobj, block->objval, sizeof (*block->objval) * ncols);
      for (j = 0; j < block->fixed.n; ++j) {
         block->actobj[block->fixed.row[j]] -= (block->fixed.coef[j] *
                                                handle->x[block->fixed.col[j]]);
      }
      status = CPXXchgobj (block->env, block->lp, ncols, block->objind,
                           block->actobj);
      if ( status != 0 ) {
         fprintf (stderr, "CHGOBJ failed (%d)\n", status);
         goto TERMINATE;
      }

      /* If the problem is unbounded we need to get an infinite ray in
       * order to be able to generate the respective Benders cut. If
       * CPLEX proves unboundedness in presolve then it will return
       * CPX_STAT_INForUNBD and no ray will be available. So we need to
       * disable presolve.
       */
      CPXXsetintparam (block->env, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
      CPXXsetintparam (block->env, CPXPARAM_Preprocessing_Reduce, 0);

      /* Solve the updated problem to optimality. */
      status = CPXXprimopt_async (block->env, block->lp, &block->async);
      if ( status != 0 ) {
         /* Optimization failed. Note that we need to join all other
          * parallel solves (ignoring their results) before we exit
          * the callback. Otherwise we will have pending asynchronous
          * calls which result in resource leaks.
          */
         fprintf (stderr, "PRIMPOPT failed\n");
         while (--b >= 0)
            (void)CPXXprimopt_join (&handle->blocks[b]->async);
         goto TERMINATE;
      }
   }

   /* Wait for the various LP solves to complete. */
   for (b = 0; b < handle->nblocks; ++b) {
      BLOCK *block = handle->blocks[b];
      int s;

      s = CPXXprimopt_join (&block->async);
      if ( s != 0 && status == 0 )
         status = s;
   }

   if ( status )
      goto TERMINATE;

   /* See if we need to generate cuts. */
   for (b = 0; b < handle->nblocks; ++b) {
      BLOCK *block = handle->blocks[b];
      CPXDIM const ncols = CPXXgetnumcols (block->env, block->lp);
      double *val = handle->cutval;
      CPXDIM *ind = handle->cutind;
      CPXDIM j;
      CPXDIM cutnz = 0;
      char cutsense;
      double cutrhs;

      if ( (ray = malloc (ncols * sizeof (*ray))) == NULL ) {
         status = CPXERR_NO_MEMORY;
         goto TERMINATE;
      }
      memset (tmp, 0, mcols * sizeof (*tmp));

      /* Depending on the status either seperate a feasibility or an
       * optimality cut.
       */
      switch (CPXXgetstat (block->env, block->lp)) {
      case CPX_STAT_UNBOUNDED:
         {
            /* The subproblem is unbounded. We need to extract a feasibility
             * cut from an unbounded ray of the problem (see also the comments
             * at the top of this file).
             */
            printf ("unbounded ");
            if ( (status = CPXXgetray (block->env, block->lp, ray)) != 0 )
               goto TERMINATE;
            cutrhs = 0;
            for (j = 0; j < ncols; ++j)
               cutrhs -= ray[j] * block->objval[j];
            cutnz = 0;
            for (j = 0; j < block->fixed.n; ++j) {
               tmp[block->fixed.col[j]] -= (block->fixed.coef[j] *
                                            ray[block->fixed.row[j]]);
            }
            for (j = 0; j < mcols; ++j) {
               if ( fabs (tmp[j]) > 1e-6 ) {
                  val[cutnz] = tmp[j];
                  ind[cutnz] = j;
                  ++cutnz;
               }
            }
            cutsense = 'L';
         }
         break;
      case CPX_STAT_OPTIMAL:
         {
            /* The subproblem has a finite optimal solution.
             * We need to check if this gives rise to an optimality cut (see
             * also the comments at the top of this file).
             */
            CPXDIM const etaind = handle->etaind + b;
            double objval;
            double eta = handle->x[etaind];
            
            printf ("optimal ");
            if ( (status = CPXXgetobjval (block->env, block->lp,
                                          &objval)) != 0 ||
                 (status = CPXXgetx (block->env, block->lp, ray, 0,
                                     ncols - 1)) != 0 )
               goto TERMINATE;

            objval += block->objshift;
            
            if ( objval > eta + 1e-6 ) {
               cutrhs = 0;
               for (j = 0; j < ncols; ++j)
                  cutrhs -= ray[j] * block->objval[j];
               cutnz = 0;
               for (j = 0; j < block->fixed.n; ++j) {
                  tmp[block->fixed.col[j]] -= (block->fixed.coef[j] *
                                               ray[block->fixed.row[j]]);
               }
               for (j = 0; j < mcols; ++j) {
                  if ( fabs (tmp[j]) > 1e-6 ) {
                     val[cutnz] = tmp[j];
                     ind[cutnz] = j;
                     ++cutnz;
                  }
               }
               ind[cutnz] = etaind;
               val[cutnz] = -1.0;
               ++cutnz;
               cutsense = 'L';
            }
         }
         break;
      default:
         fprintf (stderr, "Unexpected status %d\n",
                  CPXXgetstat (block->env, block->lp));
         status = -1;
         break;
      }
      
      /* If a cut was found then add that. */
      if ( cutnz > 0 ) {
         printf ("cut found (%lld):", (long long)cutnz);
         for (j = 0; j < cutnz; ++j)
            printf (" %+fx[%d]", handle->cutval[j], handle->cutind[j]);
         printf (" %c %f\n", cutsense, cutrhs);
         status = CPXXcutcallbackadd (env, cbdata, wherefrom, cutnz, cutrhs,
                                      cutsense, handle->cutind, handle->cutval,
                                      CPX_USECUT_FORCE);
         if ( status != 0 )
            goto TERMINATE;
      }
      else
         printf ("no cuts.\n");
      free (ray);
      ray = NULL;
   }

 TERMINATE:
   free (tmp);
   free (ray);

   return status;
}

/** Solve the problem in LP using a distributed implementation of Benders'
 * decomposition.
 * The function decomposes the problem in LP according to NBLOCKS and BLOCKS
 * and then performs a distributed Benders' decomposition solve. For each
 * sub-block it will create a remote object so that Benders' cuts can be
 * separated in parallel on different machines. To create the remote objects
 * the function will use the arguments in ARGC and ARGV.
 *
 * ENV      --  The environment that holds LP.
 * LP       --  The problem to be solved.
 * NBLOCKS  --  Number of blocks in BLOCK.
 * BLOCK    --  An array with length equal to the number of columns in LP.
 *              For a column j BLOCK[j] is negative if j is in the master or
 *              specifies the non-negative block number (in [0,NBLOCKS-1])
 *              for the block that contains j.
 * ARGC     --  Length of ARGV.
 * ARGV     --  Connection arguments for the process transport.
 *
 * ATTENTION: To get the final results the function will perform a final
 *            solve with all integer variables in the master fixed. So
 *            upon return the bounds in LP may have been changed.
 */
static int
bendersopt (CPXENVptr env, CPXLPptr lp, CPXDIM nblocks, CPXDIM const *block,
            int argc, char const *const *argv, char const *const *machines)
{
   int status;
   BLOCK **blocks = NULL;
   BLOCK *master = NULL;
   CBHANDLE handle = { NULL, NULL, 0, NULL, NULL, NULL, -1 };
   CPXDIM b;
   CPXDIM const cols = CPXXgetnumcols (env, lp);
   CPXDIM j, mcols;
   double *x = NULL;
   double *origlb = NULL, *origub = NULL, *lb = NULL, *ub = NULL;
   CPXDIM *ind = NULL;
   char *lu = NULL;
   int restorelb = 0, restoreub = 0;

   /* Extract blocks and master problem. */
   printf ("Extracting %d blocks.\n", nblocks);
   if ( (blocks = calloc (nblocks, sizeof (*blocks))) == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }
   for (b = 0; b < nblocks; ++b) {
      if ( (status = extract_block (env, lp, block, b, &blocks[b],
                                    argc, argv, machines)) != 0 )
         goto TERMINATE;
   }
   if ( (status = extract_master (env, lp, block, blocks, nblocks,
                                  &master, &handle.etaind)) != 0 )
      goto TERMINATE;

   /* Write out the master and all blocks (for debugging). */
   CPXXwriteprob (master->env, master->lp, "master.lp", NULL);
   for (b = 0; b < nblocks; ++b) {
      char name[256];
      sprintf (name, "block%03d.lp", b);
      CPXXwriteprob (blocks[b]->env, blocks[b]->lp, name, NULL);
   }

   /* Setup callback that separates Benders' cuts. */
   handle.master = master;
   handle.blocks = blocks;
   handle.nblocks = nblocks;
   if ( (handle.cutind = malloc (sizeof (*handle.cutind) * CPXXgetnumcols (master->env, master->lp))) == NULL ||
        (handle.cutval = malloc (sizeof (*handle.cutval) * CPXXgetnumcols (master->env, master->lp))) == NULL ||
        (handle.x = malloc (sizeof (*handle.x) * CPXXgetnumcols (master->env, master->lp))) == NULL )
   {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }

   status = CPXXsetintparam (master->env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
   if ( status != 0 )
      goto TERMINATE;
   status = CPXXsetlazyconstraintcallbackfunc (master->env, callback, &handle);
   if ( status != 0 )
      goto TERMINATE;

   /* Solve the problem. */
   status = CPXXmipopt (master->env, master->lp);
   if ( status != 0 )
      goto TERMINATE;

   mcols = CPXXgetnumcols (master->env, master->lp);
   if ( (x = malloc (mcols * sizeof (*x))) == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }
   if ( (status = CPXXgetx (master->env, master->lp, x, 0, mcols - 1)) != 0 )
      goto TERMINATE;

   /* Perform a last solve with the integral variables in the master
    * fixed to get solution values into the original problem instance.
    * An alternative way to get the solution into the original problem
    * object would be to add the solution in x[] as a MIP start and solve
    * the original problem with a solution of 1.
    */
   if ( (origlb = malloc (cols * sizeof (*origlb))) == NULL ||
        (origub = malloc (cols * sizeof (*origub))) == NULL ||
        (lb = malloc (cols * sizeof (*lb))) == NULL ||
        (ub = malloc (cols * sizeof (*ub))) == NULL ||
        (ind = malloc (cols * sizeof (*ind))) == NULL ||
        (lu = malloc (cols * sizeof (*lu))) == NULL )
   {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }

   if ( (status = CPXXgetlb (env, lp, origlb, 0, cols - 1)) != 0 ||
        (status = CPXXgetub (env, lp, origub, 0, cols - 1)) != 0 ||
        (status = CPXXgetctype (env, lp, lu, 0, cols - 1)) != 0 )
      goto TERMINATE;

   for (j = 0; j < cols; ++j) {
      ind[j] = j;
      if ( block[j] < 0 && (lu[j] == 'B' || lu[j] == 'I') ) {
         lb[j] = floor (x[master->map[j]] + 0.5);
         ub[j] = floor (x[master->map[j]] + 0.5);
      }
      else {
         lb[j] = origlb[j];
         ub[j] = origub[j];
      }
   }

   memset (lu, 'L', cols * sizeof (*lu));
   if ( (status = CPXXchgbds (env, lp, cols, ind, lu, lb)) != 0 )
      goto TERMINATE;
   restorelb = 1;

   memset (lu, 'U', cols * sizeof (*lu));
   if ( (status = CPXXchgbds (env, lp, cols, ind, lu, ub)) != 0 )
      goto TERMINATE;
   restoreub = 1;

   if ( (status = CPXXmipopt (env, lp)) != 0 )
      goto TERMINATE;

 TERMINATE:
   free (lu);
   free (ind);
   free (origlb);
   free (origub);
   free (lb);
   free (ub);
   free (x);
   free (handle.x);
   free (handle.cutval);
   free (handle.cutind);
   free_block(&master);
   if ( blocks != NULL ) {
      for (b = 0; b < nblocks; ++b)
         free_block (&blocks[b]);
      free (blocks);
   }
   return status;
}

int
main (int argc, char **argv)
{
   CPXENVptr env = NULL;
   CPXLPptr lp = NULL;
   CPXDIM *block = NULL;
   CPXDIM nblocks = -1;
   int status;
   int myargc, i;
   char const **myargv = NULL;
   double objval;
   CPXDIM cols;
   CPXDIM j;
   char const **machines = NULL;
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
   if ( (machines = malloc (argc * sizeof (*machines))) == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }
   nmachines = 0;
   myargc = 0;
#endif


   /* Process the command line. */
   if ( (myargv = malloc (argc * sizeof (myargv))) == NULL ) {
      status = CPXERR_NO_MEMORY;
      goto TERMINATE;
   }

   for (i = 1; i < argc; ++i) {
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
      if ( strncmp (argv[i], "-address=", 9) == 0 )
         machines[nmachines++] = argv[i];
#endif
      else
         myargv[myargc++] = argv[i];
   }

   /* Initialize the input problem. */
   env = CPXXopenCPLEX (&status);
   if ( env == NULL || status != 0 )
      goto TERMINATE;
   (void)CPXXsetintparam (env, CPXPARAM_ScreenOutput, CPX_ON);
   status = init (env, &lp, &block, &nblocks);
   if ( status != 0 )
      goto TERMINATE;
   (void)CPXXwriteprob (env, lp, "benders.lp", NULL);

   if ( nblocks > nmachines ) {
      fprintf (stderr, "We have %d blocks but only %d machines!\n",
               nblocks, nmachines);
      status = CPXERR_BAD_ARGUMENT;
      goto TERMINATE;
   }

   status = bendersopt (env, lp, nblocks, block, myargc, myargv, machines);
   if ( status != 0 )
      goto TERMINATE;

   cols = CPXXgetnumcols (env, lp);
   if ( (status = CPXXgetobjval (env, lp, &objval)) != 0 )
      goto TERMINATE;
   printf ("#### Problem solved (%f, %d).\n", objval,
           CPXXgetstat (env, lp));
   for (j = 0; j < cols; ++j) {
      double x;

      if ( (status = CPXXgetx (env, lp, &x, j, j)) != 0 )
         goto TERMINATE;
      printf ("#### \tx[%d] = %f\n", j, x);
   }

 TERMINATE:
   free (block);
   CPXXfreeprob (env, &lp);
   CPXXcloseCPLEX (&env);
   free ((char **)myargv);

#if defined(USE_MPI)
   MPI_Finalize ();
#elif defined(USE_PROCESS)
   free (machines);
#elif defined(USE_TCPIP)
   free (machines);
#endif


   return status;
}

#endif /* COMPILE_MASTER */
