/* --------------------------------------------------------------------------
 * File: xlpex5.c
 * Version 12.6.1
 * --------------------------------------------------------------------------
 * Licensed Materials - Property of IBM
 * 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
 * Copyright IBM Corporation 1997, 2014. All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with
 * IBM Corp.
 * --------------------------------------------------------------------------
 */

/* xlpex5.c - Illustrating the CPLEX message handler
   This is a modification of lpex1.c .
   We'll label each channels output with its name, and
   also put our own output into a file.  Note that the labels
   may cause the output to go past 80 characters per line. */

/* Bring in the CPLEX function declarations and the C library 
   header file stdio.h with the following single include. */

#include <ilcplex/cplexx.h>

/* Bring in the declarations for the string functions */

#include <string.h>

/* Include declaration for function at end of program */

static int
   populatebycolumn  (CPXENVptr env, CPXLPptr lp);

static void CPXPUBLIC
   ourmsgfunc     (void *handle, const char *message);



/* The problem we are optimizing will have 2 rows, 3 columns 
   and 6 nonzeros.  */

#define NUMROWS    2
#define NUMCOLS    3
#define NUMNZ      6


int
main (void)
{
   char     probname[16];  /* Problem name is max 16 characters */

   /* Declare and allocate space for the variables and arrays where we
      will store the optimization results including the status, objective
      value, variable values, dual values, row slacks and variable
      reduced costs. */

   int      solstat;
   double   objval;
   double   x[NUMCOLS];
   double   pi[NUMROWS];
   double   slack[NUMROWS];
   double   dj[NUMCOLS];


   CPXENVptr     env = NULL;
   CPXLPptr      lp = NULL;
   int           status;
   CPXDIM        i, j;
   CPXDIM        cur_numrows, cur_numcols;
   char          errmsg[CPXMESSAGEBUFSIZE];

   CPXCHANNELptr  cpxerror   = NULL;
   CPXCHANNELptr  cpxwarning = NULL;
   CPXCHANNELptr  cpxresults = NULL;
   CPXCHANNELptr  ourchannel = NULL;

   char *errorlabel = "cpxerror";
   char *warnlabel  = "cpxwarning";
   char *reslabel   = "cpxresults";
   char *ourlabel   = "Our Channel";

   CPXFILEptr fpout  = NULL;


   /* Initialize the CPLEX environment */

   env = CPXXopenCPLEX (&status);

   /* If an error occurs, the status value indicates the reason for
      failure.  A call to CPXXgeterrorstring will produce the text of
      the error message.  Note that CPXXopenCPLEX produces no output,
      so the only way to see the cause of the error is to use
      CPXXgeterrorstring.  For other CPLEX routines, the errors will
      be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON.  */

   /* Since the message handler is yet to be set up, we'll call our
      messaging function directly to print out any errors  */

   if ( env == NULL ) {
      ourmsgfunc ("Our Message", "Could not open CPLEX environment.\n");
      goto TERMINATE;
   }

   /* Now get the standard channels.  If an error, just call our
      message function directly. */

   status = CPXXgetchannels (env, &cpxresults, &cpxwarning, &cpxerror, NULL);
   if ( status ) {
      ourmsgfunc ("Our Message", "Could not get standard channels.\n");
      CPXXgeterrorstring (env, status, errmsg);
      ourmsgfunc ("Our Message", errmsg);
      goto TERMINATE;
   }

   /* Now set up the error channel first.  The label will be "cpxerror" */

   status = CPXXaddfuncdest (env, cpxerror, errorlabel, ourmsgfunc);
   if ( status ) {
      ourmsgfunc ("Our Message", "Could not set up error message handler.\n");
      CPXXgeterrorstring (env, status, errmsg);
      ourmsgfunc ("Our Message", errmsg);
   }

   /* Now that we have the error message handler set up, all CPLEX
      generated errors will go through ourmsgfunc.  So we don't have
      to use CPXXgeterrorstring to determine the text of the message.
      We can also use CPXXmsg to do any other printing.  */

   status = CPXXaddfuncdest (env, cpxwarning, warnlabel, ourmsgfunc);
   if ( status ) {
      CPXXmsg (cpxerror, "Failed to set up handler for cpxwarning.\n");
      goto TERMINATE;
   }

   status = CPXXaddfuncdest (env, cpxresults, reslabel, ourmsgfunc);
   if ( status ) {
      CPXXmsg (cpxerror, "Failed to set up handler for cpxresults.\n");
      goto TERMINATE;
   }
   
   /* Now turn on the iteration display. */

   status = CPXXsetintparam (env, CPXPARAM_Simplex_Display, 2);
   if ( status ) {
      CPXXmsg (cpxerror, "Failed to turn on simplex display level.\n");
      goto TERMINATE;
   }

   /* Create the problem. */

   strcpy (probname, "example");
   lp = CPXXcreateprob (env, &status, probname);

   /* A returned pointer of NULL may mean that not enough memory
      was available or there was some other problem.  In the case of 
      failure, an error message will have been written to the error 
      channel from inside CPLEX.  In this example, the setting of
      the parameter CPXPARAM_ScreenOutput causes the error message to
      appear on stdout.  */

   if ( lp == NULL ) {
      CPXXmsg (cpxerror, "Failed to create LP.\n");
      goto TERMINATE;
   }

   /* Now populate the problem with the data. */

   status = populatebycolumn (env, lp);

   if ( status ) {
      CPXXmsg (cpxerror, "Failed to populate problem data.\n");
      goto TERMINATE;
   }


   /* Optimize the problem and obtain solution. */

   status = CPXXlpopt (env, lp);
   if ( status ) {
      CPXXmsg (cpxerror, "Failed to optimize LP.\n");
      goto TERMINATE;
   }

   status = CPXXsolution (env, lp, &solstat, &objval, x, pi, slack, dj);
   if ( status ) {
      CPXXmsg (cpxerror, "Failed to obtain solution.\n");
      goto TERMINATE;
   }


   /* Write the output to the screen.  We will also write it to a
      file as well by setting up a file destination and a function
      destination. */

   ourchannel = CPXXaddchannel (env);
   if ( ourchannel == NULL ) {
      CPXXmsg (cpxerror, "Failed to set up our private channel.\n");
      goto TERMINATE;
   }

   fpout = CPXXfopen ("lpex5.msg", "w");
   if ( fpout == NULL ) {
      CPXXmsg (cpxerror, "Failed to open lpex5.msg file for output.\n");
      goto TERMINATE;
   }
   status = CPXXaddfpdest (env, ourchannel, fpout);
   if ( status ) {
      CPXXmsg (cpxerror, "Failed to set up output file destination.\n");
      goto TERMINATE;
   }

   status = CPXXaddfuncdest (env, ourchannel, ourlabel, ourmsgfunc);
   if ( status ) {
      CPXXmsg (cpxerror, "Failed to set up our output function.\n");
      goto TERMINATE;
   }

   /* Now any message to channel ourchannel will go into the file 
      and into the file opened above. */

   CPXXmsg (ourchannel, "\nSolution status = %d\n", solstat);
   CPXXmsg (ourchannel, "Solution value  = %f\n\n", objval);

   /* The size of the problem should be obtained by asking CPLEX what
      the actual size is, rather than using sizes from when the problem
      was built.  cur_numrows and cur_numcols store the current number 
      of rows and columns, respectively.  */

   cur_numrows = CPXXgetnumrows (env, lp);
   cur_numcols = CPXXgetnumcols (env, lp);
   for (i = 0; i < cur_numrows; i++) {
      CPXXmsg (ourchannel, "Row %d:  Slack = %10f  Pi = %10f\n", 
              i, slack[i], pi[i]);
   }

   for (j = 0; j < cur_numcols; j++) {
      CPXXmsg (ourchannel, "Column %d:  Value = %10f  Reduced cost = %10f\n",
              j, x[j], dj[j]);
   }

   /* Finally, write a copy of the problem to a file. */

   status = CPXXwriteprob (env, lp, "lpex5.lp", NULL);
   if ( status ) {
      CPXXmsg (cpxerror, "Failed to write LP to disk.\n");
      goto TERMINATE;
   }
   
   
TERMINATE:

   /* First check if ourchannel is open */

   if ( ourchannel != NULL ) {
      int  chanstat;
      chanstat = CPXXdelfuncdest (env, ourchannel, ourlabel, ourmsgfunc);
      if ( chanstat ) {
         strcpy (errmsg, "CPXXdelfuncdest failed.\n");
         ourmsgfunc ("Our Message", errmsg); 
         if (!status)  status = chanstat;
      }
      if ( fpout != NULL ) {
         chanstat = CPXXdelfpdest (env, ourchannel, fpout);
         if ( chanstat ) {
            strcpy (errmsg, "CPXXdelfpdest failed.\n");
            ourmsgfunc ("Our Message", errmsg);
            if (!status)  status = chanstat;
         }
         CPXXfclose (fpout);
      }

      chanstat = CPXXdelchannel (env, &ourchannel);
      if ( chanstat ) {
         strcpy (errmsg, "CPXXdelchannel failed.\n");
         ourmsgfunc ("Our Message", errmsg);
         if (!status)  status = chanstat;
      }
   }

   /* Free up the problem as allocated by CPXXcreateprob, if necessary */

   if ( lp != NULL ) {
      status = CPXXfreeprob (env, &lp);
      if ( status ) {
         strcpy (errmsg, "CPXXfreeprob failed.\n");
         ourmsgfunc ("Our Message", errmsg);
      }
   }

   /* Now delete our function destinations from the 3 CPLEX channels. */
   if ( cpxresults != NULL ) {
      int  chanstat;
      chanstat = CPXXdelfuncdest (env, cpxresults, reslabel, ourmsgfunc);
      if ( chanstat && !status ) {
         status = chanstat;
         strcpy (errmsg, "Failed to delete cpxresults function.\n");
         ourmsgfunc ("Our Message", errmsg);
      }
   }

   if ( cpxwarning != NULL ) {
      int  chanstat;
      chanstat = CPXXdelfuncdest (env, cpxwarning, warnlabel, ourmsgfunc);
      if ( chanstat && !status ) {
         status = chanstat;
         strcpy (errmsg, "Failed to delete cpxwarning function.\n");
         ourmsgfunc ("Our Message", errmsg);
      }
   }

   if ( cpxerror != NULL ) {
      int  chanstat;
      chanstat = CPXXdelfuncdest (env, cpxerror, errorlabel, ourmsgfunc);
      if ( chanstat && !status ) {
         status = chanstat;
         strcpy (errmsg, "Failed to delete cpxerror function.\n");
         ourmsgfunc ("Our Message", errmsg);
      }
   }

   /* Free up the CPLEX environment, if necessary */

   if ( env != NULL ) {
      status = CPXXcloseCPLEX (&env);

      /* Note that CPXXcloseCPLEX produces no output,
         so the only way to see the cause of the error is to use
         CPXXgeterrorstring.  For other CPLEX routines, the errors will
         be seen if the CPXPARAM_ScreenOutput indicator is set to CPX_ON. */

      if ( status ) {
         strcpy (errmsg, "Could not close CPLEX environment.\n");
         ourmsgfunc ("Our Message", errmsg);
         CPXXgeterrorstring (env, status, errmsg);
         ourmsgfunc ("Our Message", errmsg);
      }
   }
     
   return (status);

}  /* END main */


/* This function builds by column the linear program:

      Maximize
       obj: x1 + 2 x2 + 3 x3
      Subject To
       c1: - x1 + x2 + x3 <= 20
       c2: x1 - 3 x2 + x3 <= 30
      Bounds
       0 <= x1 <= 40
      End
 */

static int
populatebycolumn (CPXENVptr env, CPXLPptr lp)
{
   int        status    = 0;
   double     obj[NUMCOLS];
   double     lb[NUMCOLS];
   double     ub[NUMCOLS];
   char const *colname[NUMCOLS];
   CPXNNZ     matbeg[NUMCOLS];
   CPXDIM     matind[NUMNZ];
   double     matval[NUMNZ];
   double     rhs[NUMROWS];
   char       sense[NUMROWS];
   char const *rowname[NUMROWS];

   /* To build the problem by column, create the rows, and then 
      add the columns. */

   status = CPXXchgobjsen (env, lp, CPX_MAX);  /* Problem is maximization */
   if ( status )  goto TERMINATE;

   /* Now create the new rows.  First, populate the arrays. */

   rowname[0] = "c1";
   sense[0]   = 'L';
   rhs[0]     = 20.0;

   rowname[1] = "c2";
   sense[1]   = 'L';
   rhs[1]     = 30.0;

   status = CPXXnewrows (env, lp, NUMROWS, rhs, sense, NULL, rowname);
   if ( status )   goto TERMINATE;

   /* Now add the new columns.  First, populate the arrays. */

       obj[0] = 1.0;      obj[1] = 2.0;           obj[2] = 3.0;

    matbeg[0] = 0;     matbeg[1] = 2;          matbeg[2] = 4;
      
    matind[0] = 0;     matind[2] = 0;          matind[4] = 0;
    matval[0] = -1.0;  matval[2] = 1.0;        matval[4] = 1.0;
 
    matind[1] = 1;     matind[3] = 1;          matind[5] = 1;
    matval[1] = 1.0;   matval[3] = -3.0;       matval[5] = 1.0;

        lb[0] = 0.0;       lb[1] = 0.0;           lb[2]  = 0.0;
        ub[0] = 40.0;      ub[1] = CPX_INFBOUND;  ub[2]  = CPX_INFBOUND;

   colname[0] = "x1"; colname[1] = "x2";      colname[2] = "x3";

   status = CPXXaddcols (env, lp, NUMCOLS, NUMNZ, obj, matbeg, matind,
                        matval, lb, ub, colname);
   if ( status )  goto TERMINATE;

TERMINATE:

   return (status);

}  /* END populatebycolumn */


/* For our message functions, we will interpret the handle as a pointer
 * to a string, which will be the label for the channel.  We'll put
 * angle brackets <> around the message so its clear what the function is
 * sending to us.  We'll place the newlines that appear at the end of
 * a message after the > bracket.
 */

static void CPXPUBLIC
ourmsgfunc (void *handle, const char *msg)
{
   char const *label;
   size_t     lenstr;
   int        flag = 0;
   char       message[CPXMESSAGEBUFSIZE + 2];

   /* Copy message and remove trailing '\n' if any. */
   strcpy(message, msg);
   lenstr = strlen(message);
   if ( message[lenstr-1] == '\n' ) {
      message[lenstr-1] = '\0';
      flag = 1;
   }

   label = handle;
   printf ("%-15s: <%s>", label, message);
   if (flag) putchar('\n');

} /* END ourmsgfunc */
