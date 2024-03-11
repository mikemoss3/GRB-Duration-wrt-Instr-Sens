#include <fitsio.h>
#include <math.h>
#include <string.h>
#include "pil.h"
#include "headas.h"
#include "headas_gti.h"
#include "battblocks.h"

/* 
 * battblocks - task for computing various time blocks from event or
 * light curve data
 * 
 * Main HEADAS driver routines
 *
 * C. Markwardt
 *
 * $Id$
 *
 */

#define TOOLSUB battblocks
static char taskname[] = "battblocks";
static char taskver[]  = "1.18";

/* headas_main() requires that TOOLSUB be defined first */
#include "headas_main.c"

/* Offset fudge time used when applying breakfile */
#define TFUDGE (1.0e-7)

/* ============================================================= */
struct parm_struct 
{
  char *taskname;
  char *taskver;
  char infile[PIL_PATH_MAX];    /* Light curve / event file name */
  char outfile[PIL_PATH_MAX];   /* Output GTI file name */
  char durfile[PIL_PATH_MAX];   /* Output of burst duration GTI file */
  char diagfile[PIL_PATH_MAX];  /* Output diagnostic file name */
  int nspill;
  double timedel;               /* Internal tick size for event data */
  double tlookback;             /* Lookback time for local BB analysis */
  double ncp_prior;             /* Prior probability for cost */
  double txx;                   /* User Txx percentage interval 0-100 */
  double tpeak;                 /* Window to find peak [s] */
  int gaussian;                 /* Use gaussian stats? 1=yes, 0=no,-1=unknown*/

  int intype;                   /* Input file type: EVENTS or LIGHTCURVE */
  int bkgsub;                   /* Background subtract light curve? 1=yes, 0=no */

  char timecol[80];             /* Input "TIME" column */
  char countscol[80];           /* Input COUNTS/RATE column */
  char errcol[80];              /* Input error column corresponding to countscol */
  char expocol[80];             /* Input exposure/timedel column (if irregularly sampled) */

  char hduclas3[80];            /* default HDUCLAS3 keyword value */
  double coalescefrac;           /* Coalesce first/last blocks if they are small enough */
  int durerrmeth;               /* Method to compute duration uncertainty:
				   1: total variance; 2: fractional variance */
  double global_tstart,         /* Global start/stop time of data in sec, */
         global_tstop;          /*   (default -1e307,+1e307) */

  double burst_tstart,          /* Burst start/stop time, to override BB algorithm,  */
         burst_tstop;           /*  when determining burst duration */
				   
  char breakfile[PIL_PATH_MAX]; /* Requested time-breaks */
};



/* ============================================================= */
/* Utility routine to read 'r' real parameters which could be INDEF */
int PILGetRealIndef(char *parname, double *data, double dfault)
{
  int status;
  status = PILGetReal(parname, data);
  if (status == PIL_VALUE_UNDEFINED) {
    *data = dfault;
    status = 0;
  }
  return status;
}

/* ============================================================= */
int battblocks_getpar(struct parm_struct *parms) 
{
  int status = 0;
  char gausstr[100];
  char durerrstr[100];

  parms->infile[0] = 0;
  parms->outfile[0] = 0;
  parms->diagfile[0] = 0;
  parms->breakfile[0] = 0;
  parms->gaussian = -1;
  parms->tlookback = 0;
  parms->bkgsub = 0;
  parms->global_tstart = -(FLOAT_INDEF);
  parms->global_tstop  = +(FLOAT_INDEF);
  parms->burst_tstart  = FLOAT_INDEF;
  parms->burst_tstop   = FLOAT_INDEF;

  strcpy(parms->timecol, "TIME");
  strcpy(parms->countscol, "INDEF");
  strcpy(parms->errcol,  "ERROR");
  strcpy(parms->expocol, "TIMEDEL");
  strcpy(parms->hduclas3, "INDEF");

  if ((status = PILGetFname("infile", parms->infile)))
    fprintf(stderr, "Error reading the 'infile' parameter.\n");

  else if ((status = PILGetString("outfile", parms->outfile)))
    fprintf(stderr, "Error reading the 'outfile' parameter.\n");

  else if ((status = PILGetString("durfile", parms->durfile)))
    fprintf(stderr, "Error reading the 'durfile' parameter.\n");

  else if ((status = PILGetInt("nspill", &parms->nspill)))
    fprintf(stderr, "Error reading the 'nspill' parameter.\n");

  else if ((status = PILGetReal("timedel", &parms->timedel)))
    fprintf(stderr, "Error reading the 'timedel' parameter.\n");

  else if ((status = PILGetReal("tlookback", &parms->tlookback)))
    fprintf(stderr, "Error reading the 'tlookback' parameter.\n");

  else if ((status = PILGetReal("tpeak", &parms->tpeak)))
    fprintf(stderr, "Error reading the 'tpeak' parameter.\n");

  else if ((status = PILGetReal("txx", &parms->txx)))
    fprintf(stderr, "Error reading the 'txx' parameter.\n");

  else if ((status = PILGetReal("ncp_prior", &parms->ncp_prior)))
    fprintf(stderr, "Error reading the 'ncp_prior' parameter.\n");

  else if ((status = PILGetString("gaussian", gausstr)))
    fprintf(stderr, "Error reading the 'gaussian' parameter.\n");
  
  else if ((status = PILGetBool("bkgsub", &parms->bkgsub)))
    fprintf(stderr, "Error reading the 'bkgsub' parameter.\n");

  else if ((status = PILGetReal("coalescefrac", &parms->coalescefrac)))
    fprintf(stderr, "Error reading the 'coalescefrac' parameter.\n");

  else if ((status = PILGetString("timecol", parms->timecol)))
    fprintf(stderr, "Error reading the 'timecol' parameter.\n");

  else if ((status = PILGetString("countscol", parms->countscol)))
    fprintf(stderr, "Error reading the 'countscol' parameter.\n");

  else if ((status = PILGetString("errcol", parms->errcol)))
    fprintf(stderr, "Error reading the 'errcol' parameter.\n");

  else if ((status = PILGetString("expocol", parms->expocol)))
    fprintf(stderr, "Error reading the 'expocol' parameter.\n");

  else if ((status = PILGetString("hduclas3", parms->hduclas3)))
    fprintf(stderr, "Error reading the 'hduclas3' parameter.\n");

  else if ((status = PILGetString("diagfile", parms->diagfile)))
    fprintf(stderr, "Error reading the 'diagfile' parameter.\n");
  
  else if ((status = PILGetString("durerrmeth", durerrstr)))
    fprintf(stderr, "Error reading the 'durerrmeth' parameter.\n");

  else if ((status = PILGetRealIndef("burst_tstart", &parms->burst_tstart, FLOAT_INDEF)))
    fprintf(stderr, "Error reading the 'burst_tstart' parameter.\n");

  else if ((status = PILGetRealIndef("burst_tstop", &parms->burst_tstop, FLOAT_INDEF)))
    fprintf(stderr, "Error reading the 'burst_tstop' parameter.\n");

  else if ((status = PILGetRealIndef("global_tstart", &parms->global_tstart, -FLOAT_INDEF)))
    fprintf(stderr, "Error reading the 'global_tstart' parameter.\n");

  else if ((status = PILGetRealIndef("global_tstop", &parms->global_tstop, +FLOAT_INDEF)))
    fprintf(stderr, "Error reading the 'global_tstop' parameter.\n");

  else if ((status = PILGetString("breakfile", parms->breakfile)))
    fprintf(stderr, "Error reading the 'breakfile' parameter.\n");

  if (status) return status;

  /* Determine if the user wants Gaussian stats: 'yes', 'no' or 'INDEF' */
  if ((gausstr[0] == 'y') || (gausstr[0] == 'Y')) {
    parms->gaussian = 1;
  } else if ((gausstr[0] == 'n') || (gausstr[0] == 'N')) {
    parms->gaussian = 0;
  } else if (strcasecmp(gausstr, "INDEF") == 0) {
    parms->gaussian = -1;
  } else {
    fprintf(stderr, "ERROR: gaussian must be yes, no, or INDEF\n");
    return -1;
  }

  /* Default value for diagfile */
  if (strcasecmp(parms->diagfile, "NONE") == 0) {
    parms->diagfile[0] = 0;
  }

  /* Default value for durfile */
  if (strcasecmp(parms->durfile, "NONE") == 0) {
    parms->durfile[0] = 0;
  }

  /* Default value for breakfile */
  if (strcasecmp(parms->breakfile, "NONE") == 0) {
    parms->breakfile[0] = 0;
  }

  /* Error checking on tlookback */
  if (parms->tlookback < 0) {
    fprintf(stderr, "WARNING: tlookback must be greater than or equal to 0\n");
    fprintf(stderr, "         (resetting to 0)\n");
    parms->tlookback = 0;
  }

  /* Default values for column names */
  if (strcasecmp(parms->timecol, "INDEF") == 0) {
    strcpy(parms->timecol, "TIME");
  }
  if (strcasecmp(parms->errcol, "INDEF") == 0) {
    strcpy(parms->timecol, "ERROR");
  }
  /* NOTE: countscol and expocol need to be determined on the fly */

  /* Convert durerrmeth from string to integer */
  if (strcasecmp(durerrstr, "TOTVAR") == 0) {
    parms->durerrmeth = 1;
  } else if (strcasecmp(durerrstr, "FRACVAR") == 0) {
    parms->durerrmeth = 2;
  } else {
    fprintf(stderr, "ERROR: 'durerrmeth' must be either TOTVAR or FRACVAR\n");
    return -1;
  }

  return status;
}

/* ============================================================= */
void banner(struct parm_struct *parms)
{
  headas_chat(2, "******************************************\n");
  headas_chat(1, "         %s v%s\n", parms->taskname, parms->taskver);
  headas_chat(2, "------------------------------------------\n");
  headas_chat(2, "     Input Data: %s\n", parms->infile);
  headas_chat(2, "     Output GTI: %s\n", parms->outfile);
  headas_chat(2, " Events to Skip: %d     Changepoint log(prob) Prior: %f\n", 
	      parms->nspill, parms->ncp_prior);
  headas_chat(2, "  Internal Tick: %e s\n", parms->timedel);
  headas_chat(2, "  Lookback time: %g s\n", parms->tlookback);
  headas_chat(2, " Bkg. Subtract?: %s (for fluence/T50/T90 calculations)\n", parms->bkgsub ? "yes" : "no");
  headas_chat(2, "------------------------------------------\n");
}

/* ============================================================= */
void summary(struct parm_struct *parms, int ngti)
{
  headas_chat(2, "  Created GTI with %d entries\n", ngti);
  headas_chat(2, "------------------------------------------\n");
}  

/* ============================================================= */
void printcps(int ncp, int *cparray)
{
  int i;

  headas_chat(5, " ...ncp=%d;  sample of values...\n", ncp);
  for(i=0; i<ncp; i++){
    headas_chat(5, "     %d %d\n", i, cparray[i]);
    if (i >= 10) {
      headas_chat(5, "   (...continues...)\n");
      break;
    }
  }
  if (ncp > 10) {
    headas_chat(5, "     %d %d\n", ncp-1, cparray[ncp-1]);
  }
}

/* ============================================================= */
/* Read light curve data from input file:
 *  t is input variable (number of times)
 */
int readlc(fitsfile *fitsptr, struct parm_struct *parms, 
	   double *t, double *rate, double *rerror, double *dt,
	   long int *nrows, char *hduclas3)
{
  double timedeli;
  int ratecol, errcol, dtcol;
  int i, j;
  int status = 0;
  char ratecolname[FLEN_CARD] = {0};
  char hduclas30[FLEN_CARD] = {0};
  double timepixr = 0.5;
  int typecode = 0;
  long int repeat = 0, width = 0;

  if (strcasecmp(parms->countscol, "INDEF") == 0) {
    strcpy(hduclas30, "RATE");
    strcpy(ratecolname, "RATE");
    fits_get_colnum(fitsptr, CASEINSEN, ratecolname, &ratecol, &status);
    if (status) {
      status = 0;
      strcpy(hduclas30, "COUNT");
      strcpy(ratecolname, "COUNTS");
      fits_get_colnum(fitsptr, CASEINSEN, ratecolname, &ratecol, &status);
    }
    if (status) {
      fprintf(stderr, "ERROR: could not locate RATE/COUNTS column in %s\n", 
	      parms->infile);
      return status;
    }
  } else {
    strcpy(hduclas30, "INDEF");
    strcpy(ratecolname, parms->countscol);
    fits_get_colnum(fitsptr, CASEINSEN, ratecolname, &ratecol, &status);
    if (status) {
      fprintf(stderr, "ERROR: could not locate %s column in %s\n", 
	      ratecolname, parms->infile);
    }
  }
  fits_get_coltype(fitsptr, ratecol, &typecode, &repeat, &width, &status);
  if (status) {
    fprintf(stderr, "ERROR: could not get type of %s column\n", ratecolname);
    return status;
  }
  if (repeat != 1) {
    fprintf(stderr, "ERROR: light curve must be a single channel OGIP light curve\n");
    return -1;
  }
  
  fits_read_col(fitsptr, TDOUBLE, ratecol, 1, 1, *nrows, 0, rate, 0, &status);
  if (status) {
    fprintf(stderr, "ERROR: could not read %s column data\n", ratecolname);
    return status;
  }

  if (strcasecmp(parms->expocol, "INDEF") == 0) {
    fits_get_colnum(fitsptr, CASEINSEN, "TIMEDEL", &dtcol, &status);
    if (status) {
      status = 0;
      fits_get_colnum(fitsptr, CASEINSEN, "EXPOSURE", &dtcol, &status);
    }
  } else {
    fits_get_colnum(fitsptr, CASEINSEN, parms->expocol, &dtcol, &status);
  }

  if (status == 0) {
    fits_read_col(fitsptr, TDOUBLE, dtcol, 1, 1, *nrows, 0, dt, 0, &status);
    if (status) {
      fprintf(stderr, "ERROR: could not read TIMEDEL/EXPOSURE data\n");
      return status;
    }
  } else {
    status = 0;
    fits_read_key(fitsptr, TDOUBLE, "TIMEDEL", &timedeli, 0, &status);
    if (status) {
      status = 0;
      fits_read_key(fitsptr, TDOUBLE, "EXPOSURE", &timedeli, 0, &status);
      timedeli /= *nrows;
    }
    if (status) {
      fprintf(stderr, "ERROR: could not find TIMEDEL/EXPOSURE data\n");
      return status;
    }
    
    for (i=0; i<(*nrows); i++) {
      dt[i] = timedeli;
    }
  }

  
  if (strcasecmp(parms->hduclas3,"INDEF") != 0) {
    strcpy(hduclas3, parms->hduclas3);
  } else {
    fits_read_key(fitsptr, TSTRING, "HDUCLAS3", hduclas3, 0, &status);
  }

  if (status) {
    
    /* We could not find the HDUCLAS3 keyword.  Try using the default
       value. */
    status = 0;
    strcpy(hduclas3, hduclas30);
    if ((hduclas3[0] == 0) || (strcasecmp(hduclas3, "INDEF") == 0)) {
      fprintf(stderr, "ERROR: HDUCLAS3 not found, there was no way to determine the type of this file\n");
      return status;
    }
    fprintf(stderr, "WARNING: HDUCLAS3 not found, assuming %s\n", hduclas3);

    /* We found the HDUCLAS3 keyword, but it does not match the
       expected default value. */
  } else if ((strcasecmp(hduclas3, hduclas30) != 0) &&
	     (strcasecmp(hduclas30, "INDEF") != 0)) {
    fprintf(stderr, "WARNING: HDUCLAS3 does not match expected RATE/COUNTS column (found %s, expected %s).\n", hduclas3, hduclas30);

    /* We found the HDUCLAS3 keyword, but it says RATE when the column
       name is COUNTS; or vice versa.  This is not an error, but it is
       suspicious. */
  } else if (((strcasecmp(hduclas3, "RATE") == 0) && 
	      strstr(ratecolname,"COUNT"))
	     ||
	     ((strcasecmp(hduclas3, "COUNT") == 0) &&
	      strstr(ratecolname,"RATE"))) {
    fprintf(stderr, "WARNING: HDUCLAS3 is %s and the light curve column name is %s.\n",
	    hduclas3, ratecolname);
    fprintf(stderr, "         These do not appear to match.  Is this right?\n");
  }

  headas_chat(5,"  ...HDUCLAS3=%s\n", hduclas3);

  fits_read_key(fitsptr, TDOUBLE, "TIMEPIXR", &timepixr, 0, &status);
  if (status) {
    status = 0;
    fprintf(stderr, "WARNING: TIMEPIXR not found, assuming %g\n", timepixr);
  }
  headas_chat(5,"  ...timepixr=%f...\n", timepixr);

  /* Translate time to center of the bin */
  for (i=0; i<(*nrows); i++) {
    t[i] = t[i] + (0.5-timepixr)*dt[i];
  }

  /* Read error column */
  if (parms->gaussian) {
    fits_get_colnum(fitsptr, CASEINSEN, parms->errcol, &errcol, &status);
    if (status) {
      fprintf(stderr, "ERROR: For background subtracted data, the %s column must be present,\n",
	      parms->errcol);
      fprintf(stderr, "       but was not found.\n");
      return status;
    }

    fits_read_col(fitsptr, TDOUBLE, errcol, 1, 1, *nrows, 0, rerror, 0, &status);
    if (status) {
      fprintf(stderr, "ERROR: could not read ERROR data\n");
      return status;
    }
  }


  /* Scan for bad data points */
  j = 0;
  if (parms->gaussian && rerror) {
    for (i=0; i<*nrows; i++) {
      if ((rate[i] == rate[i]) && (rerror[i] == rerror[i]) &&
	  (t[i] == t[i]) && (dt[i] == dt[i])) {
	rate[j]   = rate[i];
	rerror[j] = rerror[i];
	t[j]      = t[i];
	dt[j]     = dt[i];
	j ++;
      }
    }
  } else {
    for (i=0; i<*nrows; i++) {
      if ((rate[i] == rate[i]) && 
	  (t[i] == t[i]) && (dt[i] == dt[i])) {
	rate[j]   = rate[i];
	t[j]      = t[i];
	dt[j]     = dt[i];
	j ++;
      }
    }
  }

  *nrows = j;
  if (*nrows == 0) {
    fprintf(stderr, "ERROR: All the light curve rows were NULL!\n");
    return BAD_ROW_NUM;
  }
  

  return status;
}


/* ============================================================= */
/* Write a diagnostic column */
int wtcol(fitsfile *fitsptr, int *colnum, char *ttype, char *tform, 
	  int typecode, long int nrows, int repeat, void *data, 
	  int *status)
{
  fits_insert_col(fitsptr, *colnum, ttype, tform, status);
  fits_write_col(fitsptr, typecode, *colnum, 1, 1, repeat*nrows, 
		 data, status);
  (*colnum)++;
  if (*status) fprintf(stderr, "ERROR: could not write TIME column\n");
  return *status;
}

/* ============================================================= */
/* Write diagnostic output file */
int wtdiag(char *diagfile, long int nrows, long int ncells, 
	   double *t, double *dt, double *rate, double *rerror, 
	   int *cellsizes, int *cellpops, double *cellparams, 
	   int *last_start, double *best, double *cumsum)
{
  fitsfile *fitsptr;
  int status = 0;
  int colnum = 1;

  headas_clobberfile(diagfile);
  fits_create_file(&fitsptr, diagfile, &status);
  if (status) {
    fprintf(stderr, "ERROR: could not open diagnostic file\n");
    return status;
  }

  fits_create_tbl(fitsptr, BINARY_TBL, 0, 0,
		  0 /* colnames */, 0 /* colforms */, 0 /* colunits */,
		  "BATTBLOCKS_DIAG", &status);
  if (status) {
    fprintf(stderr, "ERROR: could create diagnostic table\n");
    return status;
  }

  headas_chat(5, "   (t=%s dt=%s nrows=%d ncells=%d)\n",
	      t  ? "found" : "none", dt ? "found" : "none",
	      nrows, ncells);

  fits_update_key(fitsptr, TLONG, "NROWS", &nrows, 
		  "Number of rows in input lt curve", &status);
  fits_update_key(fitsptr, TLONG, "NCELLS", &ncells, 
		  "Number of cells in Bayesian analysis", &status);

  if (t) {
    wtcol(fitsptr, &colnum, "TIME", "1D", TDOUBLE, nrows, 1, t, &status);
    if (status) return status;
  }
  
  if (dt) { 
    wtcol(fitsptr, &colnum, "TIMEDEL", "1D", TDOUBLE, nrows, 1, 
	  dt, &status);
    if (status) return status;
  }
  
  if (rate) {
    wtcol(fitsptr, &colnum, "RATE", "1D", TDOUBLE, nrows, 1, 
	  rate, &status);
    if (status) return status;
  }
  if (rerror) {
    wtcol(fitsptr, &colnum, "ERROR", "1D", TDOUBLE, nrows, 1, 
	  rerror, &status);
    if (status) return status;
  }
  
  if (cellsizes) {
    wtcol(fitsptr, &colnum, "CELL_SIZES", "1J", TINT, ncells, 1, 
	  cellsizes, &status);
    if (status) return status;
  }
  
  if (cellpops) {
    wtcol(fitsptr, &colnum, "CELL_POPS", "1J", TINT, ncells, 1, 
	  cellpops, &status);
    if (status) return status;
  }
  
  if (cellparams) {
    char tform[100];
    sprintf(tform, "%dD", NGAUSSCELLP);
    wtcol(fitsptr, &colnum, "CELL_PARAMS", tform, TDOUBLE, ncells, 
	  NGAUSSCELLP, cellparams, &status);
      if (status) return status;
  }
  
  if (last_start) {
    wtcol(fitsptr, &colnum, "LAST_START", "1J", TINT, ncells, 1, 
	  last_start, &status);
    if (status) return status;
  }
  
  if (best) {
      wtcol(fitsptr, &colnum, "BEST_LOGPROB", "1D", TDOUBLE, ncells, 1, 
	    best, &status);
      if (status) return status;
  }

  if (cumsum) {
      wtcol(fitsptr, &colnum, "CUMSUM", "1D", TDOUBLE, nrows, 1, 
	    cumsum, &status);
      if (status) return status;
  }
  
  fits_close_file(fitsptr, &status);
  return status;
}

/* ================================
 * apply_breakfile - apply breakfile operation
 * 
 * char *breakfile - name of breakfile
 * struct gti_struct *bbgti - input bayesian block GTI
 *             (out) modified blocks based on breakfile
 * int *status - CFITSIO status code pointer
 *
 * RETURNS: CFITSIO status code
 *
 */
int apply_breakfile(char *breakfile, 
		    struct gti_struct *bbgti, int *status)
{
  int i;
  struct gti_struct breakgti;
  struct gti_struct bbtmp;

  headas_chat(5, "  ... performing breakfile operation ...\n");
  HDgti_init(&bbtmp);
  HDgti_init(&breakgti);

  HDgti_copy(&bbtmp, bbgti, status);
    
  headas_chat(5, "  ... reading breakfile ...\n");
  HDgti_read(breakfile, &breakgti, 0, 0, 0, bbgti, 0, status);
  if (*status) {
    fprintf(stderr, "ERROR: could not read breakfile %s\n", 
	    breakfile);
    return *status;
  }
  for (i=0; i<bbtmp.ngti; i++) {
    bbtmp.start[i] += TFUDGE; bbtmp.stop[i]  -= TFUDGE;
  }
  for (i=0; i<breakgti.ngti; i++) {
    breakgti.start[i] += TFUDGE; breakgti.stop[i]  -= TFUDGE;
  }
    
  headas_chat(5, "  ... merging breakfile ...\n");
  HDgti_merge(GTI_AND, bbgti, /* <---- */ &bbtmp, &breakgti, status);
  if (*status) {
    fprintf(stderr, "ERROR: GTI intersection failed during breakfile\n");
    return *status;
  }

  /* Unfudge the times */
  for (i=0; i<bbgti->ngti; i++) {
    bbgti->start[i] -= TFUDGE;
    bbgti->stop[i]  += TFUDGE;
  }

  /* Deallocate storage */
  HDgti_free(&bbtmp); 
  HDgti_free(&breakgti);

  return *status;
  /* Final result is stored in bbgti */
}


int bkgsubtract(double *t, double *dt, double *cumsum, int n,
		double tstart, double tstop, double bkgstart, double bkgstop)
{
  double dti, bkgi, bkgcum, alpha;
  int i;

  if (tstart == 0 && tstop == 0) return 0;

  if (tstart == 0) {
    tstart = tstop-1;
    bkgstart = bkgstop;
  }

  if (tstop == 0) {
    tstop = tstart+1;
    bkgstop = bkgstart;
  }

  if (tstart == tstop) {
    double avg = 0.5*(bkgstart + bkgstop);
    tstop = tstart + 1;
    bkgstart = avg;
    bkgstop = avg;
  }    

  bkgcum = 0;
  alpha = (bkgstop-bkgstart)/(tstop-tstart);
  for (i=0; i<n; i++) {
    if (dt) {
      dti = dt[i];
    } else {
      dti = t[i+1]-t[i];
      if (i == (n-1)) dti = t[i]-t[i-1];
    }

    bkgi = (t[i]-tstart)*alpha + bkgstart;
    bkgcum += bkgi*dti;
    
    cumsum[i] -= bkgcum;
  }

  return 1;
}

/* ============================================================= */
int battblocks_work(struct parm_struct *parms)
{
  fitsfile *fitsptr = 0;
  int timecol = 0;
  int status = 0;
  long int nrows = 0, istart, istop, nevts, ngood;
  int ncells = 0;
  int *cellsizes = 0, *cellpops = 0;
  int nspill = 1;
  double *t = 0, *rate = 0, *dt = 0, *rerror = 0;
  double *cellparams = 0, data_range;
  double tstart, tstop, timedel;
  double burst_tstart, burst_tstop;
  double rms;
  int i;

  /* Returned by the Bayesian block optimizer */
  double *best = 0;
  int *last_start = 0;

  /* Some diagnostic arrays, and log prior probability */
  int *cparray = 0, ncp = 0;
  double ncp_prior = 6.0;

  /* For second pass blocking */
  int *cparray2 = 0;
  int ncp2 = 0;

  int nlag = 0;  /* Lookback samples */

  double *cumsum = 0; /* For computing cumulative sum counts */
  double tpeakwid;
  /* TXXGTI[] keeps track of various burst duration start/stop times 
     (see "Initialize GTIs" below) */
  struct gti_struct txxgti[7];
  double txxlist[3] = {90.0, 50.0, 0.0}; /* Percentage points */
  double txx_uncertainty[3];
  int ntxx = 2;
  double t90 = 0.0, t50 = 0.0, tuser = 0.0, tpeak = 0.0, totdur = 0.0;

  /* Estimated background rates at start & stop of interval */
  double bkgrate_start = 0.0, bkgrate_stop = 0.0;
  double tbkg_start = 0.0, tbkg_stop = 0.0;
  double fluence;

  /* FITS GTI output variables */
  char *ttype[] = {"START", "STOP"};
  char *extname = "STDGTI";
  struct gti_struct bbgti;
  char hduclas1[FLEN_CARD];
  char hduclas2[FLEN_CARD];
  char hduclas3[FLEN_CARD];

  /* ---- Initialize GTIs */
  HDgti_init(&bbgti);     /* Bayesian block GTI */
  HDgti_init(&txxgti[0]); /* T90 GTI */
  HDgti_init(&txxgti[1]); /* T50 GTI */
  HDgti_init(&txxgti[2]); /* Txx GTI (user) */
  HDgti_init(&txxgti[3]); /* Peak n second interval */
  HDgti_init(&txxgti[4]); /* Total burst duration */
  HDgti_init(&txxgti[5]); /* Background 1 */
  HDgti_init(&txxgti[6]); /* Background 2 */
    
  /* ---- Initialize Local variables */
  timedel = parms->timedel;
  nspill = parms->nspill;
  ncp_prior = parms->ncp_prior;
  tpeakwid = parms->tpeak;
  txxlist[2] = parms->txx;
  if (txxlist[2] > 0) ntxx = 3;

  /* ---- Print informative banner */
  banner(parms);

  /* Open input data */
  fits_open_data(&fitsptr, parms->infile, READONLY, &status);
  if (status) {
    fprintf(stderr, "ERROR: could not open %s\n", parms->infile);
    return status;
  }

  /* ------------------------- */
  /* Read basic keywords */
  fits_get_num_rows(fitsptr, &nrows, &status);
  fits_get_colnum(fitsptr, CASEINSEN, parms->timecol, &timecol, &status);
  tstart = HDget_frac_time(fitsptr, "TSTART", 0, 0, &status);
  tstop  = HDget_frac_time(fitsptr, "TSTOP", 0, 0, &status);
  if (status) {
    fprintf(stderr, 
	    "WARNING: could not read TSTART/TSTOP keywords\n"
	    "         assuming they are zero!!\n");
    tstart = tstop = 0;
    status = 0;
  }
  fits_read_key(fitsptr, TSTRING, "HDUCLAS1", hduclas1, 0, &status);
  if (status) {
    fprintf(stderr, 
	    "ERROR: could not read the HDUCLAS1 keyword\n"
	    "       this keyword is required to know the light curve type.\n");
    goto SAFE_RETURN;
  }
	    
  bbgti.mjdref = HDget_frac_time(fitsptr, "MJDREF", 0, 0, &status);
  if (status) {
    fprintf(stderr, 
	    "WARNING: could not read the MJDREF keyword\n"
	    "         assuming it is zero!!\n");
    bbgti.mjdref = 0;
    status = 0;
  }
  for (i=0; i<7; i++) { txxgti[i].mjdref = bbgti.mjdref; }

  if (status) {
    int safestatus;
  SAFE_RETURN:

    safestatus = 0;
    fprintf(stderr, "ERROR: could not read keywords from %s\n", parms->infile);
    fits_close_file(fitsptr, &safestatus);
    fitsptr = 0;
    return status;
  }
  headas_chat(5, "   (tstart=%f tstop=%f HDUCLAS1=%s)\n",
	      tstart, tstop, hduclas1);

  /* ------------------------- */
  /* Estimate lookback samples, based on lookback time, and mean
     sampling rate */
  nlag = 0;
  if (parms->tlookback > 0) {
    nlag = ((double) nrows * parms->tlookback / (tstop-tstart));
    headas_chat(5, "   (tlookback=%f nlookback=%d)\n", 
		parms->tlookback, nlag);
  }

  /* ------------------------- */
  /* Question: is it EVENTS or is it a light curve? */
  if ((strcasecmp(hduclas1, "EVENTS") == 0) || 
      (strcasecmp(hduclas1, "EVENT") == 0)) {
    parms->intype = EVENTS;
    headas_chat(5, " ...EVENTS data found...\n");
  } else if ( (strcasecmp(hduclas1, "LIGHTCURVE") == 0) || 
	      (strcasecmp(hduclas1, "LIGHT CURVE") == 0) ) {
    parms->intype = LC;
    fits_read_key(fitsptr, TSTRING, "HDUCLAS2", hduclas2, 0, &status);
    if (status) goto SAFE_RETURN;
    headas_chat(5, " ...LIGHTCURVE data found (HDUCLAS2=%s)...\n",hduclas2);

    /* Determine whether it's a background subtracted light curve; if
       yes, then set parms->gaussian if it hasn't been set already. */
    if ((strcasecmp(hduclas2, "NET") == 0) && (parms->gaussian == -1)) {
      parms->gaussian = 1;
    }
    if (parms->gaussian == -1) {
      parms->gaussian = 0;
    }
    headas_chat(5, " ...gaussian=%s...\n", parms->gaussian?"yes":"no");

  } else {
    fprintf(stderr, "ERROR: %s must be either EVENTS or LIGHTCURVE\n",
	    parms->infile);
    fits_close_file(fitsptr, &status);
    fitsptr = 0;
    return -1;
  }

  /* ------------------------- */
  headas_chat(5, " ...allocating memory for TIME...\n");
  t = (double *) malloc(sizeof(double)*nrows);
  if (t == 0) {
    fprintf(stderr, "ERROR: could not allocate memory for TIME\n");
    fits_close_file(fitsptr, &status);
    fitsptr = 0;
    return MEMORY_ALLOCATION;
  }

  fits_read_col(fitsptr, TDOUBLE, timecol, 1, 1, nrows, 0, t, 0, &status);

  if (parms->intype == EVENTS) {
    /* -----------------------   EVENT DATA   ------------------ */

    nlag /= nspill; 

    /* Close the input file */
    fits_close_file(fitsptr, &status);
    fitsptr = 0;

    /* Skip the data which precedes TSTART */
    headas_chat(5, " ...crude filtering of TIME...\n");
    for (istart=0; (istart < nrows) && (t[istart] < tstart); istart++);
    for (nevts = nrows - istart; (nevts > 0) && (t[istart+nevts-1] > tstop); nevts--);
    headas_chat(5, "   (nrows=%d nevts=%d istart=%d)\n", nrows, nevts, istart);

    /* Convert the events to "cells", possibly decimating by nspill.
       I.e., only each nspill'th event is selected. */
    evt2cells(t+istart, nevts, &cellsizes, &cellpops, &ncells, 
	      tstart, tstop, timedel, nspill);
    headas_chat(5, " ...calling evtbayes() w/ ncells=%d...\n", ncells);
    
    /* Estimate the Bayesian block change points, based on the event
       cost function. */
    cparray = evtbayes(cellsizes, cellpops, ncells, ncp_prior, &ncp,
		       &best, &last_start, nlag);

    /* Multiply by nspill, so that cparray[] again refers to the input
       event list */
    for (i=0; i<ncp; i++) {
      cparray[i] *= nspill;
      if (cparray[i] >= nevts) cparray[i] = nevts-1;
    }
    printcps(ncp, cparray);

    /* Second pass if we did the nibble the first time */
    if (parms->tlookback != 0) {

      headas_chat(5, " ...Reiterating on blocked data..\n");

      /* Rebin events into the new Bayesian block cells */
      free(cellsizes); free(cellpops); cellsizes = cellpops = 0;
      rebinevts(t+istart, nevts, &cellsizes, &cellpops, cparray, ncp, timedel);
      ncells = ncp-1;

      headas_chat(5, " ...calling evtbayes() w/ ncells=%d...\n", ncells);
      if (best) free(best);       best = 0;
      if (last_start) free(last_start); last_start = 0;

      /* Find a new change point array */
      cparray2 = evtbayes(cellsizes, cellpops, ncells, ncp_prior, &ncp2,
			  &best, &last_start, 0);

      /* Convert back to the original indexing system */
      for (i=0; i<ncp2-1; i++) {
	cparray2[i] = cparray[cparray2[i]];
      }
      cparray2[ncp2-1] = cparray[ncp-1];

      free(cparray); cparray = cparray2;
      ncp = ncp2;

      printcps(ncp, cparray);
    }


    /* ================= Coalesce first or last blocks if they are outliers */
    if (ncp > 2 && parms->coalescefrac > 0) {
      double t0start = t[istart+cparray[0]];
      double t1start = t[istart+cparray[1]];
      double t1stop  = t[istart+cparray[2]];
      double t0stop = t1start;

      /* Remove the first changepoint */
      if ((t0stop-t0start) < (t1stop-t1start)*parms->coalescefrac) {
	headas_chat(5," ...coalescing first two blocks (%f vs %f secs)...\n",
		    t0stop-t0start, t1stop-t1start);
	for (i=1; i<ncp-1; i++) cparray[i] = cparray[i+1];
	ncp--;
      }
    }

    if (ncp > 2 && parms->coalescefrac > 0) {
      double t0start = t[istart+cparray[ncp-3]];
      double t1start = t[istart+cparray[ncp-2]];
      double t1stop  = t[istart+cparray[ncp-1]];
      double t0stop  = t1start;

      /* Remove the last changepoint */
      if ((t1stop-t1start) < (t0stop-t0start)*parms->coalescefrac) {
	headas_chat(5," ...coalescing last two blocks (%f vs %f secs)...\n",
		    t0stop-t0start, t1stop-t1start);
	cparray[ncp-2] = cparray[ncp-1];
	ncp--;
      }
    }


    /* Prepare the GTI for the Bayesian blocks... */
    HDgti_grow(&bbgti, ncp-1, &status);
    if (status) {
      fprintf(stderr, "ERROR: could not allocate memory for GTI\n");
      return status;
    }

    /* ... and fill it */
    for (i=0; i<ncp-1; i++) {
      bbgti.start[i] = t[istart+cparray[i]];
      bbgti.stop[i]  = t[istart+cparray[i+1]];
    }
    bbgti.ngti = ncp-1;

    /* Estimate starting and stopping background rates */
    if (bbgti.stop[0] > bbgti.start[0]) {
      bkgrate_start = (cparray[1] - cparray[0]) / 
	(bbgti.stop[0] - bbgti.start[0]);
      tbkg_start = (bbgti.stop[0] + bbgti.start[0])/2;
    }
    if (bbgti.stop[ncp-2] > bbgti.start[ncp-2]) {
      bkgrate_stop  = (cparray[ncp-1] - cparray[ncp-2]) / 
	(bbgti.stop[ncp-2] - bbgti.start[ncp-2]);
      tbkg_stop = (bbgti.stop[ncp-2] + bbgti.start[ncp-2])/2;
    }
    
    /* Compute the cumulative total number of counts, matched to the
       time array */
    cumsum = evtcumsum(nrows);

    /* Perform background subtraction, in preparation for the fluence
       and duration-type activities that are coming up */
    if (parms->bkgsub && cumsum) {
      bkgsubtract(t, 0, cumsum, nrows, 
		  tbkg_start, tbkg_stop, bkgrate_start, bkgrate_stop);
    }

  } else {
    /* -----------------------   LIGHTCURVE DATA   ------------------ */

    rate   = (double *) malloc(sizeof(double)*nrows);
    rerror = (double *) malloc(sizeof(double)*nrows);
    dt     = (double *) malloc(sizeof(double)*nrows);
    if ((rate == 0) || (rerror == 0) || (dt == 0)) {
      fprintf(stderr, "ERROR: could not allocate memory for light curve data\n");
      status = MEMORY_ALLOCATION;
      goto CLEANUP;
    }

    /* Read light curve data */
    status = readlc(fitsptr, parms, t, rate, rerror, dt, &nrows, hduclas3);
    {
      int mystatus = 0;
      fits_close_file(fitsptr, &mystatus);
      fitsptr = 0;
    }
    if (status) goto CLEANUP;

    if (parms->gaussian == 0) {
      /* ============ POISSON binned data */
      if (strcasecmp(hduclas3, "RATE") == 0) {
	/* Convert RATE to counts */
	for (i=0; i<nrows; i++) rate[i] *= dt[i];
      } else if (strcasecmp(hduclas3, "COUNT") == 0) {
	/* Do nothing, already in counts */
      } else {
	fprintf(stderr, "ERROR: Invalid HDUCLAS3 found (='%s')\n", hduclas3);
	goto CLEANUP;
      }

      /* Bayesian blocks for pure Poisson counting data */
      headas_chat(5, " ...Calling lc2cells...\n");
      lc2cells(t, rate, dt, nrows, &cellsizes, &cellpops, &ncells, 
	       tstart, tstop, timedel);
      headas_chat(5, " ...calling lcbayes() w/ ncells=%d...\n", ncells);

      /* Retrieve the change points array (i.e. the edges of the
         Bayesian blocks) */
      cparray = lcbayes(cellsizes, cellpops, ncells, ncp_prior, &ncp,
			&best, &last_start, nlag);
      printcps(ncp, cparray);

      /* ================= Re-iterate if we nibbled */
      if (parms->tlookback != 0) {
	headas_chat(5, " ...Reiterating on blocked data..\n");

	/* Bayesian blocks for pure Poisson counting data */
	headas_chat(5, " ...Calling rebinlc...\n");

	/* Rebin according to the initial set of change points */
	free(cellsizes); free(cellpops); cellsizes = cellpops = 0;
	rebinlc(dt, rate, &cellsizes, &cellpops, cparray, ncp, timedel);
	ncells = ncp-1;

	headas_chat(5, " ...calling lcbayes() w/ ncells=%d...\n", ncells);
	if (best) free(best);             best = 0;
	if (last_start) free(last_start); last_start = 0;

	/* Estimate new Bayesian blocks */
	cparray2 = lcbayes(cellsizes, cellpops, ncells, ncp_prior, &ncp2,
			   &best, &last_start, 0);
      }

      /* Compute the cumulative total number of counts, matched to the
         time array */
      cumsum = lccumsum(rate, nrows);



      /* End of Poisson block */
    } else {
      /* ============ GAUSSIAN binned data */

      /* No need to rescale */

      headas_chat(5, " ...Calling lcgauss2cells...\n");
      /* Convert the bins to "cells", which are really the same thing */
      lcgauss2cells(t, rate, rerror, nrows, &cellparams, &ncells, 
		    &data_range, tstart, tstop, timedel);

      headas_chat(5, " ...calling lcgaussbayes() w/ ncells=%d...\n", ncells);
      /* Retrieve the change points array, i.e. the edges of the
         Bayesian blocks. */
      cparray = lcgaussbayes(cellparams, ncells, data_range,
			     ncp_prior, &ncp, 
			     &best, &last_start, nlag);
      printcps(ncp, cparray);

      /* ================= Re-iterate if we nibbled */
      if (parms->tlookback != 0) {
	headas_chat(5, " ...Reiterating on blocked data..\n");

	/* Bayesian blocks for Gaussian data */
	headas_chat(5, " ...Calling rebinlcgauss...\n");
	free(cellparams); cellparams = 0;
	rebinlcgauss(rate, rerror, &cellparams, &data_range,
		     cparray, ncp, timedel);
	ncells = ncp-1;

	headas_chat(5, " ...calling lcgaussbayes() w/ ncells=%d...\n", ncells);
	if (best) free(best);             best = 0;
	if (last_start) free(last_start); last_start = 0;
	cparray2 = lcgaussbayes(cellparams, ncells, data_range,
				ncp_prior, &ncp2, 
				&best, &last_start, 0);

      }

      /* Compute the cumulative total number of counts, matched to the
         time array */
      cumsum = lcgausscumsum(dt, rate, rerror, nrows, hduclas3);

    }  /* End of gaussian block */


    /* If we nibbled and then ran the second iteration, then replace
       the cparray[]. */
    if (parms->tlookback != 0) {
      for (i=0; i<ncp2-1; i++) {
	cparray2[i] = cparray[cparray2[i]];
      }
      cparray2[ncp2-1] = cparray[ncp-1];

      free(cparray); cparray = cparray2;
      ncp = ncp2;

      printcps(ncp, cparray);
    }

    /* ================= Coalesce first or last blocks if they are outliers */
    if (ncp > 2 && parms->coalescefrac > 0) {
      double t0start = t[cparray[0]] - 0.5*dt[cparray[0]];
      double t1start = t[cparray[1]] - 0.5*dt[cparray[1]];
      double t0stop = t1start;
      double t1stop;

      if (cparray[2] == nrows) {
	t1stop = t[nrows-1] + 0.5*dt[nrows-1];
      } else {
	t1stop = t[cparray[2]] - 0.5*dt[cparray[2]];
      }

      /* Remove the first changepoint */
      if ((t0stop-t0start) < (t1stop-t1start)*parms->coalescefrac) {
	headas_chat(5," ...coalescing first two blocks (%f vs %f secs)...\n",
		    t0stop-t0start, t1stop-t1start);
	for (i=1; i<ncp-1; i++) cparray[i] = cparray[i+1];
	ncp--;
      }
    }

    if (ncp > 2 && parms->coalescefrac > 0) {
      double t0start = t[cparray[ncp-3]] - 0.5*dt[cparray[ncp-3]];
      double t1start = t[cparray[ncp-2]] - 0.5*dt[cparray[ncp-2]];
      double t0stop = t1start;
      double t1stop;

      if (cparray[ncp-1] == nrows) {
	t1stop = t[nrows-1] + 0.5*dt[nrows-1];
      } else {
	t1stop  = t[cparray[ncp-1]] - 0.5*dt[cparray[ncp-1]];
      }

      /* Remove the last changepoint */
      if ((t1stop-t1start) < (t0stop-t0start)*parms->coalescefrac) {
	headas_chat(5," ...coalescing last two blocks (%f vs %f secs)...\n",
		    t0stop-t0start, t1stop-t1start);
	cparray[ncp-2] = cparray[ncp-1];
	ncp--;
      }
    }


    /* ================================= */
    /* Allocate storage for the GTI */
    HDgti_grow(&bbgti, ncp-1, &status);
    if (status) {
      fprintf(stderr, "ERROR: could not allocate memory for GTI\n");
      return status;
    }

    /* Fill the GTI.  The last point is tricky, since it refers to the
       N+1'th cell. */
    for (i=0; i<ncp-1; i++) {
      bbgti.start[i] = t[cparray[i]] - 0.5*dt[cparray[i]];
      if (cparray[i+1] == nrows) {
	bbgti.stop[i]  = t[nrows-1] + 0.5*dt[nrows-1];
      } else {
	bbgti.stop[i]  = t[cparray[i+1]] - 0.5*dt[cparray[i+1]];
      }
    }
    bbgti.ngti = ncp-1;

    /* ==============  Estimate the background rates */
    {
      double cts = 0.0;

      for (i=cparray[0]; i<cparray[1]; i++) cts += rate[i];

      if (parms->gaussian == 0) {
	if ((bbgti.stop[0]-bbgti.start[0]) > 0) 
	  bkgrate_start = cts / (bbgti.stop[0]-bbgti.start[0]);
      } else {
	if (cparray[1] > cparray[0]) 
	  bkgrate_start = cts / (cparray[1] - cparray[0]);
      }
      tbkg_start = (bbgti.stop[0] + bbgti.start[0])/2;
      
      cts = 0.0;
      
      for (i=cparray[ncp-2]; i<cparray[ncp-1]; i++) cts += rate[i];

      if (parms->gaussian == 0) {
	if ((bbgti.stop[ncp-2]-bbgti.start[ncp-2]) > 0) 
	  bkgrate_stop = cts / (bbgti.stop[ncp-2]-bbgti.start[ncp-2]);
      } else {
	if (cparray[ncp-1] > cparray[ncp-2])
	  bkgrate_stop = cts / (cparray[ncp-1] - cparray[ncp-2]);
      }
      tbkg_stop = (bbgti.stop[ncp-2] + bbgti.start[ncp-2])/2;

      if (parms->bkgsub && cumsum) {
	bkgsubtract(t, dt, cumsum, nrows, 
		    tbkg_start, tbkg_stop, bkgrate_start, bkgrate_stop);
      }
	
    }


  } /* End of LIGHTCURVE processing */

  /* Estimate the burst duration, from the end of the first BB to the
     beginning of the last BB */
  burst_tstart = burst_tstop = FLOAT_INDEF;
  if (bbgti.ngti >= 3) {
    burst_tstart = bbgti.stop[0];
    burst_tstop  = bbgti.start[bbgti.ngti-1];
  }
  /* ... or based on command line override */
  if (parms->burst_tstart != FLOAT_INDEF) { burst_tstart = parms->burst_tstart; }
  if (parms->burst_tstop  != FLOAT_INDEF) { burst_tstop  = parms->burst_tstop; }
  headas_chat(5,"  ...burst start/stop: %f/%f...\n",
	      burst_tstart, burst_tstop);
  

  /* Fill in the out-going T100 good time interval (txxgti[4]) */
  HDgti_grow(&txxgti[4], 1, &status);
  if (status) { 
    fprintf(stderr, "ERROR: gti grow failed for total duration\n");
    goto CLEANUP;
  }
  txxgti[4].ngti = 1;
  txxgti[4].start[0] = burst_tstart;
  txxgti[4].stop[0]  = burst_tstop;

  /* ------------------------- */
  /* Estimate the burst interval */
  if (burstspan(t, dt, nrows, burst_tstart, burst_tstop, &istart, &istop) <= 0) {
    fprintf(stderr, "ERROR: Not enough Bayesian blocks were found to\n");
    fprintf(stderr, "       select a burst interval\n");
    ngood = 0;
    rms = 1.0;

  } else {
    double diff = cumsum[istop] - cumsum[istart];
    if (diff <= 0) diff = 1.0;   /* XXX tragedy alert! */
    
    ngood = istop - istart + 1;

    /* Estimate the standard error of the mean for efficiencies */
    if (parms->gaussian == 1) {
      /* Special case: gaussian errors tacked at the end */
      double *cumerr2 = cumsum + nrows;  
      
      rms = cumerr2[istop] - cumerr2[istart];
      rms = sqrt(rms) / diff;
    } else {
      rms = 1.0 / sqrt(diff);
    }

    headas_chat(5,"  (sigma rms = %f)\n", rms);
  }

  /* ------------------------- */
  /* Estimate the fluence */
  if (cumsum && ngood > 0) {
    fluence = cumsum[istop] - cumsum[istart];
    headas_chat(5, "  ...fluence start/stop: %d/%d %f/%f %f/%f\n",
		istart, istop, t[istart], t[istop], cumsum[istart], cumsum[istop]);
  }

  /* ------------------------- */
  /* Estimate the burst duration intervals */
  if (parms->durfile[0] && cumsum && ngood > 0) {
    burstdur(t+istart, dt?(dt+istart):0, cumsum+istart, ngood, txxlist, ntxx, 
	     tpeakwid, &txxgti[0], &txxgti[3], txx_uncertainty, rms, 
	     parms->durerrmeth);

    /* Compute durations */
    if (txxgti[0].ngti > 0) {
      t90 = txxgti[0].stop[txxgti[0].ngti-1] - txxgti[0].start[0];
    }
    if (txxgti[1].ngti > 0) {
      t50 = txxgti[1].stop[txxgti[1].ngti-1] - txxgti[1].start[0];
    }
    if (txxgti[2].ngti > 0) {
      tuser = txxgti[2].stop[txxgti[2].ngti-1] - txxgti[2].start[0];
    }
    if (txxgti[3].ngti > 0) {
      tpeak = 0.5*(txxgti[3].stop[txxgti[3].ngti-1] + txxgti[3].start[0]);
    }

  }

  if (txxgti[4].ngti > 0) {
    totdur = txxgti[4].stop[txxgti[4].ngti-1] - txxgti[4].start[0];
  }


  /* ======================== */
  /* If user has requested to break up the BB intervals at specific times,
     then honor that request now. */
  if (parms->breakfile[0]) {
    apply_breakfile(parms->breakfile, &bbgti, &status);
    if (status) {
      fprintf(stderr, "ERROR: breakfile operation failed.\n");
      return status;
    }
    ncp = bbgti.ngti + 1;
    /* Final result is stored in bbgti */
  }

  /* ------------------------- */
  /* Write Bayesian Block GTIs */
  headas_clobberfile(parms->outfile);
  if (fits_create_file(&fitsptr, parms->outfile, &status)) {
    fprintf(stderr, "ERROR: could not create %s\n", parms->outfile);
    return status;
  }


  /* Write the BB intervals to the GTI extension */
  HDgti_write(fitsptr, &bbgti, extname, ttype[0], ttype[1], &status);
  if (status) {
    int safestatus = 0;
    fprintf(stderr, "ERROR: could not write Bayesian GTI to %s\n", parms->outfile);
    fits_close_file(fitsptr, &safestatus);
    return status;
  }

  /* ------------------------- */
  /* Write keywords which describe the duration measures (if they were
     computed) */
  if (t90 > 0) {
    fits_update_key(fitsptr, TDOUBLE, "T90DUR", &t90,
		    "[s] Burst T90 duration measure", &status);
    fits_update_key(fitsptr, TDOUBLE, "T90ERR", txx_uncertainty+0,
		    "[s] Crude T90DUR error estimate", &status);
    headas_chat(1, "  Estimated T90 duration: %g s +/- %g s\n", 
		t90, txx_uncertainty[0]);
  }
  if (t50 > 0) {
    fits_update_key(fitsptr, TDOUBLE, "T50DUR", &t50,
		    "[s] Burst T50 duration measure", &status);
    fits_update_key(fitsptr, TDOUBLE, "T50ERR", txx_uncertainty+1,
		    "[s] Crude T50DUR error estimate", &status);
    headas_chat(1, "  Estimated T50 duration: %g s +/- %g s\n", 
		t50, txx_uncertainty[1]);
  }
  if (tuser > 0) {
    fits_update_key(fitsptr, TDOUBLE, "TUSERDUR", &tuser,
		    "[s] Burst TXX duration measure", &status);
    fits_update_key(fitsptr, TDOUBLE, "TUSERPCT", &parms->txx,
		    "Percentage of burst enclosed by TUSERDUR", &status);
    fits_update_key(fitsptr, TDOUBLE, "TUSERERR", txx_uncertainty+2,
		    "[s] Crude TUSERDUR error estimate", &status);
    headas_chat(1, "  Estimated T%2.2d duration: %g s +/- %g s\n", 
		(int) rint(parms->txx), tuser, txx_uncertainty[2]);
  }
  if (tpeak > 0) {
    fits_update_key(fitsptr, TDOUBLE, "TPEAK", &tpeak,
		    "[s] MET center of peak interval (GTI_PEAK)", &status);
    headas_chat(1, "  Estimated Peak Interval: MET %20.8f +/- %f s\n",
		tpeak, tpeakwid/2);
  }
  if (tbkg_start != 0) {
    fits_update_key(fitsptr, TDOUBLE, "TBKG1", &tbkg_start,
		    "[s] MET center of first background interval", &status);
    fits_update_key(fitsptr, TDOUBLE, "BKGRATE1", &bkgrate_start,
		    "[count/s] Estimated background rate (1st)", &status);
    headas_chat(1, "  Estimated background rate 1: %12.3f near MET %f s\n",
		bkgrate_start, tbkg_start);
    if (totdur > 0) {
      int safestatus = 0;
      HDgti_grow(&txxgti[5], 1, &safestatus);
      txxgti[5].ngti     = 1;
      txxgti[5].start[0] = parms->global_tstart;
      txxgti[5].stop[0]  = txxgti[4].start[0];
    }
  }
  if (tbkg_stop != 0) {
    fits_update_key(fitsptr, TDOUBLE, "TBKG2", &tbkg_stop,
		    "[s] MET center of second background interval", &status);
    fits_update_key(fitsptr, TDOUBLE, "BKGRATE2", &bkgrate_stop,
		    "[count/s] Estimated background rate (2nd)", &status);
    headas_chat(1, "  Estimated background rate 2: %12.3f near MET %f s\n",
		bkgrate_stop, tbkg_stop);
    if (totdur > 0) {
      int safestatus = 0;
      HDgti_grow(&txxgti[6], 1, &safestatus);
      txxgti[6].ngti     = 1;
      txxgti[6].start[0] = txxgti[4].stop[0];
      txxgti[6].stop[0]  = parms->global_tstop;
    }
  }
  if (totdur > 0) {
    fits_update_key(fitsptr, TDOUBLE, "TOTDUR", &totdur,
		    "[s] Duration of GRB (for data selection)", &status);
    fits_update_key(fitsptr, TDOUBLE, "TOTSTART", &(txxgti[4].start[0]),
		    "[s] Start of GRB (for data selection)", &status);
    fits_update_key(fitsptr, TDOUBLE, "TOTSTOP", &(txxgti[4].stop[0]),
		    "[s] Stop of GRB (for data selection)", &status);
    headas_chat(1, "  Estimated total duration: %g s (for data selection)\n",
		totdur);
    headas_chat(1, "       (from MET %f to MET %f)\n",
		txxgti[4].start[0], txxgti[4].stop[0]);
  }
    

  /* ------------------------- */
  /* Write background subtraction information */
  {
    int backapp = (parms->bkgsub || (parms->gaussian == 1));
    fits_update_key(fitsptr, TLOGICAL, "BACKAPP", &backapp,
		    "Background subtraction applied?", &status);
  }

  /* ------------------------- */
  /* Write fluence information */
  if (cumsum) {
    fits_update_key(fitsptr, TDOUBLE, "CNTFLU", &fluence,
		    "[count] Total counts fluence", &status);
    headas_chat(1, "  Estimated total fluence: %f count\n", fluence);
  }


  if (status) {
    int safestatus = 0;
    fits_close_file(fitsptr, &safestatus);
    fitsptr = 0;
    fprintf(stderr, "ERROR: could not write GTI to %s\n", parms->outfile);
    return status;
  }
  status = HDpar_stamp(fitsptr, 0, &status);
  fits_close_file(fitsptr, &status);
  fitsptr = 0;

  /* ------------------------- */
  /* Write burst duration measures to GTI extensions */
  if (parms->durfile[0]) {
    int created = 0;
    char *gtinames[] = {
      "GTI_T90", "GTI_T50", "GTI_TXX", "GTI_PEAK", "GTI_TOT",
      "GTI_BKG1", "GTI_BKG2"
    };

    for (i=0; i<7; i++) if (txxgti[i].ngti > 0) {

	/* Create the file if needed */
	if (!created) {
	  headas_chat(5, "  ...creating duration file %s...\n", 
		      parms->durfile);
	  headas_clobberfile(parms->durfile);
	  if (fits_create_file(&fitsptr, parms->durfile, &status)) {
	    fprintf(stderr, "ERROR: could not create %s\n", parms->durfile);
	    return status;
	  }
	  created = 1;
	}

	headas_chat(5, "  ...writing extension %s (%d entries)...\n",
		    gtinames[i], txxgti[i].ngti);
	HDgti_write(fitsptr, &txxgti[i], gtinames[i], 
		    ttype[0], ttype[1], &status);

	if (strncmp(gtinames[i],"GTI_BKG",7) == 0) {
	  fits_write_comment(fitsptr, 
	     "NOTE: values of -/+ 1.0E307 indicate start/end of observation", 
			     &status);
	}

	if (status) {
	  int safestatus = 0;
	  fits_close_file(fitsptr, &safestatus);
	  fitsptr = 0;
	  fprintf(stderr, "ERROR: could not write %s to %s\n", 
		  gtinames[i], parms->durfile);
	  return status;
	}
    }

    if (created) {
      fits_close_file(fitsptr, &status);
      fitsptr = 0;
    }
  }

  /* ------------------------- */
  summary(parms, ncp-1);

  /* Write diagnostic summary file */
  if (parms->diagfile[0]) {
    status = wtdiag(parms->diagfile, nrows, ncells, 
		    t, dt, rate, rerror, cellsizes, cellpops, cellparams, 
		    last_start, best, cumsum);
  }

 CLEANUP:
  /* Close file if it's open */
  if (fitsptr) {
    int safestatus = 0;
    fits_close_file(fitsptr, &safestatus);
    fitsptr = 0;
  }
  /* Free GTI storage */
  HDgti_free(&bbgti);

  /* Free malloc'd storage */
  if (cparray) free(cparray); cparray = 0;
  if (cellsizes) free(cellsizes); cellsizes = 0;
  if (cellpops) free(cellpops); cellpops = 0;
  if (t) free(t); t = 0;
  if (dt) free(dt); dt = 0;
  if (rate) free(rate); rate = 0;
  if (rerror) free(rerror); rerror = 0;
  if (cellparams) free(cellparams); cellparams = 0;
  if (best) free(best); best = 0;
  if (last_start) free(last_start); last_start = 0;
  if (cumsum) free(cumsum); cumsum = 0;
  
  return status;
}

/* ============================================================= */
int battblocks(void)
{
  int status = 0;

  struct parm_struct parms;

  /* Register taskname and version. */

  set_toolname(taskname);
  set_toolversion(taskver);

  if ((status = battblocks_getpar(&parms)) != 0) {
    fprintf(stderr, "Could not read parameter file\n");
    return status;
  }

  parms.taskname = &taskname[0];
  parms.taskver  = &taskver[0];

  return battblocks_work(&parms);

}
