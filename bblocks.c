#include <fitsio.h>
#include <math.h>
#include <string.h>

/* Bayesian blocks for pure Poisson counting data */
lc2cells(t, rate, dt, nrows, &cellsizes, &cellpops, &ncells, tstart, tstop, timedel);

/* Retrieve the change points array (i.e. the edges of the Bayesian blocks) */
cparray = lcbayes(cellsizes, cellpops, ncells, ncp_prior, &ncp, &best, &last_start, nlag);

/* Compute the cumulative total number of counts, matched to the time array */
cumsum = lccumsum(rate, nrows);


/* ================= Coalesce first or last blocks if they are outliers */
if (ncp > 2 && parms->coalescefrac > 0) {
  double t0start = t[cparray[0]] - 0.5*dt[cparray[0]];
  double t1start = t[cparray[1]] - 0.5*dt[cparray[1]];
  double t0stop = t1start;
  double t1stop;

  if (cparray[2] == nrows)
  	{ 
  		t1stop = t[nrows-1] + 0.5*dt[nrows-1];
	} else { 
		t1stop = t[cparray[2]] - 0.5*dt[cparray[2]];
	}

  /* Remove the first changepoint */
  if ((t0stop-t0start) < (t1stop-t1start)*parms->coalescefrac)
  {
	for (i=1; i<ncp-1; i++) cparray[i] = cparray[i+1];
	ncp--;
  }
}

if (ncp > 2 && parms->coalescefrac > 0) {
  double t0start = t[cparray[ncp-3]] - 0.5*dt[cparray[ncp-3]];
  double t1start = t[cparray[ncp-2]] - 0.5*dt[cparray[ncp-2]];
  double t0stop = t1start;
  double t1stop;

  if (cparray[ncp-1] == nrows) 
  {
	t1stop = t[nrows-1] + 0.5*dt[nrows-1];
  } else {
	t1stop  = t[cparray[ncp-1]] - 0.5*dt[cparray[ncp-1]];
  }

  /* Remove the last changepoint */
  if ((t1stop-t1start) < (t0stop-t0start)*parms->coalescefrac)
  {
	cparray[ncp-2] = cparray[ncp-1];
	ncp--;
  }
}

/* Fill the GTI.  The last point is tricky, since it refers to the N+1'th cell. */
for (i=0; i<ncp-1; i++) 
{
	bbgti.start[i] = t[cparray[i]] - 0.5*dt[cparray[i]];
	if (cparray[i+1] == nrows) 
	{
		bbgti.stop[i]  = t[nrows-1] + 0.5*dt[nrows-1];
	} else {
		bbgti.stop[i]  = t[cparray[i+1]] - 0.5*dt[cparray[i+1]];
	}
}
bbgti.ngti = ncp-1;

/* Estimate the burst duration, from the end of the first BB to the beginning of the last BB */
burst_tstart = bbgti.stop[0];
burst_tstop  = bbgti.start[bbgti.ngti-1];

/* Find istart and istop */
burstspan(t, dt, nrows, burst_tstart, burst_tstop, &istart, &istop)

/* Estimate the fluence */
fluence = cumsum[istop] - cumsum[istart];




int burstspan(double *t, double *dt, int nrows,
	      double tstart, double tstop,
	      long int *istart0, long int *istop0)
{
  int istart, istop;
  int i;

  /* Find start and stop times:
     1. start time is end of first BB, which is assumed to be the
     pre-burst background.
     2. stop is the start of the last BB, which is assumed to be the
     post-burst background.
  */
  /* 
  tstart = bbgti->stop[0];
  tstop  = bbgti->start[bbgti->ngti-1];
  */
  
  for (i=0; i<nrows; i++) {
    if (t[i] >= tstart) break;
  }
  istart = i-1;
  if (i == 0) istart = 0;
  for (i=istart+1; i<nrows; i++) {
    if (t[i] >= tstop) break;
  }
  istop = i;
  if (i == nrows) istop = nrows-1;

  headas_chat(5, "  ...time range start=%f(%d) stop=%f(%d)...\n",
	      t[istart], istart, t[istop], istop);

  *istart0 = istart;
  *istop0 = istop;
  
  return (istop-istart);

/* ============================================================= */
/* Form the cumulative sum of the light curve, between two points.
 *  The data is assumed to be expressed in counts already. 
 *  Times are assumed to be center-bin.
 */
double *lccumsum(double *counts, int ntimes)
{
  int i;
  double *cumcounts = 0;

  cumcounts = (double *)malloc(sizeof(double)*ntimes);
  if (cumcounts == 0) return 0;

  cumcounts[0] = counts[0];
  for (i=1; i<ntimes; i++) {
    cumcounts[i] = cumcounts[i-1] + counts[i];
  }

  return cumcounts;
}

/* ============================================================= */
/* Determine the Bayesian block change points, based on the 
   Poisson binned cost function */
int *lcbayes(int *cellsizes, int *cellpops, int ncells, 
	     double ncp_prior, int *ncparray,
	     double **bestlogprob, int **lastcellstart, int nlag)
{
  int *cumsizes = 0, *cumpops = 0, *last_start;
  int *cparray = 0;
  double *merged = 0, *best = 0;
  double temp;
  int i, j, imaxer, ncp, index, icp;
  int istart = 0, ioldstart = 0;

  if (bestlogprob) *bestlogprob = 0;
  if (lastcellstart) *lastcellstart = 0;

  cumsizes = (int *) malloc(sizeof(int)*ncells);
  cumpops  = (int *) malloc(sizeof(int)*ncells);
  last_start = (int *) malloc(sizeof(int)*ncells);
  merged   = (double *) malloc(sizeof(double)*ncells);
  best     = (double *) malloc(sizeof(double)*ncells);
  if ((cumsizes == 0) || (cumpops == 0) || (merged == 0) || 
      (best == 0) || (last_start == 0)) {
    if (cumsizes) free(cumsizes);
    if (cumpops) free(cumpops);
    if (merged) free(merged);
    if (best) free(best);
    if (last_start) free(last_start);
    return 0;
  }
  
  for (i=0; i<ncells; i++) {
    cumsizes[i] = 0;
    cumpops[i] = 0;
  }

  istart = 0; ioldstart = 0;
  for (i=0; i<ncells; i++) {
    /* Approximation to the "nibble" algorithm */
    if (nlag > 0) {
      istart = i - nlag;
      if (istart < 0) istart = 0;
    }

    /* If we are nibbling, then we must shift the best[] array, so
       that the normalized probability is unity before the starting
       element. */
    if ( (istart > 0) && (istart != ioldstart) ) {
      for (j=istart; j<i; j++) {
	best[j] -= best[istart-1];
      }
    }

    /* Accumulate the parameters */
    for(j=istart; j<i; j++) {
      cumsizes[j] += cellsizes[i];
      cumpops[j]  += cellpops[i];
    }
    cumsizes[i] = cellsizes[i];
    cumpops[i] = cellpops[i];

    /* Compute the cost function for the cumulants */
    logprob_lc(cumpops+istart, cumsizes+istart, i+1-istart, 
	       merged+istart, ncp_prior);

    /* Where is the maximum probability in the joint best|merged
       arrays? */
    imaxer = istart;
    best[i] = merged[istart];
    if (i > 0) {
      for(j=istart+1; j<=i; j++) {
	temp = best[j-1]+merged[j];
	if (temp > best[i]) {
	  best[i] = temp;
	  imaxer = j;
	}
      }
    }

    /* Record the new best position */
    last_start[i] = imaxer;

    /* Keep track of the previous nibble starting point */
    ioldstart = istart;
  }

#if 0
  /* Debugging output to a file */
  { 
    FILE *out;
    out = fopen("test.dat", "w");
    for (i=0; i<ncells; i++) {
      fprintf(out, "%d %d %f %f\n", i, last_start[i], best[i], merged[i]);
    }
    fclose(out);
  }
#endif

  /* Count number of change points */
  ncp = 2;
  index = last_start[ncells-1];
  while (index > 1) {
    ncp ++;
    index = last_start[index-1];
  }

  /* Create output array of change points */
  cparray = (int *) malloc(sizeof(int)*ncp);
  if (cparray == 0) {
    ncp = 0;
    goto CLEANUP;
  }
  
  icp = ncp-1;
  cparray[icp--] = ncells;
  index = last_start[ncells-1];
  while (index > 1) {
    cparray[icp--] = index;
    index = last_start[index-1];
  }
  cparray[0] = 0;

 CLEANUP:
  if (cumsizes) free(cumsizes);
  if (cumpops) free(cumpops);
  if (merged) free(merged);
  if (last_start) {
    if (lastcellstart) *lastcellstart = last_start;
    else free(last_start);
  }
  if (bestlogprob) {
    if (bestlogprob) *bestlogprob = best;
    else free(best);
  }

  *ncparray = ncp;
  return cparray;
}

/* ============================================================= */
/*
 * logprob_lc - Compute log posterior probability for binned data
 *
 * int *cellpops - populations of cells (i.e. number of events per cell)
 * int *cellsizes - widths of cells, in units of timedel
 * int ncells - number of cells
 * double *logprob - upon return, the log probability
 * double ncp_prior - log(prob) prior
 *
 * RETURNS: CFITSIO status value 
 *
;-----------------------------------------------------------------
; See:  J.D. Scargle, 1998, ApJ, 504, 405
;
; Log posterior (Bayes factor) for constant-rate Poisson data:
;   * flat prior on Poisson rate parameter (unnormalized)
;   * geometric prior on number of changepoints
;
; Input: cell_sizes -- size (length in 1D) of each cell (array)
;        cell_pops  -- number of events in each cell (array)
;        ncp_prior  -- log parameter for number of changepoints
;
; The first two inputs are computed from make_cells.  The third
; input, ncp_prior, acts like a smoothing parameter, weighting
; against separate blocks.  The optimal value for ncp_prior was
; determined empirically, adjusting its value until virtually
; no "spikes" remained, while retaining significant block
; structure across the time series.
;
; datatype: 3 --> binned data, binned posterior
;-----------------------------------------------------------------
*/

int logprob_lc(int *cellpops, int *cellsizes, int ncells, 
		double *logprob, double ncp_prior)
{
  int i;

  for (i=0; i<ncells; i++) {
    logprob[i] = lgamma(cellpops[i]+1) - (cellpops[i]+1)*log(cellsizes[i]);
    logprob[i] -= ncp_prior;
  }

  return 0;
}

/* ============================================================= */
/* Convert a light curve to cells, by direct transcription */
int lc2cells(double *t, double *counts, double *dt, int ntimes, 
	     int **cellsizes, int **cellpops, int *ncells, 
	     double tstart, double tstop, double timedel)
{
  int *cpops = 0, *csize = 0;
  int nc, i;

  if ((t == 0) || (counts == 0) || (dt == 0) || (ntimes <= 0) || 
      (cellsizes == 0) || (cellpops == 0) || (tstop - tstart <= 0)) {
    return 0;
  }

  nc = ntimes;

  cpops = (int *) malloc(sizeof(int)*nc);
  csize = (int *) malloc(sizeof(int)*nc);
  if ((cpops == 0) || (csize == 0)) {
    if (cpops) free(cpops);
    if (csize) free(csize);
    return 0;
  }

  for (i=0; i<nc; i++) {
    csize[i] = dt[i] / timedel;
    cpops[i] = counts[i];
  }

  *cellsizes = csize;
  *cellpops  = cpops;
  *ncells    = nc;
  
  return 0;
}