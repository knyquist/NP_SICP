#include "SICalgorithm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

step_fit NoStepSIC(double* pos,
		   int n,
		   double nu,
		   double So) {
  step_fit no_step;
  no_step = InitializeFitToZeros(no_step,1);
  for(int i=0;i<n;i++) no_step.means[0] += pos[i];
  no_step.means[0] = no_step.means[0]/n;
  for(int i=0;i<n;i++) no_step.chisq += pow(pos[i]-no_step.means[0],2.);
  no_step.chisq = (/*1840.5*/So+no_step.chisq)/(n+nu-2/*204.5*/);//nu=206.5; Sn=1840.5 for var=9,varvar=0.8
  no_step.SIC = 2*log(n)+(n+nu+2/*208.5*/)*log(no_step.chisq)+/*1840.5*/So/no_step.chisq;
  return no_step;
}

step_fit AddStepSIC(double* pos, 
		    int n, 
		    int* step_indices,
		    int n_dwell,
		    step_fit prev_fit,
		    double nu,
		    double So) {
  step_fit curr_fit; // current interation of the fit
  curr_fit = InitializeFitToZeros(curr_fit,n_dwell+1);
  step_fit win_fit; //running best fit during this AddStepSIC iteration
  win_fit = InitializeFitToZeros(win_fit,n_dwell+1); win_fit.SIC=-1.0;
  int win_step_location=0; //keep track of the best trial step location
  for(int i=1;i<n-1;i++) {
    if(step_indices[i]==0) {
      step_indices[i] = 1; // trial step location
      /* get means of the dwells */
      int j=i-1; double tmpleft=0.;
      while((step_indices[j]==0)&&(j>=0)){
	tmpleft+=pos[j];
	j-=1;
      }
      tmpleft=tmpleft/(i-j-1);
      j=i+1; double tmpright=0.;
      while((step_indices[j]==0)&&(j<=n-1)){
	tmpright+=pos[j];
	j+=1;
      }
      tmpright=tmpright/(j-i-1); j=0;
      /* get the new chisq (also remove the previous part) */
      double prev_mean=GetPreviousDwellMean(prev_fit,step_indices,i);
      curr_fit.chisq=prev_fit.chisq*(n+/*34.3*/nu-2)-/*580.*/So;
      j=i-1;
      while((step_indices[j]==0)&&(j>=0)) {
	curr_fit.chisq -= pow(pos[j]-prev_mean,2.);
	curr_fit.chisq += pow(pos[j]-tmpleft,2.);
	j-=1;
      }
      j=i;
      do {
     	curr_fit.chisq -= pow(pos[j]-prev_mean,2.);
	curr_fit.chisq += pow(pos[j]-tmpright,2.);
	j+=1;
      }
      while((step_indices[j]==0)&&(j<n));
      curr_fit.chisq=(/*580.*/So+curr_fit.chisq)/(n+nu-2/*34.3*/);
      curr_fit.SIC=((n_dwell)+2)*log(n)+(n+nu+2/*38.3*/)*log(curr_fit.chisq)+/*580.*/So/curr_fit.chisq;//n_dwell is now the number of steps for this added step.  it is updated after this completes but before the next iteration.
      if((curr_fit.SIC<win_fit.SIC)||(win_fit.SIC==-1.0)) {
	win_fit = SetFitToZeros(win_fit,n_dwell+1);
	win_step_location=i;
	for(int k=0;k<n_dwell+1;k++) curr_fit.means[k]=0.;
	int dum=0;
	for(int k=0;k<n_dwell;k++) {
	  if(prev_fit.means[k]!=prev_mean) {
	    curr_fit.means[dum] = prev_fit.means[k];dum+=1;
	  }
	  else if(prev_fit.means[k]==prev_mean) {
	    curr_fit.means[dum]=tmpleft;dum+=1;
	    curr_fit.means[dum]=tmpright;dum+=1; 
	  }
	}
	dum=0;
	for(int i=0;i<n_dwell+1;i++) win_fit.means[i]=curr_fit.means[i];
	win_fit.SIC = curr_fit.SIC;
	win_fit.chisq = curr_fit.chisq;
      }
      curr_fit = SetFitToZeros(curr_fit,n_dwell+1);
      step_indices[i]=0;
    }
    else if(step_indices[i]==1) {
      //do nothing
    }
    else {
      puts("Error: step_indices array has value that is not 0 or 1");
      exit(1);
    }
  }
  int dum=0; //update step_locations to include the best step
  for(int k=0;k<n;k++) {
    if(step_indices[k]==1) {
      win_fit.step_locations[dum]=k;
      dum+=1;
    }
    if(win_step_location==k) {
      win_fit.step_locations[dum]=win_step_location;
      dum+=1;
    }
  }
  return win_fit;
}

step_fit InitializeFitToZeros(step_fit fit, 
			      int length_of_mean_list) {
  fit.means = (double*) malloc(length_of_mean_list*sizeof(double));
  fit.step_locations = (int*) malloc((length_of_mean_list-1)*sizeof(int));
  for(int i=0;i<length_of_mean_list;i++) {
    fit.means[i] = 0.;
    if(i>0) fit.step_locations[i-1] = 0;
  }
  fit.SIC = 0.;
  fit.chisq = 0;
  return fit;
}

step_fit SetFitToZeros(step_fit fit, 
		       int length_of_mean_list) {
  for(int i=0;i<length_of_mean_list;i++) {
    fit.means[i] = 0.;
    if(i>0) fit.step_locations[i-1] = 0;
  }
  fit.SIC = 0.;
  fit.chisq = 0;
  return fit;
}

step_fit MakeFitOneDwellBigger(step_fit fit, 
			       int newlength) {
  fit.means = (double*) realloc(fit.means,(newlength*sizeof(double)));
  fit.step_locations = (int*) realloc(fit.step_locations,((newlength-1)*sizeof(int)));
  if((fit.means!=NULL)||(fit.step_locations!=NULL)) {
    return fit;
  }
  else {
    free(fit.means); free(fit.step_locations);
    puts("Error (re)allocating memory");
    exit(1);
  }
}

int* UpdateStepIndices(step_fit fit,
		       int n_dwells,
		       int* step_indices) {
  int n_steps = n_dwells-1;
  for(int i=0;i<n_steps;i++) {
    step_indices[fit.step_locations[i]] = 1;
  }
  return step_indices;
}

double GetPreviousDwellMean(step_fit prev_fit,
			    int* step_indices,
			    int index) {
  int loc=0; 
  for(int i=0;i<index;i++) {
    loc+=step_indices[i];
  }
  return prev_fit.means[loc];
}


