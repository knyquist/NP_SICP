/* Kristofor Nyquist 
Step fitting algorithm in C using bias on expected noise.
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "InputTrace.h"
#include "SICalgorithm.h"

int main(int argc, char *argv[] )
{
  /* bias parameters */

  /* Get the trace information */
  int numberofpoints = GetTraceLength(argv[1]);
  double** arr = InputTrace(argv[1], numberofpoints);
  double* time;     time     = (double*) malloc(numberofpoints*sizeof(double));
  double* force;    force    = (double*) malloc(numberofpoints*sizeof(double));
  double* position; position = (double*) malloc(numberofpoints*sizeof(double));
  for(int i=0;i<numberofpoints;i++) { //dump data to time and position arrays
    time[i] = arr[0][i];
    force[i] = arr[1][i];
    position[i] = arr[2][i];
  }

  /* Get the biasing parameters from parameter file 
     (this is supplied by argv[2])                  */
  double nu = 0.; double So = 0.; double bparams[2] = {}; //nu first, then So
  char bnames[2]={};
  FILE * fi; fi=fopen(argv[2],"r");
  int i=0; while(i<4) { //there are only two params, nu and So
    fscanf(fi,"%s %lf", &bnames[i], &bparams[i]); i+=1;
  } fclose(fi);
  if(bnames[0]=='n' && bnames[1]=='S') {
    nu=bparams[0]; So=bparams[1]; } 
  else { return(0); printf("Error: biasing parameters could not be properly initialized\n"); }
 
  // Initialize with No Step SIC calculation 
  step_fit this_fit; this_fit = InitializeFitToZeros(this_fit,1);
  this_fit = NoStepSIC(position,
		       numberofpoints,
		       nu,
		       So);
  
  // Start Adding Steps, terminate when SIC no longer minimized 
  int n_dwells = 1; step_fit final_fit;
  step_fit prev_fit; prev_fit = InitializeFitToZeros(prev_fit,n_dwells);
  int* step_indices; step_indices = (int*) malloc(numberofpoints*sizeof(int));  
  for(int i=0;i<numberofpoints;i++) step_indices[i] = 0;
  do {
    // the last iteration's fit becomes the previous fit
    prev_fit = InitializeFitToZeros(prev_fit,n_dwells-1);
    prev_fit = MakeFitOneDwellBigger(prev_fit,n_dwells);
    for(int i=0;i<n_dwells;i++) {
      prev_fit.means[i]=this_fit.means[i];
      if(i!=n_dwells) prev_fit.step_locations[i] = this_fit.step_locations[i];
      }
    prev_fit.SIC = this_fit.SIC;
    prev_fit.chisq = this_fit.chisq;
    step_indices = UpdateStepIndices(prev_fit,n_dwells,step_indices);
    // resize this_fit's lists to allow space for the next fit
    this_fit = MakeFitOneDwellBigger(this_fit,n_dwells+1);
    this_fit = SetFitToZeros(this_fit,n_dwells+1); //populate with zeros
    this_fit = AddStepSIC(position,
			  numberofpoints,
			  step_indices,
			  n_dwells,
			  prev_fit,
			  nu,
			  So);
    n_dwells+=1;
  }
  while(this_fit.SIC < prev_fit.SIC);
  //printf("%f\n",prev_fit.chisq);
  //output the fit!
  printf("%f %f\n", time[0], prev_fit.means[0]);
  for(int i=0;i<n_dwells-2;i++) {
    printf("%f %f\n",time[prev_fit.step_locations[i]],prev_fit.means[i]);
    printf("%f %f\n",time[prev_fit.step_locations[i]],prev_fit.means[i+1]);
  }
  printf("%f %f\n",time[numberofpoints-1],prev_fit.means[n_dwells-2]);
  return 0;
}


