#ifndef SICALGORITHM_H_INCLUDED
#define SICALGORITHM_H_INCLUDED

const double Pi = 3.1415926535897;

struct step_fit {
  double* means; //list of means of dwells
  int* step_locations;
  double SIC; //SIC value for this fit
  double chisq; //chisq value for this fit
};

step_fit NoStepSIC(double*,
		   int);

step_fit AddStepSIC(double*, 
		    int, 
		    int*,
		    int,
		    step_fit);

step_fit InitializeFitToZeros(step_fit,
			      int);

step_fit SetFitToZeros(step_fit,
		       int);

step_fit MakeFitOneDwellBigger(step_fit,
			   int);

int* UpdateStepIndices(step_fit,
		       int,
		       int*);

double GetPreviousDwellMean(step_fit,
			    int*,
			    int);

#endif
