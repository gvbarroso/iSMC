/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 16/01/2019
 *
 */


#include "ParameterCategoryTransitions.h"

using namespace std;
using namespace bpp;

//#define DEBUG_CT


void ParameterCategoryTransitions::setUpCategoryTransitionMatrix() { 
    
  //NOTE: numberOfCategories_ has information about the model (ie, it's not just gamma categories)
    
  if(hetRateModel_ == "Hotspot") {
    parameterTransitionMatrix_[0][0] = 1. - getParameterValue(parameterPrefix_ + "_backHot");
    parameterTransitionMatrix_[0][1] = getParameterValue(parameterPrefix_ + "_backHot");
    parameterTransitionMatrix_[1][0] = getParameterValue(parameterPrefix_ + "_hotBack");
    parameterTransitionMatrix_[1][1] = 1. - getParameterValue(parameterPrefix_ + "_hotBack");
  }
  else if(hetRateModel_ == "Gamma") {
    double selfTransitionProb = getParameterValue(parameterPrefix_ + "_ii");
    for(size_t i = 0; i < numberOfCategories_; ++i) {
      for(size_t j = 0; j < numberOfCategories_; ++j) {
        if(i == j) { 
          parameterTransitionMatrix_[i][j] = selfTransitionProb;
        }
        else { 
          parameterTransitionMatrix_[i][j] = (1. - selfTransitionProb) / static_cast< double >(numberOfCategories_ - 1);
        }
      }
    }
  }
  else if(hetRateModel_ == "Gamma+Hotspot") { 
    //because of constraints in the full Markov Chain, we need to use reparametrisation when building the TM  
    double kappa = getParameterValue(parameterPrefix_ + "_k");
    //cout << "kappa = " << kappa << "; "; 
    double hotBack = getParameterValue(parameterPrefix_ + "_v");
    //cout << "hotBack = " << hotBack << "; "; 
    double delta = getParameterValue(parameterPrefix_ + "_d"); 
    //cout << "delta = " << delta << "; " << endl << endl; 
    //0 = hotspot; 1:(numberOfCategories_) = gamma categories
    for(size_t i = 0; i < numberOfCategories_; ++i) {
      for(size_t j = 0; j < numberOfCategories_; ++j) {
        if(i == 0) {
          if(j == 0) { //hotspot -> hotspot
            parameterTransitionMatrix_[i][j] = 1. - hotBack; 
            //cout << parameterTransitionMatrix_[i][j] << " ";
          }
          else { //hotspot -> background category
            parameterTransitionMatrix_[i][j] = hotBack / static_cast< double >(numberOfCategories_ - 1);
            //cout << parameterTransitionMatrix_[i][j] << " ";
          }  
        }
        else { 
          if(j == 0) { //background category -> hotspot 
            parameterTransitionMatrix_[i][j] = kappa * delta;
            //cout << parameterTransitionMatrix_[i][j] << " ";
          }
          else { //background category -> same background category
            if(i == j) { //same background categories
              parameterTransitionMatrix_[i][j] = 1. - delta;
              //cout << parameterTransitionMatrix_[i][j] << " ";
            }
            else { //different background categories
              parameterTransitionMatrix_[i][j] = (1. - kappa) * delta / static_cast< double >(numberOfCategories_ - 2);
              //cout << parameterTransitionMatrix_[i][j] << " ";
            }
          }
        }         
      }
      //cout << endl;
    }
    //cout << endl;
  }
  else {
    throw Exception("Mis-specified model of spatial rate heterogeneity");
  }
  #ifdef DEBUG_CT
  cout << endl << "computing param. TM row sums..." << endl;
  for(size_t i = 0; i < parameterTransitionMatrix_.size(); ++i) {
    Vdouble temp = parameterTransitionMatrix_[i];
    sort(temp.begin(), temp.end());
    double rowSum = VectorTools::sum(temp);
    if(!(abs(1. - rowSum) < 1e-9)) {
      cout << "Sum of entries for row " << i << " = " << setprecision(20) << rowSum << endl;
    }
  }
  cout << "   done." << endl;
  #endif
}
  
void ParameterCategoryTransitions::includeCategoryTransitionParameters_() {
  ParameterList categoryTransitions;
  if(hetRateModel_ == "Hotspot") {
    double backHot = 2e-6;
    string transitionIntoHotspots = parameterPrefix_ + "_backHot";
    categoryTransitions.addParameter(new Parameter(transitionIntoHotspots, backHot, make_shared<IntervalConstraint>(0., 1., true, true)));
    double hotBack = 5e-4;
    string transitionIntoBackground = parameterPrefix_ + "_hotBack";
    categoryTransitions.addParameter(new Parameter(transitionIntoBackground, hotBack, make_shared<IntervalConstraint>(0., 1., true, true)));
  }
  else if(hetRateModel_ == "Gamma") {
    double selfTransitionProb = 1. - 1e-5;
    if(numberOfCategories_ == 1) {
      selfTransitionProb = 1.;
    }    
    string paramName = parameterPrefix_ + "_ii";
    categoryTransitions.addParameter(new Parameter(paramName, selfTransitionProb, make_shared<IntervalConstraint>(0., 1., true, true)));
  }
  else if(hetRateModel_ == "Gamma+Hotspot") {
    string kappa = parameterPrefix_ + "_k"; //transition into hotspots
    double kappaVal = 2e-6;
    categoryTransitions.addParameter(new Parameter(kappa, kappaVal, make_shared<IntervalConstraint>(0., 1., true, true)));
    //transition into gamma and transition between gamma categories have been re-parametrized 
    string transitionIntoBackground = parameterPrefix_ + "_v"; 
    double hotBack = 5e-4;
    categoryTransitions.addParameter(new Parameter(transitionIntoBackground, hotBack, make_shared<IntervalConstraint>(0., 1., true, true)));
    string delta = parameterPrefix_ + "_d";
    double deltaVal = .9999;
    categoryTransitions.addParameter(new Parameter(delta, deltaVal, make_shared<IntervalConstraint>(0., 1., true, true)));
  }
  addParameters_(categoryTransitions);
}
