/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 21/07/2018
 * Last modified: 25/02/2020
 *
 */


#include "SplinesModel.h"

using namespace std;
using namespace bpp;


void SplinesModel::fireParameterChanged(const ParameterList& parameters) {
    
  for(size_t i = 0; i < mmsmc_ -> getParameterScalings().size(); ++i) {
      
    if(parameters.getCommonParametersWith(mmsmc_ -> getParameterScalings()[i] -> getIndependentParameters()).size() > 0) {
      mmsmc_ -> getParameterScalings()[i] -> matchParametersValues(parameters);
      mmsmc_ -> getParameterScalings()[i] -> discretize(); 
    }
    
    if(parameters.getCommonParametersWith(mmsmc_ -> getParameterTransitions()[i] -> getIndependentParameters()).size() > 0) {
      mmsmc_ -> getParameterTransitions()[i] -> matchParametersValues(parameters);
      mmsmc_ -> getParameterTransitions()[i] -> setUpCategoryTransitionMatrix();
    }
  }
  
  if(parameters.getCommonParametersWith(mmsmc_ -> getParameters()).size() > 0) {
    mmsmc_ -> matchParametersValues(parameters);
  }
  
  if(parameters.getCommonParametersWith(getParameters()).size() > 0) {
    mapLambdasFromSplines(mmsmc_ -> getLambdaVector()); 
  }
  
  mmsmc_ -> computeAverageCoalescenceTime();
  mmsmcep_ -> setUpExpectedMatrix();
  mmsmctp_ -> setUpExpectedMatrix();
  
  logLikelihood_ =  mPsmc_ -> fetchCompositeLogLikelihood();
  
  computeAic();
}

ParameterList SplinesModel::fetchNonSplinesParameters() {
    
  ParameterList params(getParameters());
  
  //we don't pass splines parameters to enforce flat demography
  for(size_t i = 0; i <= (numberOfKnots_ + 1); ++i) {
    string nameIntercept = "y" + TextTools::toString(i);
    string nameFirstDerivative = "y" + TextTools::toString(i) + "_prime";
    params.deleteParameter(nameIntercept);
    params.deleteParameter(nameFirstDerivative);
  }
  
  return params;
}

ParameterList SplinesModel::fetchModelParameters() { 
  
  ParameterList params(getParameters());
  
  if(splinesType_ == "Sigmoidal") {
    //we don't pass the derivatives to the optimizer so that they stay fixed to 0.
    for(size_t i = 0; i <= (numberOfKnots_ + 1); ++i) {
      string nameFirstDerivative = "y" + TextTools::toString(i) + "_prime";
      params.deleteParameter(nameFirstDerivative);
    }
  }
  
  for(size_t i = 0; i < ignoreParams_.size(); ++i)
    params.deleteParameter(ignoreParams_[i]);
      
  return params;
}
