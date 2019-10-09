/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 18/04/2018
 *
 */


#include "RelativeStopCondition.h"
#include <algorithm>

using namespace std;
using namespace bpp;

  
double RelativeStopCondition::getCurrentTolerance() const {
  if(callCount_ > burnin_) {
    return NumTools::abs< double >((newFunctionValue_ - lastFunctionValue_) / lastFunctionValue_);
  }
  else {
    return max(tolerance_, 1.);
  }
}
  
bool RelativeStopCondition::checkFunctionTolerance_() const {
  if(reachedFunctionTolerance_) {
   return true; //once tolerance_ is reached, it will always return true in later steps  
  }
  else {
    callCount_++;
    lastFunctionValue_ = newFunctionValue_;
    newFunctionValue_ = optimizer_ -> getFunctionValue();
    if(callCount_ <= burnin_) {
      return false;
    }
    //uses a relative improvement stop condition
    double tol = NumTools::abs< double >((newFunctionValue_ - lastFunctionValue_) / lastFunctionValue_);
    if(tol < tolerance_ && !reachedFunctionTolerance_) {
      //cout << endl << "Reached function numerical tolerance!" << endl;
      reachedFunctionTolerance_ = true;
    }
  }
  return reachedFunctionTolerance_;
}
  
bool RelativeStopCondition::checkParametersTolerance_() const {
  lastParameters_.matchParametersValues(newParameters_);  
  newParameters_.matchParametersValues(optimizer_ -> getParameters());
  vector< string > parameterNames(newParameters_.getParameterNames());
  for(size_t i = 0; i < newParameters_.size(); ++i) {
    double newVal = newParameters_.getParameterValue(parameterNames[i]);
    double lastVal = lastParameters_.getParameterValue(parameterNames[i]);
    if(lastVal == 0) {  //avoids dividing by zero
      lastVal = 1e-6;
    }
    double tol = NumTools::abs< double >((newVal - lastVal) / lastVal);
    if(tol < parameterValueTolerance_) {
      if(!reachedParametersTolerance_[i]) {
        reachedParametersTolerance_[i] = true;  
        if(parameterValueTolerance_ < numeric_limits< double >::max()) { 
          cout << endl << "Numerical tolerance reached for parameter: " << parameterNames[i] << endl;    
        }
      }
      //we can choose to remove the parameter from the lists (new, last and "freeze" the optimizer?) 
    }
  }
  //returns true IFF all parameters have reached their tolerance
  return count(reachedParametersTolerance_.begin(), reachedParametersTolerance_.end(), false) < 1;
}
