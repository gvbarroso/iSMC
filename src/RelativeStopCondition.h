/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 15/04/2018
 *
 */

#ifndef _RELATIVESTOPCONDITION_H_
#define _RELATIVESTOPCONDITION_H_

#include <string>
#include <vector>
#include <limits>
#include <math.h>

#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/OptimizationStopCondition.h>
#include <Bpp/Numeric/NumTools.h>


class RelativeStopCondition:
  public bpp::FunctionStopCondition { 
private:
  mutable double parameterValueTolerance_;
  mutable bpp::ParameterList newParameters_;
  mutable bpp::ParameterList lastParameters_;
  mutable bool reachedFunctionTolerance_;
  mutable std::vector< bool > reachedParametersTolerance_;  
  mutable double lastFunctionValue_;
  mutable double newFunctionValue_;

public:
  RelativeStopCondition(const bpp::OptimizerInterface* optimizer,
                        double likelihoodTolerance,
                        double paramValueTolerance):
  bpp::FunctionStopCondition(optimizer, likelihoodTolerance),
  parameterValueTolerance_(paramValueTolerance),
  newParameters_(optimizer -> getParameters()),
  lastParameters_(optimizer -> getParameters()),
  reachedFunctionTolerance_(false),
  reachedParametersTolerance_(optimizer -> getParameters().size()),
  lastFunctionValue_(-log(0.)),
  newFunctionValue_(-log(0.))
  {
    init();
    for(size_t i = 0; i < reachedParametersTolerance_.size(); ++i) {
      reachedParametersTolerance_[i] = false;
    }
  }
       
  RelativeStopCondition* clone() const { return new RelativeStopCondition(*this); }
  
  void init() {
    FunctionStopCondition::init();
    newFunctionValue_ = -log(0.);
    if(optimizer_ -> getFunction() != 0) {
      newFunctionValue_ = optimizer_ -> getFunctionValue();
    }
  }
  
  bool isToleranceReached() const { //combines conditions from parameter values and the function being optimized (likelihood)
    bool reachedAllParameters = checkParametersTolerance_();
    bool reachedFunction = checkFunctionTolerance_();
    return(reachedFunction && reachedAllParameters); 
  }
  
  double getCurrentTolerance() const;
  
private: 
  bool checkFunctionTolerance_() const;
  
  bool checkParametersTolerance_() const;
  
};

#endif
