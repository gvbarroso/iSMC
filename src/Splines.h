/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 03/07/2019
 *
 */


#ifndef _SPLINES_H_
#define _SPLINES_H_

#include <string>
#include <vector>

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/ParameterList.h>

#include "MarkovModulatedSmc.h"
#include "MmSmcEmissionProbabilities.h"
#include "MmSmcTransitionProbabilities.h"
#include "MultipleMmPsmc.h"


class Splines:
  public bpp::AbstractParameterAliasable {   
      
protected:
  std::shared_ptr< MarkovModulatedSmc > mmsmc_;
  size_t numberOfKnots_;
  bpp::Vdouble fixedKnots_;
  std::string splinesType_; 
  
public:
  Splines(std::shared_ptr< MarkovModulatedSmc > mmsmc,
          size_t numberOfKnots,
          const std::string& splinesType): 
  AbstractParameterAliasable(""),
  mmsmc_(mmsmc),
  numberOfKnots_(numberOfKnots),
  fixedKnots_(0),
  splinesType_(splinesType)
  {
    setUpSplinesKnots_();
    includeSplinesParameters_();
  }
    
public:
  Splines* clone() const { return new Splines(*this); }

  void setParameters(const bpp::ParameterList& parameters) {
    AbstractParameterAliasable::setParametersValues(parameters);
  }
    
  void fireParameterChanged(const bpp::ParameterList& parameters);
  
  const bpp::Vdouble& getSplinesKnots() const {
    return fixedKnots_;
  }
  bpp::Vdouble& getSplinesKnots() {
    return fixedKnots_;
  }  
  
  size_t getNumberOfKnots() {
    return numberOfKnots_;
  }
  
  void mapLambdasFromSplines(bpp::ParameterList& Lambdas);
    
private:
  void setUpSplinesKnots_() {    
      
    size_t span = (mmsmc_ -> getNumberOfIntervals() + numberOfKnots_)  / (numberOfKnots_ + 1);
    
    for(size_t i = 0; i < mmsmc_ -> getNumberOfIntervals() - 1; ++i) {
      if((i + 1) % span == 0) { 
        fixedKnots_.push_back(mmsmc_ -> getTimeIntervals()[i]);
      }
    }
  }
  
  void includeSplinesParameters_();   
  
  double getPowerOne_(double leftKnot, double rightKnot, size_t index,
                      double leftIntercept, double rightIntercept,
                      double leftDerivative, double rightDerivative);
  
  double getPowerTwo_(double leftKnot, double rightKnot, size_t index,
                      double leftIntercept, double rightIntercept,
                      double leftDerivative, double rightDerivative);
  
  double getPowerThree_(double leftKnot, double rightKnot, size_t index,
                        double leftIntercept, double rightIntercept,
                        double leftDerivative, double rightDerivative);
  
  double getConstant_(double leftKnot, double rightKnot, size_t index,
                      double leftIntercept, double rightIntercept,
                      double leftDerivative, double rightDerivative);
  
};

#endif
