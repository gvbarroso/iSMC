/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 01/01/2019
 *
 */


#ifndef _SMCTRANSITIONPROBABILITIES_H_
#define _SMCTRANSITIONPROBABILITIES_H_

#include <limits>
#include <math.h>

#include <Bpp/Numeric/ParameterList.h>

#include "AbstractParametrizedMatrix.h"
#include "SequentiallyMarkovCoalescent.h"


class SmcTransitionProbabilities:
  public AbstractParametrizedMatrix {
private:
  std::shared_ptr< SequentiallyMarkovCoalescent > smc_;

public:
  SmcTransitionProbabilities(std::shared_ptr< SequentiallyMarkovCoalescent > smc):
  AbstractParametrizedMatrix(),
  smc_(smc)
  { 
    expectedMatrix_.resize(smc -> getNumberOfHiddenStates(), bpp::Vdouble(smc -> getNumberOfHiddenStates()));
    setUpExpectedMatrix();
  }
 
public:
  void setUpExpectedMatrix();
  
protected:
  //see Schiffels & Durbin (2014)
  double exponentiatedIntegral_(double time1, double time2, const bpp::ParameterList& lambdaVec);
  
  //see Schiffels & Durbin (2014)
  double exponentiatedIntegralM_(double time1, double time2, const bpp::ParameterList& lambdaVec);
  
  //see Schiffels & Durbin (2014) 
  double computeTransitionProbabilityDecreasedTime_(double timeBeta, double timeAlpha, const bpp::ParameterList& lambVec, double rho);
  
  //see Schiffels & Durbin (2014) 
  double computeTransitionProbabilityIncreasedTime_(double timeBeta, double timeAlpha, const bpp::ParameterList& lambVec, double rho);

};

#endif