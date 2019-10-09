/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 21/07/2018
 * Last modified: 04/07/2019
 *
 */


#ifndef _SPLINESMODEL_H_
#define _SPLINESMODEL_H_

#include <Bpp/Numeric/Function/Functions.h>
#include "Splines.h"


class SplinesModel:
  public Splines,
  public bpp::Function {   
protected:
  //bpp::ParameterList params_; //
  std::shared_ptr< MmSmcEmissionProbabilities > mmsmcep_;
  std::shared_ptr< MmSmcTransitionProbabilities > mmsmctp_;
  std::shared_ptr< MultipleMmPsmc > mPsmc_;
  double logLikelihood_;
  double aic_;
  
public:
  SplinesModel(std::shared_ptr< MarkovModulatedSmc > mmsmc,
               std::shared_ptr< MmSmcEmissionProbabilities > mmsmcep,
               std::shared_ptr< MmSmcTransitionProbabilities > mmsmctp,
               std::shared_ptr< MultipleMmPsmc > mPsmc,
               const bpp::ParameterList& params,
               size_t numberOfKnots,
               const std::string& splinesType): 
  Splines(mmsmc, numberOfKnots, splinesType),
  mmsmcep_(mmsmcep),
  mmsmctp_(mmsmctp),
  mPsmc_(mPsmc),
  logLikelihood_(-1.),
  aic_(-1.)
  {
    includeParameters_(params);
  }
  
  SplinesModel(std::shared_ptr< MarkovModulatedSmc > mmsmc,
               std::shared_ptr< MmSmcEmissionProbabilities > mmsmcep,
               std::shared_ptr< MmSmcTransitionProbabilities > mmsmctp,
               std::shared_ptr< MultipleMmPsmc > mPsmc,
               size_t numberOfKnots,
               const std::string& splinesType): 
  Splines(mmsmc, numberOfKnots, splinesType),
  mmsmcep_(mmsmcep),
  mmsmctp_(mmsmctp),
  mPsmc_(mPsmc),
  logLikelihood_(-1.),
  aic_(-1.)
  { }

  SplinesModel* clone() const { return new SplinesModel(*this); }

  void setParameters(const bpp::ParameterList& parameters) {
    SplinesModel::setParametersValues(parameters);
  }

  double getValue() const { return - logLikelihood_; }
  
  void fireParameterChanged(const bpp::ParameterList& parameters);
  
  double getLogLikelihood() { 
    return logLikelihood_;
  }
  
  void computeAic() { 
    double numberOfParameters = static_cast< double >(getNumberOfIndependentParameters());
    aic_ = 2. * numberOfParameters - 2. * logLikelihood_;
  }

  double getAic() { 
    return aic_;
  }

  //returns SMC-related parameters (rho, theta, rho.alpha, r_ii...);
  //useful to fit a model enforcing flat demography (e.g., to test the effect of neglecting demography)
  bpp::ParameterList fetchNonSplinesParameters();

  bpp::ParameterList fetchModelParameters(); 

};

#endif




