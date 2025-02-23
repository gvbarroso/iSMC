/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 03/07/2019
 *
 */


#ifndef _SMCOPTIMIZATIONWRAPPER_H_
#define _SMCOPTIMIZATIONWRAPPER_H_

#include <string>
#include <vector>
#include <limits>

#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Phyl/PseudoNewtonOptimizer.h>

#include "RelativeStopCondition.h"
#include "BaumWelch.h"
#include "SplinesModel.h"
#include "MarkovModulatedSmc.h"
#include "MmSmcEmissionProbabilities.h"
#include "MmSmcTransitionProbabilities.h"
#include "MultipleMmPsmc.h"
#include "OptionsContainer.h"


class SmcOptimizationWrapper {
    
private:
    
  std::shared_ptr< MarkovModulatedSmc > mmsmc_;
  std::shared_ptr< MmSmcEmissionProbabilities > mmsmcep_;
  std::shared_ptr< MmSmcTransitionProbabilities > mmsmctp_;
  std::shared_ptr< MultipleMmPsmc > mPsmc_;
  std::shared_ptr< OptionsContainer > smcOptions_;
  std::vector< std::shared_ptr < SplinesModel > > listOfModels_;
  bpp::VVdouble listOfTestLikelihoods_; //for LOOCV
  bpp::ParameterList bestParameters_;
  double bestAic_;
  
public:
  SmcOptimizationWrapper(std::shared_ptr< MarkovModulatedSmc > mmsmc,
                         std::shared_ptr< MmSmcEmissionProbabilities > mmsmcep,
                         std::shared_ptr< MmSmcTransitionProbabilities > mmsmctp,
                         std::shared_ptr< MultipleMmPsmc > mPsmc,
                         std::shared_ptr< OptionsContainer > smcOptions):
  mmsmc_(mmsmc),
  mmsmcep_(mmsmcep),
  mmsmctp_(mmsmctp),
  mPsmc_(mPsmc),
  smcOptions_(smcOptions),
  listOfModels_(0),
  listOfTestLikelihoods_(0, bpp::Vdouble(0)),
  bestParameters_(),
  bestAic_(std::numeric_limits<double>::max()) 
  {
    //standard SMC parameters
    bestParameters_.addParameters(mmsmc -> getParameters());
    bestParameters_.addParameters(mmsmc -> getLambdaVector());
    
    //add Markov-modulation parameters
    for(size_t i = 0; i < mmsmc -> getParameterScalings().size(); ++i) {
        
      //if hotspot model, only bring hotspot intensity to optimization (ie discard PMF probs.)
      if(mmsmc -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Hotspot") {
        bestParameters_.addParameter(mmsmc -> getParameterScalings()[i] -> parameter("V2"));  
      }
      
      else if(mmsmc -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Gamma") {
          
        if(mmsmc -> getParameterScalings()[i] -> getNumberOfCategories() > 1) {
          bestParameters_.addParameters(mmsmc -> getParameterScalings()[i] -> getIndependentParameters()); 
        }
      }
      
      else if(mmsmc -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Gamma+Hotspot") {
          
        if(mmsmc -> getParameterScalings()[i] -> getNumberOfCategories() > 1) {
          bestParameters_.addParameters(mmsmc -> getParameterScalings()[i] -> getIndependentParameters()); 
        }
        
        else {
          bestParameters_.addParameter(mmsmc -> getParameterScalings()[i] -> parameter("heat"));
        }
      }
    }
    
    for(size_t i = 0; i < mmsmc -> getParameterTransitions().size(); ++i) {
        
      if(mmsmc -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Gamma") {
          
        if(mmsmc -> getParameterScalings()[i] -> getNumberOfCategories() > 1) {
          bestParameters_.addParameters(mmsmc -> getParameterTransitions()[i] -> getParameters());
        }
      }
      
      else {
        bestParameters_.addParameters(mmsmc -> getParameterTransitions()[i] -> getParameters());
      }
    }
  }
  
public:
    
  bpp::ParameterList& getBestParameters() { 
    return bestParameters_;
  }
  const bpp::ParameterList& getBestParameters() const { 
    return bestParameters_;
  }
    
  std::vector< std::shared_ptr < SplinesModel > >& getListOfModels() {
    return listOfModels_;
  }
  const std::vector< std::shared_ptr < SplinesModel > >& getListOfModels() const {
    return listOfModels_;
  }
  
  const bpp::VVdouble& getListOfTestLikelihoods() const {
    return listOfTestLikelihoods_;
  }
  
  double getAic() { 
    return bestAic_;
  }
  
  std::shared_ptr< SplinesModel > selectBestModel();

  void stepwiseExpectationMaximization();
  
  void optimizeParameters();
  
  //to resume optimisation after a problem:
  void optimizeParameters(const bpp::ParameterList& backupParams);
  
  void writeEstimatesToFile(const bpp::ParameterList& params, double AIC);

  void writeDemographyToFile(double theta);

private:
  void fireUpdateBestValues_(SplinesModel* bestSm, const bpp::ParameterList &params);
  
  void createAndFitSplinesModels_(bpp::ParameterList& nonSplinesParams);
    
  void fitModel_(std::shared_ptr<SplinesModel> sm);
  
};

#endif
