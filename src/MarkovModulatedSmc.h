/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 31/12/2018
 *
 */


#ifndef _MARKOVMODULATEDSMC_H_
#define _MARKOVMODULATEDSMC_H_

#include <string>
#include <vector>

#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Text/TextTools.h>

#include "SequentiallyMarkovCoalescent.h"
#include "HmmStatesLibrary.h"
#include "ParameterCategoryTransitions.h"
#include "OptionsContainer.h"

class MarkovModulatedSmc:
  public SequentiallyMarkovCoalescent {
private:
  std::shared_ptr< HmmStatesLibrary > hmmSates_;
  std::vector<std::shared_ptr<bpp::DiscreteDistributionInterface>> paramScalings_;
  std::vector<std::shared_ptr<ParameterCategoryTransitions>> categoryTransitions_;

public:    
  MarkovModulatedSmc(std::shared_ptr<OptionsContainer> smcOptions,
                     const std::vector<std::vector<unsigned char>>& snpCallsFromAllDataSets,
                     std::shared_ptr<HmmStatesLibrary> hmmSates,
                     std::vector<std::shared_ptr<bpp::DiscreteDistributionInterface>> paramScalings,
                     std::vector<std::shared_ptr<ParameterCategoryTransitions>> categoryTransitions, 
                     const ParameterAlphabet& parameterAlphabet): 
  SequentiallyMarkovCoalescent(smcOptions, snpCallsFromAllDataSets),
  hmmSates_(hmmSates),
  paramScalings_(paramScalings),
  categoryTransitions_(categoryTransitions)
  { }
  MarkovModulatedSmc(unsigned int numIntervals,
                     const std::string& timeDisc,
                     double tMax,
                     const std::vector< std::vector < unsigned char > >& snpCallsFromAllDataSets,
                     std::shared_ptr< HmmStatesLibrary > hmmSates,
                     std::vector< std::shared_ptr< bpp::DiscreteDistributionInterface > > paramScalings,
                     std::vector< std::shared_ptr< ParameterCategoryTransitions > > categoryTransitions, 
                     const ParameterAlphabet& parameterAlphabet): 
  SequentiallyMarkovCoalescent(numIntervals, timeDisc, tMax, snpCallsFromAllDataSets),
  hmmSates_(hmmSates),
  paramScalings_(paramScalings),
  categoryTransitions_(categoryTransitions)
  { }
  MarkovModulatedSmc(unsigned int numIntervals,
                     const std::string& timeDisc,
                     double tMax,
                     const bpp::ParameterList& rho_theta, 
                     const bpp::ParameterList& lambdas,
                     std::shared_ptr< HmmStatesLibrary > hmmSates,
                     std::vector< std::shared_ptr< bpp::DiscreteDistributionInterface > > paramScalings,
                     std::vector< std::shared_ptr< ParameterCategoryTransitions > > categoryTransitions, 
                     const ParameterAlphabet& parameterAlphabet): 
  SequentiallyMarkovCoalescent(numIntervals, timeDisc, tMax, rho_theta, lambdas),
  hmmSates_(hmmSates),
  paramScalings_(paramScalings),
  categoryTransitions_(categoryTransitions)
  { }
  MarkovModulatedSmc(unsigned int numIntervals,
                     const std::string& timeDisc,
                     double tMax,
                     const bpp::ParameterList& optimParams, //all optimised parameters
                     std::shared_ptr< HmmStatesLibrary > hmmSates,
                     std::vector< std::shared_ptr< bpp::DiscreteDistributionInterface > > paramScalings,
                     std::vector< std::shared_ptr< ParameterCategoryTransitions > > categoryTransitions, 
                     const ParameterAlphabet& parameterAlphabet): 
  SequentiallyMarkovCoalescent(numIntervals, timeDisc, tMax, optimParams),
  hmmSates_(hmmSates),
  paramScalings_(paramScalings),
  categoryTransitions_(categoryTransitions)
  {
    for(size_t i = 0; i < paramScalings_.size(); ++i) {
      if(optimParams.getCommonParametersWith(paramScalings_[i] -> getIndependentParameters()).size() > 0) {
        paramScalings_[i] -> matchParametersValues(optimParams);
        paramScalings_[i] -> discretize(); 
      }
      if(optimParams.getCommonParametersWith(categoryTransitions_[i] -> getIndependentParameters()).size() > 0) {
        categoryTransitions_[i] -> matchParametersValues(optimParams);
        categoryTransitions_[i] -> setUpCategoryTransitionMatrix();
      }
    }
  }
  
  MarkovModulatedSmc* clone() const { return new MarkovModulatedSmc(*this); } 

public:
  std::vector< std::shared_ptr< bpp::DiscreteDistributionInterface > >& getParameterScalings() {
    return paramScalings_;
  }
  const std::vector< std::shared_ptr< bpp::DiscreteDistributionInterface > >& getParameterScalings() const {
    return paramScalings_;
  }
  
  void setParameterScalings(std::shared_ptr< bpp::DiscreteDistributionInterface > newScaling, size_t index) {
    paramScalings_[index] = newScaling;  
  }
  
  std::vector< std::shared_ptr< ParameterCategoryTransitions > >& getParameterTransitions() {
    return categoryTransitions_;
  }
  const std::vector< std::shared_ptr< ParameterCategoryTransitions > >& getParameterTransitions() const {
    return categoryTransitions_;
  }
  
  void setParameterTransitions(std::shared_ptr< ParameterCategoryTransitions > newTrans, size_t index) {
    categoryTransitions_[index] = newTrans;
  }
  
  std::shared_ptr< HmmStatesLibrary >& getHmmStatesLibrary() { 
    return hmmSates_;
  }

  const std::shared_ptr< HmmStatesLibrary >& getHmmStatesLibrary() const { 
    return hmmSates_;
  }
  
  size_t getNumberOfHiddenStates() {
    return hmmSates_ -> getNumberOfHiddenStates();
  }
  
  size_t getNumberOfModulatedParameters() {
    return hmmSates_ -> getNumberOfModulatedParams();
  }

  size_t countNumberOfModulatedParameters() {
    return static_cast< size_t >(std::count_if(std::begin(categoryTransitions_), std::end(categoryTransitions_), isParamHeterogeneous));
  }

  static bool isParamHeterogeneous(const std::shared_ptr< ParameterCategoryTransitions >& paramTrans) { 
    return paramTrans -> getNumberOfCategories() > 1;
  }
  
  void scaleCoalescenceRates(size_t NeCategory);
  
};

#endif
