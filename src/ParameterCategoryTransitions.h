/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 01/01/2019
 *
 */

#ifndef _PARAMETERCATEGORYTRANSITIONS_H_
#define _PARAMETERCATEGORYTRANSITIONS_H_

#include <string>
#include <vector>

#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/VectorTools.h>

class ParameterCategoryTransitions:
  public bpp::AbstractParameterAliasable { 
protected:
  bpp::VVdouble parameterTransitionMatrix_;
  bpp::Vdouble paramInitProbs_;
  unsigned int numberOfCategories_; //NOTE has info. about any model (gamma, hotspot, gamma+hotspot)
  std::string parameterPrefix_;
  std::string hetRateModel_;
  
public:
  ParameterCategoryTransitions(unsigned int numGammaCategories, // >= 1
                               const std::string& hetRateModel, 
                               const std::string& spatialParamPrefix):
  AbstractParameterAliasable(""),
  parameterTransitionMatrix_(0, bpp::Vdouble(0)),
  paramInitProbs_(0),
  numberOfCategories_(numGammaCategories),
  parameterPrefix_(spatialParamPrefix),
  hetRateModel_(hetRateModel)
  {
    //updates numberOfCategories_ based on het. model
    if(hetRateModel == "Gamma+Hotspot") {
      ++numberOfCategories_; //adds category for hotspot heat
    }
    else if(hetRateModel == "Hotspot") {
      numberOfCategories_ = 2; //background, hotspot
    }

    paramInitProbs_.resize(numberOfCategories_);
    for(size_t i = 0; i < numberOfCategories_; ++i) {
      paramInitProbs_[i] = 1. / static_cast< double >(numberOfCategories_);
    }

    parameterTransitionMatrix_.resize(numberOfCategories_);
    for(size_t i = 0; i < parameterTransitionMatrix_.size(); ++i) {
      parameterTransitionMatrix_[i].resize(numberOfCategories_);  
    }

    includeCategoryTransitionParameters_();
    setUpCategoryTransitionMatrix();
  }
    
  ParameterCategoryTransitions* clone() const { return new ParameterCategoryTransitions(*this); }

  virtual ~ParameterCategoryTransitions() {}

public:
  const bpp::VVdouble& getCategoryTransitions() const { 
    return parameterTransitionMatrix_;
  }
  bpp::VVdouble& getCategoryTransitions() { 
    return parameterTransitionMatrix_;
  }
  
  const bpp::Vdouble& getInitializationProbabilities() const { 
    return paramInitProbs_;
  }
  bpp::Vdouble& getInitializationProbabilities() { 
    return paramInitProbs_;
  }
 
  unsigned int getNumberOfCategories() { 
    return numberOfCategories_;
  }
    
  const std::string& getSpatialParameterPrefix() const { 
    return parameterPrefix_;
  }

  const std::string& getHeterogeneousRateModel() const {
    return hetRateModel_;
  }

  void setUpCategoryTransitionMatrix();
     
private:
  void includeCategoryTransitionParameters_();
  
};

#endif
