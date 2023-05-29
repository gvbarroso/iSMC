/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 08/05/2019
 *
 */

#ifndef _MATRIXOPTIMIZATIONFUNCTION_H_
#define _MATRIXOPTIMIZATIONFUNCTION_H_

#include <string>
#include <vector>

#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/ParameterList.h>

#include "MarkovModulatedSmc.h"
#include "ParametrizedMatrixI.h"


class MatrixOptimizationFunction:
  public bpp::FunctionInterface,
  public bpp::AbstractParameterAliasable {
protected:
  double totalDistance_;
  bpp::Vdouble matrixDistancesVector_;
  bpp::VVVdouble referenceMatrices_;
  std::shared_ptr< MarkovModulatedSmc > mmsmc_;
  std::shared_ptr< ParametrizedMatrixI > hmi_;
  
public:  
  MatrixOptimizationFunction(std::shared_ptr< MarkovModulatedSmc > mmsmc,
                             std::shared_ptr< ParametrizedMatrixI > hmi,
                             const bpp::ParameterList& params):
    AbstractParameterAliasable(""),
    totalDistance_(0.),
    matrixDistancesVector_(0),
    referenceMatrices_(0, bpp::VVdouble(0, bpp::Vdouble(0))),
    mmsmc_(mmsmc),
    hmi_(hmi)
  {
    addParameters_(params);
  }
 
  MatrixOptimizationFunction* clone() const { return new MatrixOptimizationFunction(*this); } 
  
public:
  void setParameters(const bpp::ParameterList& parameters) {
    AbstractParameterAliasable::setParametersValues(parameters);
  }

  double getValue() const { 
    return totalDistance_;
  }

  void fireParameterChanged(const bpp::ParameterList& parameters);
  
  //e.g. to import baum-welch-optimised ("proposed") matrices:
  void setReferenceMatrices(const bpp::VVVdouble& refMatrices) {
    referenceMatrices_ = refMatrices;
    matrixDistancesVector_.resize(referenceMatrices_.size());
  }
    
private:
  void computeAllMatrixDistances_();
    
  double computeSingleMatrixDistance_(const bpp::VVdouble& focalMatrix);
  
  double computeSumOfMatrixDistances_();
  
  double computeMeanMatrixDistance_() {
    computeAllMatrixDistances_();  
    double sumOfDistances = computeSumOfMatrixDistances_();
    return sumOfDistances / static_cast< double >(matrixDistancesVector_.size());  
  }
  
  double computeMedianMatrixDistance_();
  
  double computeMaximumMatrixDistance_() {
    computeAllMatrixDistances_();
    bpp::Vdouble::iterator it = std::max_element(matrixDistancesVector_.begin(), matrixDistancesVector_.end()); 
    return *it;
  }

  double computeMinimumMatrixDistance_() {
    computeAllMatrixDistances_();
    bpp::Vdouble::iterator it = std::min_element(matrixDistancesVector_.begin(), matrixDistancesVector_.end()); 
    return *it;
  }
  
};

#endif
