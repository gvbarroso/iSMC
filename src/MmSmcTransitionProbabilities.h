/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 08/01/2019
 *
 */


#ifndef _MMSMCTRANSITIONPROBABILITIES_H_
#define _MMSMCTRANSITIONPROBABILITIES_H_

#include <string>
#include <vector>
#include <utility>

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

#include "SmcTransitionProbabilities.h"
#include "MarkovModulatedSmc.h"


class MmSmcTransitionProbabilities:
  public SmcTransitionProbabilities {
protected:
  std::shared_ptr< MarkovModulatedSmc > mmsmc_;

public:
  MmSmcTransitionProbabilities(std::shared_ptr< MarkovModulatedSmc > mmsmc):
  SmcTransitionProbabilities(mmsmc),
  mmsmc_(mmsmc)
  {
    expectedMatrix_.resize(mmsmc -> getNumberOfHiddenStates(), bpp::Vdouble(mmsmc -> getNumberOfHiddenStates()));
    setUpExpectedMatrix();
  }

public:
  void setUpExpectedMatrix();

  bpp::VVVdouble fetchCompositeTransitionMatrix(bool missingData);
  
  bpp::VVVdouble fetchCompositeTransitionMatrix(bool missingData, size_t numRhoCateg, std::shared_ptr< bpp::DiscreteDistributionInterface > rhoScaling);

private:  
  void setTransitionProbabilitiesEqualTime_(const bpp::Vdouble& rowSumVector, const std::vector< std::vector < unsigned char > >& hiddenStates);
  
  //takes a particular lambdaVector and rho value and returns a transitionMatrix
  bpp::VVdouble computeTreeTransitionMatrix_(double rho, const bpp::ParameterList& lambdaVector);
  
  std::vector< bool > compareHiddenStates_(const std::vector< unsigned char >& hs1, const std::vector< unsigned char >& hs2);

};

#endif
