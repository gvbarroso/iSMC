/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 19/12/2018
 *
 */


#ifndef _MMSMCEMISSIONPROBABILITIES_H_
#define _MMSMCEMISSIONPROBABILITIES_H_

#include <vector>

#include "SmcEmissionProbabilities.h"
#include "MarkovModulatedSmc.h"


class MmSmcEmissionProbabilities:                 
  public SmcEmissionProbabilities {
protected:
  std::shared_ptr< MarkovModulatedSmc > mmsmc_;  
  
public:
  MmSmcEmissionProbabilities(std::shared_ptr< MarkovModulatedSmc > mmsmc,
                             const std::vector< std::vector< unsigned char > >& snpCalls):
  SmcEmissionProbabilities(mmsmc, snpCalls),
  mmsmc_(mmsmc)
  {
    expectedMatrix_.resize(mmsmc -> getNumberOfHiddenStates(), bpp::Vdouble(numberOfObservedStates_));
    setUpExpectedMatrix();
  }
  MmSmcEmissionProbabilities(std::shared_ptr< MarkovModulatedSmc > mmsmc,
                             unsigned int numObsStates):
  SmcEmissionProbabilities(mmsmc, numObsStates),
  mmsmc_(mmsmc)
  {
    expectedMatrix_.resize(mmsmc -> getNumberOfHiddenStates(), bpp::Vdouble(numberOfObservedStates_));
    setUpExpectedMatrix();
  }

public:
  void setUpExpectedMatrix();
  
  bpp::VVVdouble fetchCompositeEmissionMatrix();
      
  //HMM layers: emission probabilities of tree transitions given rho values
  bpp::VVdouble fetchRhoEmissions(const bpp::VVVdouble& compositeTransitionMatrix); 

  //HMM layers: emission probabilities of state (0, 1, 2) emissions given theta values
  bpp::VVdouble fetchThetaEmissions(const bpp::VVVdouble& compositeEmissionMatrix);

};

#endif