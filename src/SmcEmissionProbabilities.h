/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 02/01/2019
 *
 */


#ifndef _SMCEMISSIONPROBABILITIES_H_
#define _SMCEMISSIONPROBABILITIES_H_

#include <string>
#include <vector>
#include <iterator>
#include <limits>
#include <math.h>  
#include <set>

#include "SequentiallyMarkovCoalescent.h"
#include "AbstractParametrizedMatrix.h"


class SmcEmissionProbabilities:                 
  public AbstractParametrizedMatrix {

protected:
  std::shared_ptr< SequentiallyMarkovCoalescent > smc_;
  size_t numberOfObservedStates_;
  
public:
  SmcEmissionProbabilities(std::shared_ptr< SequentiallyMarkovCoalescent > smc,
                           const std::vector< std::vector< unsigned char > >& snpCalls):
  AbstractParametrizedMatrix(),
  smc_(smc),
  numberOfObservedStates_(0)  
  {
    computeNumberOfHmmObservedStates_(snpCalls);  
    expectedMatrix_.resize(smc -> getNumberOfIntervals(), bpp::Vdouble(numberOfObservedStates_));
    setUpExpectedMatrix();
  }
  SmcEmissionProbabilities(std::shared_ptr< SequentiallyMarkovCoalescent > smc,
                           size_t numObsStates):
  AbstractParametrizedMatrix(),
  smc_(smc),
  numberOfObservedStates_(numObsStates)  
  {
    expectedMatrix_.resize(smc -> getNumberOfIntervals(), bpp::Vdouble(numObsStates));
    setUpExpectedMatrix();
  }

public:
  size_t getNumberOfObservedStates() { 
    return numberOfObservedStates_;
  }
  
  void setUpExpectedMatrix();
  
  //to get the number of observed states in a particular fragment (e.g. when writing data structures)
  size_t fetchNumberOfObservedStates(const std::vector < unsigned char >& seq) {
    return std::set< unsigned char >(begin(seq), end(seq)).size();
  }

protected:
  double computeEmissionProbabilityUsingIntegration_(double lowerTime, double upperTime, size_t obs, double theta);
  
  //to get the number of observed states in the whole sequence, in all diploids
  void computeNumberOfHmmObservedStates_(const std::vector< std::vector < unsigned char > >& allSnpCalls);
  
};

#endif