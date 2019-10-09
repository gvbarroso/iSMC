/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 02/01/2019
 *
 */

#include "SmcEmissionProbabilities.h"

using namespace bpp;
using namespace std;


void SmcEmissionProbabilities::setUpExpectedMatrix() { 
    
  for(size_t i = 0; i < smc_ -> getNumberOfIntervals(); ++i) {
      
    double lowerTimeBoundary = smc_ -> getTimeIntervals()[i];
    double upperTimeBoundary = smc_ -> fetchUpperTimeBoundary(lowerTimeBoundary);
    
    if(upperTimeBoundary == std::numeric_limits< double >::max()) { 
      //if the upper time boundary is infinity we use 2x the avg. coal. time instead
      upperTimeBoundary = 2. * smc_ -> getAverageCoalescenceTime(smc_ -> getTimeIntervals()[i]);
    }
    
    //0 = homozygote; 1 = heterozygote; 2 = missing data:
    for(size_t k = 0; k < numberOfObservedStates_; ++k) {
      expectedMatrix_[i][k] = computeEmissionProbabilityUsingIntegration_(lowerTimeBoundary, upperTimeBoundary, k, smc_ -> getParameterValue("theta"));
    }
  }
}  

double SmcEmissionProbabilities::computeEmissionProbabilityUsingIntegration_(double lowerTime, double upperTime, size_t obs, double theta) {
  if(obs == 0) { //if homozygote:
    return (exp(- lowerTime *  theta) / theta - exp(- upperTime * theta) / theta) / (upperTime - lowerTime);
  } 
  else if(obs == 1) { //if heterozygote:
    return 1. - (exp(- lowerTime *  theta) / theta - exp( - upperTime * theta) / theta) / (upperTime - lowerTime);
  }
  else if(obs == 2) { //if missing data:
    return 1.;
  } 
  else {
    throw Exception("iSMC::Attempted to compute emission prob. for mis-specified obs!");
  } 
}
  
void SmcEmissionProbabilities::computeNumberOfHmmObservedStates_(const vector< vector < unsigned char > >& allSnpCalls) { 

  size_t numDiploids = allSnpCalls.size();
  vector< size_t > allUniqueCounts(numDiploids);

  for(size_t i = 0; i < numDiploids; ++i) { //for all diploids
    allUniqueCounts[i] = fetchNumberOfObservedStates(allSnpCalls[i]);
  }

  if(!equal(begin(allUniqueCounts) + 1, end(allUniqueCounts), begin(allUniqueCounts))) {
    throw Exception("iSMC::Not all diploids have same number of observed states! Check *_seqs.txt file");
  }

  numberOfObservedStates_ = allUniqueCounts.front();
}


