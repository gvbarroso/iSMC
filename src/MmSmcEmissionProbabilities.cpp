/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 22/12/2018
 *
 */


#include "MmSmcEmissionProbabilities.h"

using namespace std;
using namespace bpp;

  
void MmSmcEmissionProbabilities::setUpExpectedMatrix() {  
    
  for(size_t i = 0; i < mmsmc_ -> getNumberOfHiddenStates(); ++i) {
      
    unsigned int timeI = mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[i][0]; //time index
    unsigned int thetaI = mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[i][1]; //theta index
    
    double theta = mmsmc_ -> getParameterValue("theta");
    
    //scales theta if theta-modulated model
    if(mmsmc_ -> getParameterTransitions()[0] -> getNumberOfCategories() > 1) {
        
      //Due to a discrepancy between # categories inside Param Scalings and Param Transitions,
      //the Gamma+Hotspot model needs to be handled differently
      if(mmsmc_ -> getParameterTransitions()[0] -> getHeterogeneousRateModel() == "Gamma+Hotspot") {
          
        if(thetaI == 0) { //hotspot  
          theta *= mmsmc_ -> getParameterScalings()[0] -> getParameterValue("heat");
        }
        
        else { //gamma categories (-1 to acess the correct indices inside the gamma dist.)
          theta *= mmsmc_ -> getParameterScalings()[0] -> getCategories()[thetaI - 1];    
        }
      }
      
      else {
        theta *= mmsmc_ -> getParameterScalings()[0] -> getCategories()[thetaI];
      }
      
    }
    
    //time interval (~tree) and its upper boundary:
    double lowerTimeBoundary = mmsmc_ -> getTimeIntervals()[timeI];
    double upperTimeBoundary = mmsmc_ -> fetchUpperTimeBoundary(lowerTimeBoundary);
    
    if(upperTimeBoundary == numeric_limits< double >::max()) { 
      //if the upper time boundary is infinity we use 2x the avg. coal. time instead
      upperTimeBoundary = 2. * mmsmc_ -> getAverageCoalescenceTime(lowerTimeBoundary);
    }
    
    //for every type of observation (0 = homozygote; 1 = heterozygote; 2 = missing data):
    for(size_t k = 0; k < numberOfObservedStates_; ++k) {
      expectedMatrix_[i][k] = computeEmissionProbabilityUsingIntegration_(lowerTimeBoundary, upperTimeBoundary, k, theta);
    }
    
  }
}  

//composition of transition matrices for different values of theta (but not scaled by theta trans probs as in MMHMM)
VVVdouble MmSmcEmissionProbabilities::fetchCompositeEmissionMatrix() {
    
  size_t numThetaCateg = mmsmc_ -> getParameterTransitions()[0] -> getNumberOfCategories();
  size_t numIntervals = mmsmc_ -> getNumberOfIntervals();
  
  //composition of transition matrices
  VVVdouble cem(numThetaCateg, VVdouble(numIntervals, Vdouble(numberOfObservedStates_)));
  
  for(size_t x = 0; x < numThetaCateg; ++x) {
      
    double theta = mmsmc_ -> getParameterValue("theta") * mmsmc_ -> getParameterScalings()[0] -> getCategories()[x];
    
    //gets (x + 1)th SMC emission matrix:
    for(size_t i = 0; i < numIntervals; ++i) {
        
      //time interval (~tree) and its upper boundary:
      double lowerTimeBoundary = mmsmc_ -> getTimeIntervals()[i];
      double upperTimeBoundary = mmsmc_ -> fetchUpperTimeBoundary(lowerTimeBoundary);
      
      if(upperTimeBoundary == numeric_limits< double >::max()) { 
        //if the upper time boundary is infinity we use 2x the avg. coal. time instead
        upperTimeBoundary = 2. * mmsmc_ -> getAverageCoalescenceTime(lowerTimeBoundary);
      }
      
      //for every type of observation (0 = homozygote; 1 = heterozygote; 2 = missing data):
      for(size_t k = 0; k < numberOfObservedStates_; ++k) {
        cem[x][i][k] = computeEmissionProbabilityUsingIntegration_(lowerTimeBoundary, upperTimeBoundary, k, theta);
      }
      
    }
  }
  
  return cem;
}

VVdouble MmSmcEmissionProbabilities::fetchRhoEmissions(const VVVdouble& compositeTransitionMatrix) { 
    
  size_t numStates = compositeTransitionMatrix.front().size(); //has info on the existence of missing data
  size_t numRhoCateg = mmsmc_ -> getParameterTransitions()[1] -> getNumberOfCategories();
  
  VVdouble rhoEmissionMatrix(numRhoCateg, Vdouble(numStates * numStates));
  
  for(size_t i = 0; i < numRhoCateg; ++i) {
      
    for(size_t j = 0; j < mmsmc_ -> getHmmStatesLibrary() -> getTreeTransitionMap().size(); ++j) {
        
      pair< size_t, size_t > treePair = mmsmc_ -> getHmmStatesLibrary() -> getTreeTransitionMap()[j];
      
      //gets SMC transitions
      size_t from = get<0>(treePair);
      size_t to = get<1>(treePair);
      
      rhoEmissionMatrix[i][j] = compositeTransitionMatrix[i][from][to]; //
      //cout << "cat. " << i << " rho emiss. from " << from << " to " << to << " = " << rhoEmissionMatrix[i][j] << endl;
    }
  }
  
  return rhoEmissionMatrix;
}

VVdouble MmSmcEmissionProbabilities::fetchThetaEmissions(const VVVdouble& compositeEmissionMatrix) { 
    
  size_t numIntervals = mmsmc_ -> getNumberOfIntervals(); 
  size_t numThetaCateg = mmsmc_ -> getParameterTransitions()[0] -> getNumberOfCategories();
  
  VVdouble thetaEmissionMatrix(numThetaCateg, Vdouble(numIntervals * numberOfObservedStates_));
  
  for(size_t i = 0; i < numThetaCateg; ++i) {
      
    for(size_t j = 0; j < mmsmc_ -> getHmmStatesLibrary() -> getTreeToSnpMap().size(); ++j) {
        
      pair< size_t, size_t > treeToSnp =  mmsmc_ -> getHmmStatesLibrary() -> getTreeToSnpMap()[j];
      
      //gets SMC emissions
      size_t tree = get<0>(treeToSnp);
      size_t state = get<1>(treeToSnp);
      
      thetaEmissionMatrix[i][j] = compositeEmissionMatrix[i][tree][state]; //NOTE this probably has built-in masking
      //cout << "index: " << j << "; tree: " << tree << "; state: " << state << "; prob = " << thetaEmissionMatrix[i][j] << endl;
    }
  }
  
  return thetaEmissionMatrix;
}