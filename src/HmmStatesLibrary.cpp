/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 22/12/2018
 *
 */


#include "HmmStatesLibrary.h"

using namespace bpp;
using namespace std;
using namespace boost::bimaps;


void HmmStatesLibrary::computeNumberOfHiddenStates_(const vector< shared_ptr< DiscreteDistributionInterface > >& paramScalings, size_t numIntervals) {
  size_t numberOfCategoriesCombinations = 1;
  for(size_t i = 0; i < paramScalings.size(); ++i) {
    if(paramScalings[i] -> hasParameter("heat")) { //if Gamma+Hotspot model
      numberOfCategoriesCombinations *= (paramScalings[i] -> getNumberOfCategories() + 1);  
    }
    else if(paramScalings[i] -> hasParameter("V2")) { //if Hotspot model
      numberOfCategoriesCombinations *= 2;    
    }
    else if(paramScalings[i] -> hasParameter("alpha")) { //if Gamma model
      numberOfCategoriesCombinations *= paramScalings[i] -> getNumberOfCategories();  
    }
  }
  numberOfHiddenStates_ = numberOfCategoriesCombinations * numIntervals;
  //cout << "computed no. HS = " << numberOfHiddenStates_ << endl;
}

//NOTE delete?
void HmmStatesLibrary::computeNumberOfModulatedParameters_(const vector< shared_ptr< DiscreteDistributionInterface > >& paramScalings) {
  for(size_t i = 0; i < paramScalings.size(); ++i) {
    if(paramScalings[i] -> getNumberOfCategories() > 1) { 
      ++numberOfModulatedParams_;
    }
    else if(paramScalings[i] -> hasParameter("heat")) { //if Gamma+Hotspot model with only 1 gamma category
      ++numberOfModulatedParams_;
    }
  }
} 

void HmmStatesLibrary::arrangeHiddenStatesCombinations_(const vector< shared_ptr< DiscreteDistributionInterface > >& paramScalings, size_t numIntervals) {
  //hidden states order: time, theta, rho and ne indices
  //regardless of how many of these are actually allowed to be heterogeneous
  hiddenStatesCombinations_.resize(numberOfHiddenStates_);
  for(size_t i = 0; i < numberOfHiddenStates_; ++i) {
    hiddenStatesCombinations_[i].resize(1 + paramScalings.size());
  }
  vector< unsigned char > conversionSeries = generateConversionSeries_(paramScalings);
  for(unsigned int i = 0; i < numberOfHiddenStates_; ++i) {
    hiddenStatesCombinations_[i] = buildStateComb_(conversionSeries, i); //stores the index of all parameters (stateComb)
  }
}
      
vector< unsigned char > HmmStatesLibrary::generateConversionSeries_(const vector< shared_ptr< DiscreteDistributionInterface > >& paramScalings) {
  size_t numParam = paramScalings.size() + 1; //number of spatial parameters + 1 (for time)
  vector< unsigned char > conversionSeries(numParam);
  size_t mult = 1;
  conversionSeries[numParam - 1] = static_cast< unsigned char >(mult);
  size_t numCat = 0; //number of categories of a given parameter
  for(size_t i = numParam - 1; i > 0; --i) {
    if(paramScalings[i - 1] -> hasParameter("heat")) { //if Gamma+Hotspots model
      numCat = paramScalings[i - 1] -> getNumberOfCategories() + 1;
    }
    else if(paramScalings[i - 1] -> hasParameter("V2")) { //if Hotspot model
      numCat = 2; 
    }
    else { //if Gamma model
      numCat = paramScalings[i - 1] -> getNumberOfCategories();  
    }
    mult *= numCat;
    conversionSeries[i - 1] = static_cast< unsigned char >(mult);
  }
  return conversionSeries;
}

vector< unsigned char > HmmStatesLibrary::buildStateComb_(const vector< unsigned char >& conversionSeries, unsigned int hsIndex) {
  size_t numParam = conversionSeries.size(); //in the normal case this is equal to 4 (time, theta, rho, Ne)
  vector< unsigned char > stateComb(numParam);
  //dealing with time
  unsigned int quo = hsIndex / conversionSeries[0];
  stateComb[0] = static_cast< unsigned char >(quo); 
  unsigned int rem = hsIndex % conversionSeries[0];
  quo = rem;
  //spatial parameters: 1 = theta, 2 = rho, 3 = neZero
  for(size_t i = 1; i < numParam; ++i) {
    stateComb[i] = static_cast< unsigned char >(quo / conversionSeries[i]);
    rem = quo % conversionSeries[i];
    quo = rem;
  }
  return stateComb;
}

void HmmStatesLibrary::initRhoEmissionsAlphabet(size_t numIntervals, bool missingData) {
  if(missingData) { //if there is missing data in the sequences, we add another "interval"
    ++numIntervals;
  }
  //cout << "HmmStatesLibrary::initRhoEmissionsAlphabet" << endl;
  treeTransitionMap_.resize(numIntervals * numIntervals);  
  size_t index = 0;
  for(size_t i = 0; i < numIntervals; ++i) {
    for(size_t j = 0; j < numIntervals; ++j) {
      treeTransitionMap_[index] = make_pair(i, j);
      //cout << "tts: {" << index << ", " << i << ", " << j << "}" << endl;
      ++index;
    }        
  }
}

void HmmStatesLibrary::initThetaEmissionsAlphabet(size_t numIntervals, size_t numObsStates) { 
  treeToSnpMap_.resize(numIntervals * numObsStates); //built-in mask
  size_t index = 0;
  for(size_t i = 0; i < numIntervals; ++i) {
    for(size_t j = 0; j < numObsStates; ++j) {
      treeToSnpMap_[index] = make_pair(i, j);
      //cout << "tes: {" << index << ", " << i << ", " << j << "}" << endl;
      ++index;
    }        
  }
}
