/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 22/12/2018
 *
 */


#ifndef _HIDDENSTATESLIBRARY_H_
#define _HIDDENSTATESLIBRARY_H_

//#define DEBUG_HS

#include <string>
#include <cstring>
#include <algorithm> 
#include <cstdlib>
#include <utility>
#include <vector>

#include <boost/bimap.hpp>

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

typedef boost::bimap< size_t, std::string > ParameterAlphabet;

class HmmStatesLibrary {
private:
  size_t numberOfHiddenStates_;
  size_t numberOfModulatedParams_;
  std::vector< std::vector< unsigned char > > hiddenStatesCombinations_; //0 -> time; 1 -> theta; 2 -> rho; 3 -> ne  
  ParameterAlphabet parameterAlphabet_;
  //these maps are used for joint posterior decoding (index of observation == index in the vector):
  std::vector< std::pair< size_t, size_t > > treeTransitionMap_; //rho
  std::vector< std::pair< size_t, size_t > > treeToSnpMap_; //theta
  
public:
  HmmStatesLibrary(size_t numIntervals,
                   const std::vector<std::shared_ptr<bpp::DiscreteDistributionInterface>>& paramScalings,
                   const ParameterAlphabet& parameterAlphabet):
  numberOfHiddenStates_(0),
  numberOfModulatedParams_(0),
  hiddenStatesCombinations_(0, std::vector< unsigned char >(0)),
  parameterAlphabet_(parameterAlphabet),
  treeTransitionMap_(0),
  treeToSnpMap_(0)
  {
    computeNumberOfHiddenStates_(paramScalings, numIntervals);
    computeNumberOfModulatedParameters_(paramScalings);
    arrangeHiddenStatesCombinations_(paramScalings, numIntervals);
    #ifdef DEBUG_HS
    std::cout << "HS seq = G T R N" << std::endl;
    for(size_t i = 0; i < numberOfHiddenStates_; ++i) {
      std::cout << "HS [" << i << "] = ";  
      for(size_t j = 0; j < hiddenStatesCombinations_[i].size(); ++j) {
        std::cout << static_cast< unsigned int > (hiddenStatesCombinations_[i][j]) << " ";
      }
      std::cout << std::endl;
    }
    #endif
  }
  
public:
  size_t getNumberOfHiddenStates() {
    return numberOfHiddenStates_;
  }
  
  size_t getNumberOfModulatedParams() {
    return numberOfModulatedParams_;
  }

  std::vector< std::vector< unsigned char > >& getHiddenStates(){ 
    return hiddenStatesCombinations_;
  }
  const std::vector< std::vector< unsigned char > >& getHiddenStates() const { 
    return hiddenStatesCombinations_;
  }
  
  std::vector< unsigned char >& fetchHiddenState(size_t hiddenStateIndex) { 
    return hiddenStatesCombinations_[hiddenStateIndex];
  }
  const std::vector< unsigned char >& fetchHiddenState(size_t hiddenStateIndex) const { 
    return hiddenStatesCombinations_[hiddenStateIndex];
  } 
  
  const ParameterAlphabet& getParameterAlphabet() { 
    return parameterAlphabet_;
  }
  
  void initRhoEmissionsAlphabet(size_t numIntervals, bool missingData); 

  void initThetaEmissionsAlphabet(size_t numIntervals, size_t numObsStates); 
  
  std::vector< std::pair< size_t, size_t > > getTreeTransitionMap() {
    return treeTransitionMap_;
  }
  const std::vector< std::pair< size_t, size_t > > & getTreeTransitionMap() const {
    return treeTransitionMap_;
  }
  
  std::vector< std::pair< size_t, size_t > > getTreeToSnpMap() {
    return treeToSnpMap_;
  }
  const std::vector< std::pair< size_t, size_t > >& getTreeToSnpMap() const {
    return treeToSnpMap_;
  }
  
private:
  void computeNumberOfHiddenStates_(const std::vector<std::shared_ptr<bpp::DiscreteDistributionInterface>>& paramScalings, size_t numIntervals);

  void computeNumberOfModulatedParameters_(const std::vector<std::shared_ptr<bpp::DiscreteDistributionInterface>>& paramScalings);
  
  void arrangeHiddenStatesCombinations_(const std::vector< std::shared_ptr<bpp::DiscreteDistributionInterface>>& paramScalings, size_t numIntervals);
  
  std::vector< unsigned char > generateConversionSeries_(const std::vector<std::shared_ptr<bpp::DiscreteDistributionInterface>>& paramScalings);

  std::vector< unsigned char > buildStateComb_(const std::vector< unsigned char >& conversionSeries, unsigned int hsIndex);

};

#endif
