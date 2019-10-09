/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 08/05/2019
 *
 */


#ifndef _MULTIPLEMMPSMC_H_
#define _MULTIPLEMMPSMC_H_

#include <string>
#include <vector>

#include "MmPsmc.h"


class MultipleMmPsmc { //organizes the (independent) MmPsmc objects
private:
  std::vector< std::shared_ptr< MmPsmc > > psmcVec_; //the entire dataset
  size_t testIndex_; //used for Leave-One-Out cross-validation
  
public:
  MultipleMmPsmc(const std::vector< std::shared_ptr< MmPsmc > >& psmcVector):
  psmcVec_(psmcVector),
  testIndex_(psmcVector.size()) //init to out-of-bounds index
  { }
  
public:
  const std::vector< std::shared_ptr< MmPsmc > >& getWholePsmcVector() const {                              
    return psmcVec_;                                                                                   
  }

  std::vector< std::shared_ptr< MmPsmc > > fetchTrainingPsmcVector();
  
  size_t getTestIndex() {
    return testIndex_;
  }
  
  void setTestIndex(size_t newIndex) {
    testIndex_ = newIndex;
  }
  

  double getTestLogLikelihood() { 
    return psmcVec_[testIndex_] -> getBiHaploidLogLikelihood();
  }

  void jointForward(const std::vector< std::vector< size_t > >& obsVector, const bpp::VVdouble& rateTrans, const bpp::Vdouble& pi,
                    const bpp::VVdouble& rateEmiss, bpp::Vdouble& forwardScales, bpp::VVdouble& rateForwardMatrix);
  
  void jointBackward(const std::vector< std::vector< size_t > >& obsVector, const bpp::VVdouble& rateTrans,
                     const bpp::VVdouble& rateEmiss, const bpp::Vdouble& forwardScales, bpp::VVdouble& rateBackwardMatrix);

  void updatePsmcBackwardMatrices();
 
  double fetchCompositeLogLikelihood(); 
  
  double fetchMedianLogLikelihood(); 
  
  double fetchMinimumLogLikelihood();
  
private:
  void computeBatchLogLikelihood_(const std::vector< std::shared_ptr< MmPsmc > >& focalPsmcVec, size_t numAvailThreads);

  void computeIndividualsLogLikelihood_();

};

#endif
