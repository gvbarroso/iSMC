/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 11/12/2022
 *
 */

#ifndef _PSMC_H_
#define _PSMC_H_

#include <vector>
#include <string>
#include <math.h>  
#include <utility>  

#ifdef SIMPLEZIPHMM
  #include <SimpleZipHMM/forwarder.hpp>
  #include <SimpleZipHMM/posterior_decoding.hpp>
  #include <SimpleZipHMM/hmm_io.hpp>
#else
  #include <zipHMM/forwarder.hpp>
  #include <zipHMM/posterior_decoding.hpp>
  #include <zipHMM/hmm_io.hpp>
#endif //SIMPLEZIPHMM

#include "SequentiallyMarkovCoalescent.h"
#include "MarkovModulatedSmc.h"
#include "SmcEmissionProbabilities.h"
#include "MmSmcEmissionProbabilities.h"
#include "SmcTransitionProbabilities.h"
#include "MmSmcTransitionProbabilities.h"
#include "HiddenMarkovModelI.h"
#include "PolymorphismData.h"

class Psmc:
  public HiddenMarkovModelI {
protected:
    
  std::string biHaploidName_;
  std::vector< unsigned char > biHaploidSnpCalling_;
  
  //for posterior decoding -> we slice biHaploidSnpCalling_ into selectedFragment_,
  //keeping biHaploidSnpCalling_ intact for iterating
  std::vector< unsigned char > selectedFragment_; 
  bpp::Vdouble posteriorAverageByPosition_;
  
  //breakpoints where the Markov chain will reset
  //(eg, to separate chromosomes or avoid long regions with missing data)
  //vector of pairs {start, end} of fragments that are INCLUDED in the analysis
  Breakpoints seqBreakpoints_; 
  
  double biHaploidLogLikelihood_;

  std::shared_ptr< SequentiallyMarkovCoalescent > smc_;
  std::shared_ptr< SmcTransitionProbabilities > smctp_;
  std::shared_ptr< SmcEmissionProbabilities > smcep_;

  std::vector<zipHMM::Forwarder> forwarders_;
  
public:
  Psmc(const std::string& biHaploidIndividualName,
       const std::vector< unsigned char >& biHaploidSnpCall,
       const Breakpoints& sequenceFragments,
       unsigned int missingBlockLength,
       std::shared_ptr< SequentiallyMarkovCoalescent > smc,
       std::shared_ptr< SmcTransitionProbabilities > smctp,
       std::shared_ptr< SmcEmissionProbabilities > smcep):
  HiddenMarkovModelI(),
  biHaploidName_(biHaploidIndividualName),
  biHaploidSnpCalling_(biHaploidSnpCall),
  selectedFragment_(0),
  posteriorAverageByPosition_(0),
  seqBreakpoints_(sequenceFragments),
  biHaploidLogLikelihood_(0.),
  smc_(smc),
  smctp_(smctp),
  smcep_(smcep),
  forwarders_(0)
  {
    setUpInitializationProbabilities_(); 
    std::cout << "Pre-processing pair of genomes " << biHaploidName_ << "..." << std::endl;
    if(missingBlockLength > 0) {
      seqBreakpoints_ = selectInformativeRegions_(missingBlockLength);
    }

    prepareDataStructures();
    std::cout << "done." << std::endl;
  }
  
  //for posterior decoding
  Psmc(const std::string& biHaploidIndividualName,
       const std::vector< unsigned char >& biHaploidSnpCall,
       const std::vector< std::pair< size_t, size_t > >& sequenceFragments,
       std::shared_ptr< SequentiallyMarkovCoalescent > smc,
       std::shared_ptr< SmcTransitionProbabilities > smctp,
       std::shared_ptr< SmcEmissionProbabilities > smcep):
  HiddenMarkovModelI(),
  biHaploidName_(biHaploidIndividualName),
  biHaploidSnpCalling_(biHaploidSnpCall),
  selectedFragment_(0),
  posteriorAverageByPosition_(0),
  seqBreakpoints_(sequenceFragments),
  biHaploidLogLikelihood_(0.),
  smc_(smc),
  smctp_(smctp),
  smcep_(smcep),
  forwarders_(0)
  {
    setUpInitializationProbabilities_();
  }
  
  virtual ~Psmc() {}
   
public:
  //not a const reference yet because version in daughter class is more elaborated
  virtual std::vector< unsigned char > getLocalTreesIndices() { 
    return reconstructedHiddenStates_;
  }
  
  const size_t getSequenceLength() const { 
    return biHaploidSnpCalling_.size();
  }
  
  const std::vector< unsigned char >& getBiHaploidSnpCalling() const { 
    return biHaploidSnpCalling_;
  }
  std::vector< unsigned char >& getBiHaploidSnpCalling() { 
    return biHaploidSnpCalling_;
  }
  
  const std::string& getBiHaploidName() const { 
    return biHaploidName_;
  }
  
  double getBiHaploidLogLikelihood() { 
    return biHaploidLogLikelihood_;
  }
  
  const Breakpoints& getSequenceBreakpoints() const {
    return seqBreakpoints_;
  }
  
  const std::vector< unsigned char >& getReconstructedStates() const { 
    return reconstructedHiddenStates_;
  }
  std::vector< unsigned char >& getReconstructedStates() { 
    return reconstructedHiddenStates_;
  }
    
  const bpp::Vdouble& getInitializationProbabilities() const { 
    return hiddenStatesInitializationProbabilities_;
  }
  bpp::Vdouble& getInitializationProbabilities() { 
    return hiddenStatesInitializationProbabilities_;
  }
    
  const bpp::VVdouble& getForwardMatrix() const { 
    return forwardMatrix_;
  }
  bpp::VVdouble& getForwardMatrix() { 
    return forwardMatrix_;
  }
  
  const bpp::Vdouble& getScalingFactorsForward() const { 
    return scalingFactorsForward_;
  }
  bpp::Vdouble& getScalingFactorsForward() { 
    return scalingFactorsForward_;
  }
  
  const bpp::VVdouble& getBackwardMatrix() const { 
    return backwardMatrix_;
  }
  bpp::VVdouble& getBackwardMatrix() { 
    return backwardMatrix_;
  }
    
  const bpp::VVdouble& getPosteriorMatrix() const { 
    return posteriorProbMatrix_;
  }
  bpp::VVdouble& getPosteriorMatrix() {
    return posteriorProbMatrix_;
  }
  
  void computeBiHaploidLogLikelihood(size_t numAvailThreads);
  
  void selectFragment(size_t genomicStart, size_t genomicEnd);
  
  std::vector< unsigned char > fetchFragment(size_t genomicStart, size_t genomicEnd);
  
  //writes data structures to file, based on seqBreakpoints_ 
  void writeDataStructures();

  void prepareDataStructures();
  
  double forwardAlgorithm();
    
  void backwardAlgorithm();
  
  void performPosteriorDecoding(size_t genomicStart, size_t genomicEnd);
  
  void posteriorDecodingUsingZipHMM(size_t genomicStart, size_t genomicEnd);
  
  void computePosteriorProbabilities();
  
protected:
  double computeBatchLogLikelihood_(const Breakpoints& bpBatch);

  //adds breakpoints to skip long regions with missing data, returns "extended" (updated) sequenceFragments
  Breakpoints selectInformativeRegions_(size_t thresholdLength);
    
  bpp::Vdouble fetchPerSiteShannonEquitability_(const bpp::VVdouble& postProbMatrix);
  
  bpp::Vdouble fetchPerSiteShannonEquitability_(const zipHMM::Matrix& postProbMatrix);

  void computePerSiteShannonEquitability_(const zipHMM::Matrix& postProbMatrix, const std::string& sequenceFileName);
  
  void computePerSiteShannonEquitability_(const bpp::VVdouble& postProbMatrix, const std::string& sequenceFileName);

  void setUpInitializationProbabilities_() {
      
    hiddenStatesInitializationProbabilities_.resize(smc_ -> getNumberOfHiddenStates());
    
    for(size_t i = 0; i < smc_ -> getNumberOfHiddenStates(); ++i) {
      hiddenStatesInitializationProbabilities_[i] = 1. / static_cast< double >(smc_ -> getNumberOfHiddenStates());
    }
  } 
  
};

#endif
