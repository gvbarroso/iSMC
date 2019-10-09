/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 09/09/2019
 *
 */


#ifndef _MMPSMC_H_
#define _MMPSMC_H_

#include <iostream>
#include <vector>
#include <string>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/bimap.hpp>

#include <Bpp/Numeric/ParameterList.h>

#ifdef SIMPLEZIPHMM
  #include <SimpleZipHMM/forwarder.hpp>
  #include <SimpleZipHMM/viterbi.hpp>
  #include <SimpleZipHMM/posterior_decoding.hpp>
  #include <SimpleZipHMM/hmm_io.hpp>
#else
  #include <zipHMM/viterbi.hpp>
  #include <zipHMM/forwarder.hpp>
  #include <zipHMM/posterior_decoding.hpp>
  #include <zipHMM/hmm_io.hpp>
#endif //SIMPLEZIPHMM

#include "MarkovModulatedSmc.h"
#include "MmSmcEmissionProbabilities.h"
#include "MmSmcTransitionProbabilities.h"
#include "Psmc.h"


class MmPsmc: 
  public Psmc {
protected:
  std::shared_ptr< MarkovModulatedSmc > mmsmc_;
  std::shared_ptr< MmSmcTransitionProbabilities > mmtp_;
  std::shared_ptr< MmSmcEmissionProbabilities > mmep_;

public:
  MmPsmc(const std::string& biHaploidIndividualName,
         const std::vector< unsigned char >& biHaploidSnpCall,
         const std::vector< std::pair< size_t, size_t > >& sequenceFragments,
         unsigned int missingBlockLength,
         std::shared_ptr< MarkovModulatedSmc > mmsmc,
         std::shared_ptr< MmSmcTransitionProbabilities > mmtp,
         std::shared_ptr< MmSmcEmissionProbabilities > mmep):
  Psmc(biHaploidIndividualName, biHaploidSnpCall, sequenceFragments, missingBlockLength, mmsmc, mmtp, mmep),
  mmsmc_(mmsmc),
  mmtp_(mmtp),
  mmep_(mmep)
  { }
  
  //used only for posterior decoding
  MmPsmc(const std::string& biHaploidIndividualName,
         const std::vector< unsigned char >& biHaploidSnpCall,
         const std::vector< std::pair< size_t, size_t > >& sequenceFragments,
         std::shared_ptr< MarkovModulatedSmc > mmsmc,
         std::shared_ptr< MmSmcTransitionProbabilities > mmtp,
         std::shared_ptr< MmSmcEmissionProbabilities > mmep):
  Psmc(biHaploidIndividualName, biHaploidSnpCall, sequenceFragments, mmsmc, mmtp, mmep),
  mmsmc_(mmsmc),
  mmtp_(mmtp),
  mmep_(mmep)
  { }
  
public: 
  std::vector< size_t > fetchTransitionSequence(std::vector< size_t >& treeSequence, size_t start, size_t end, bool missingData);
  
  std::vector< size_t > fetchEmissionSequence(const std::vector< size_t >& treeSequence, size_t start, size_t end);

  std::vector< unsigned char > fetchLocalTreesIndices(); 
  
  std::vector< size_t > fetchLocalTrees(const zipHMM::Matrix& postProbMatrix);
  
  std::vector< size_t > fetchLocalTrees(const zipHMM::Matrix& postProbMatrix, bpp::Vdouble& treeProbs);
  
  void decodeLocalTmrca(const zipHMM::Matrix& postProbMatrix, const std::string& fileName);

  std::vector< size_t > fetchLocalRateIndices(const zipHMM::Matrix& postProbMatrix, const std::string& rate);

  std::vector< unsigned char > fetchLocalRateIndices(const std::string& rate);

  void computeLocalAverageRate(const std::string& rate);
  
  void computeLocalAverageRate(const zipHMM::Matrix& postProbMatrix, const std::string& fileName);

  void decodeLocalRate(const bpp::VVdouble& ratePostProbs, const std::string& fileName, const std::string& rate);
  
  void decodeLocalRate(const zipHMM::Matrix& ratePostProbs, const std::string& fileName, const std::string& rate);

  void computeLocalAverageTmrca();

  //memory-efficient version, computes and directly writes to file
  void computeLocalAverageTmrca(const zipHMM::Matrix& postProbMatrix, const std::string& fileName);
  
  //memory-efficient version, computes and directly writes to file
  void posteriorDecodingUsingZipHMM(size_t genomicStart, size_t genomicEnd, const std::string& fileName, bool restricted, zipHMM::Matrix& postProbMatrix, bpp::Vdouble& pi);
  
  void computePosteriorProbs(size_t genomicStart, size_t genomicEnd, zipHMM::Matrix& postProbMatrix, bpp::Vdouble& pi);

  //gets posterior probability matrix considering only rate categories as hidden states
  bpp::VVdouble getRatePosteriorMatrix(const bpp::VVdouble& postProbMatrix, const std::string& rate);
  
  //gets posterior probability matrix considering only rate categories as hidden states
  bpp::VVdouble getRatePosteriorMatrix(const zipHMM::Matrix& postProbMatrix, const std::string& rate);
  
  void printTmrcaPosteriorMatrix(const zipHMM::Matrix& postProbMatrix, const std::string& fileName);
  
  void printRatePerInterval(const zipHMM::Matrix& postProbMatrix, const std::string& rate, const std::string& fileName);
  
  void printDecodedStates(const zipHMM::Matrix& postProbMatrix, const std::string& fileName);
  
private:
  void computeGammaRateMaxProb_(const bpp::VVdouble& ratePostProbs, const std::string& rate, boost::iostreams::filtering_ostream& file);
  
  void computeGammaRateMaxProb_(const zipHMM::Matrix& ratePostProbs, const std::string& rate, boost::iostreams::filtering_ostream& file);

  void computeHotspotRateMaxProb_(const bpp::VVdouble& ratePostProbs, const std::string& rate, boost::iostreams::filtering_ostream& file);
  
  void computeHotspotRateMaxProb_(const zipHMM::Matrix& postProbMatrix, const std::string& rate, boost::iostreams::filtering_ostream& file);

  void computeGammaWithHotspotRateMaxProb_(const bpp::VVdouble& ratePostProbs, const std::string& rate, boost::iostreams::filtering_ostream& file);
  
  void computeGammaWithHotspotRateMaxProb_(const zipHMM::Matrix& ratePostProbs, const std::string& rate, boost::iostreams::filtering_ostream& file);

  void computeHotspotAverageRate_(const std::string& rate);

  void computeHotspotAverageRate_(const zipHMM::Matrix& postProbMatrix, const std::string& rate, boost::iostreams::filtering_ostream& file);
  
  void computeMonoModulatedHotspotAverageRate_(const zipHMM::Matrix& postProbMatrix, const std::string& rate, boost::iostreams::filtering_ostream& file);

  void computeGammaAverageRate_(const std::string& rate);

  void computeGammaAverageRate_(const zipHMM::Matrix& postProbMatrix, const std::string& rate, boost::iostreams::filtering_ostream& file);
  
  bpp::VVdouble computeAverageRatePerInterval_(const zipHMM::Matrix& postProbMatrix, const std::string& rate);
  
  bpp::Vdouble computeAverageRateWithinDecodedTrees_(const zipHMM::Matrix& postProbMatrix, const std::vector< size_t >& decodedTmrca,
                                                     const bpp::Vdouble& maxTreeProbs, const std::string& rate);
  
  void computeGammaWithHotspotAverageRate_(const std::string& rate);

  void computeGammaWithHotspotAverageRate_(const zipHMM::Matrix& postProbMatrix, const std::string& rate, boost::iostreams::filtering_ostream& file);
  
  void computeMonoModulatedGammaWithHotspotAverageRate_(const zipHMM::Matrix& postProbMatrix, const std::string& rate, boost::iostreams::filtering_ostream& file);

  bpp::VVdouble getTreePosteriorMatrix_(const bpp::VVdouble& postProbMatrix);

  bpp::VVdouble getTreePosteriorMatrix_(const zipHMM::Matrix& postProbMatrix);
  
  bpp::VVdouble getRatePosteriorMatrix_(const bpp::VVdouble& postProbMatrix, const std::string& rate);

  bpp::VVdouble getRatePosteriorMatrix_(const zipHMM::Matrix& postProbMatrix, const std::string& rate);

  void maskTreeSequence_(std::vector< size_t >& treeSeq, size_t start, size_t end);

};

#endif
