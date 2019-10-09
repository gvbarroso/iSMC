/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 09/07/2018
 * Last modified: 03/09/2019
 *
 */


#ifndef _SMCDECODINGWRAPPER_
#define _SMCDECODINGWRAPPER_

#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>  
#include <thread>
#include <utility>  

#include <boost/bimap.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "MarkovModulatedSmc.h"
#include "MmSmcEmissionProbabilities.h"
#include "MmSmcTransitionProbabilities.h"
#include "MultipleMmPsmc.h"
#include "OptionsContainer.h"
#include "BaumWelch.h"


class SmcDecodingWrapper {
 
private:
  std::shared_ptr< MarkovModulatedSmc > mmsmc_;
  std::shared_ptr< MultipleMmPsmc > mPsmc_;
  std::shared_ptr< OptionsContainer > smcOptions_;
  std::shared_ptr< MmSmcTransitionProbabilities > tp_;
  std::shared_ptr< MmSmcEmissionProbabilities > ep_;
  
public:
  SmcDecodingWrapper(std::shared_ptr< MarkovModulatedSmc > mmsmc,
                     std::shared_ptr< MultipleMmPsmc > mPsmc,
                     std::shared_ptr< OptionsContainer > smcOptions,
                     std::shared_ptr< MmSmcTransitionProbabilities > tp,
                     std::shared_ptr< MmSmcEmissionProbabilities > ep):
  mmsmc_(mmsmc),
  mPsmc_(mPsmc),
  smcOptions_(smcOptions),
  tp_(tp),
  ep_(ep)
  { }
  
public:  
  void writePosteriorCoordinatesToFile();
  
  void writeDiploidLabelsToFile();
  
  void decodeDiploid();
  
  void computeAverageLandscapes();
      
  bpp::Vdouble readLandscapeFromFile(const std::string& fileName);

  void writeSingleNucleotideLandscapeToFile(const bpp::Vdouble& landscape, const std::string& name);
      
  void decodeDataset(bool missing);
  
  void jointlyDecodeUsingBaumWelch(bpp::VVdouble& rateEmiss, bpp::VVdouble& rateTrans, const std::string& rate, bool missing);
  
private:
  std::vector< std::vector< size_t > > fetchLocalTmrcaDist_(bpp::VVdouble& piMmSmc, size_t fragStart, size_t fragEnd, size_t fragId);  

  std::vector< std::vector < size_t > > fetchWholeTmrcaDist_();

  void jointlyDecodeFragment_(bpp::VVdouble& pi, size_t fragStart, size_t fragEnd, size_t blockId, size_t fragId, size_t numAvailThreads, bool missing);
  
  void decodeBatchOfBreakpoints_(const std::vector< std::pair< size_t, size_t > >& batchBps, const std::vector< size_t >& indexVector);
  
  void computeAverageRateOverDiploids_(const std::string& rate);
  
  void decodeBatchOfBreakpoints_(size_t numAvailThreads, bool missing,
                                 const std::vector< std::pair< size_t, size_t > >& batchBreakpoints,
                                 const std::vector< size_t >& blockIdVector);
};

#endif
