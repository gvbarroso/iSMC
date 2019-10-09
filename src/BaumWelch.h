/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 02/01/2019
 *
 */


#include <string>
#include <vector>
#include <math.h>
#include <thread>

#ifndef _BAUMWELCH_H_
#define _BAUMWELCH_H_

#include "Psmc.h"


class BaumWelch {
  
public:
  BaumWelch() {}
  
  virtual ~BaumWelch() {}
    
public:
  //for HMM layers
  void maximiseRateTransMat(const bpp::Vdouble& fwdScales, const bpp::VVdouble& fwdMat, const bpp::VVdouble& bckMat,
                            bpp::VVdouble& transMat, const bpp::VVdouble& emissMat,
                            const std::vector< std::vector< size_t > >& obsVec);

  //for HMM layers
  void maximiseRateEmissMat(const bpp::VVdouble& fwdMat, const bpp::VVdouble& bckMat, bpp::VVdouble& rateEmissMat,
                            const std::vector< std::vector< size_t > >& obsVec);

  //using data from a single diploid
  bpp::VVdouble proposeTransMat(const bpp::Vdouble& fwdScales, const bpp::VVdouble& fwdMat, const bpp::VVdouble& bckMat,
                                const bpp::VVdouble& transMat, const bpp::VVdouble& emissMat,
                                const std::vector< size_t >& seq);

  //using data from a single diploid
  bpp::VVdouble proposeEmissMat(const bpp::VVdouble& fwdMat, const bpp::VVdouble& bckMat, const bpp::VVdouble& emissMat, const std::vector< size_t >& seq);

  //to run BaumWelch in fragments
  void incrementTransitionMatrix(const zipHMM::Matrix& fwdMat, const zipHMM::Matrix& bckMat,
                                 bpp::VVdouble& newTransMat, const bpp::VVdouble& oldTransMat,
                                 const bpp::VVdouble& emissMat, const bpp::Vdouble& fwdScales,
                                 const std::vector< size_t >& frag);
  
};

#endif
