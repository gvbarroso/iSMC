/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 24/04/2019
 *
 */


#ifndef _SEQUENTIALLYMARKOVCOALESCENT_H_
#define _SEQUENTIALLYMARKOVCOALESCENT_H_

#include <string>
#include <vector>
#include <limits>

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/NumTools.h>

#include "OptionsContainer.h"


class SequentiallyMarkovCoalescent:
  public bpp::AbstractParameterAliasable {
protected:
  bpp::Vdouble timeIntervals_; 
  bpp::Vdouble averageCoalescenceTime_;
  bpp::ParameterList lambdaVector_; //stored as a ParameterList for speed

public:
  SequentiallyMarkovCoalescent(std::shared_ptr< OptionsContainer > smcOptions,
                               const std::vector< std::vector < unsigned char > >& sequences): 
  AbstractParameterAliasable(""),
  timeIntervals_(0),
  averageCoalescenceTime_(smcOptions -> getNumberOfIntervals()),
  lambdaVector_()
  { 
    addParameter_(new bpp::Parameter("theta", 1., bpp::Parameter::R_PLUS_STAR));
    computeMeanThetaAcrossDataSet_(sequences);
    addParameter_(new bpp::Parameter("rho", getParameterValue("theta") / 2., bpp::Parameter::R_PLUS_STAR));
    setTimeIntervals(smcOptions -> getNumberOfIntervals(),
                     smcOptions -> getTimeDisc(),
                     smcOptions -> getTmax());
  }
  SequentiallyMarkovCoalescent(unsigned int numIntervals,
                               const std::string& timeDisc,
                               double tMax,
                               const std::vector< std::vector < unsigned char > >& sequences): 
  AbstractParameterAliasable(""),
  timeIntervals_(numIntervals),
  averageCoalescenceTime_(numIntervals),
  lambdaVector_()
  { 
    addParameter_(new bpp::Parameter("theta", 1., bpp::Parameter::R_PLUS_STAR));
    computeMeanThetaAcrossDataSet_(sequences);
    addParameter_(new bpp::Parameter("rho", getParameterValue("theta") / 2., bpp::Parameter::R_PLUS_STAR));
    
    setTimeIntervals(numIntervals, timeDisc, tMax);
  }
  SequentiallyMarkovCoalescent(unsigned int numIntervals,
                               const std::string& timeDisc,
                               double tMax,
                               const bpp::ParameterList& rho_theta,
                               const bpp::ParameterList& lambdas): 
  AbstractParameterAliasable(""),
  timeIntervals_(numIntervals),
  averageCoalescenceTime_(numIntervals),
  lambdaVector_()
  { 
    addParameters_(rho_theta);
    setTimeIntervals(numIntervals, timeDisc, tMax);
    lambdaVector_.matchParametersValues(lambdas);
  }
  SequentiallyMarkovCoalescent(unsigned int numIntervals,
                               const std::string& timeDisc,
                               double tMax,
                               const bpp::ParameterList& optimParams): 
  AbstractParameterAliasable(""),
  timeIntervals_(numIntervals),
  averageCoalescenceTime_(numIntervals),
  lambdaVector_()
  { 
    addParameter_(new bpp::Parameter(optimParams.getParameter("theta")));
    addParameter_(new bpp::Parameter(optimParams.getParameter("rho")));
    
    setTimeIntervals(numIntervals, timeDisc, tMax);
  }

  SequentiallyMarkovCoalescent* clone() const { return new SequentiallyMarkovCoalescent(*this); } 

  virtual ~SequentiallyMarkovCoalescent() {}

public:
  bpp::Vdouble& getAverageCoalescenceTimeVector() { 
    return averageCoalescenceTime_;
  }
  const bpp::Vdouble& getAverageCoalescenceTimeVector() const { 
    return averageCoalescenceTime_;
  }
    
  bpp::Vdouble& getTimeIntervals() { 
    return timeIntervals_;
  } 
  const bpp::Vdouble& getTimeIntervals() const { 
    return timeIntervals_;
  } 
  
  virtual size_t getNumberOfHiddenStates() {
    return timeIntervals_.size();
  }
  
  size_t getNumberOfIntervals() { 
    return timeIntervals_.size();
  }
    
  const bpp::ParameterList& getLambdaVector() const {
    return lambdaVector_;
  }
  bpp::ParameterList& getLambdaVector() {
    return lambdaVector_;
  }
  
  void setTimeIntervals(unsigned int numIntervals, const std::string& discMethod, double tMax);

  void setTimeIntervals(unsigned int numIntervals, const std::string& discMethod, double tMax, const bpp::ParameterList& lambdas);
  
  double fetchLowerTimeBoundary(double time);
  
  double fetchUpperTimeBoundary(double time);
  
  unsigned int fetchTimeIndex(double time);
  
  void computeAverageCoalescenceTime();

  double getAverageCoalescenceTime(double time);
  
  bpp::ParameterList fetchIndependentLambdas();
  
  static bpp::ParameterList readParametersFromFile(const std::string& fileName);
    
private:
  void discretizeTimeUsingExpBoundaries_();
  
  void discretizeTimeUsingTmax_(double tMax);

  double fetchBiHaploidTheta_(const std::vector< unsigned char >& sequence);
  
  void computeMeanThetaAcrossDataSet_(const std::vector< std::vector < unsigned char > >& sequences);
  
};

#endif
