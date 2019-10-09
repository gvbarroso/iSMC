/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 16/04/2018
 *
 */

#ifndef _HIDDENMARKOVMODELI_H_
#define _HIDDENMARKOVMODELI_H_

#include <vector>


class HiddenMarkovModelI { 
protected:
  bpp::Vdouble hiddenStatesInitializationProbabilities_;
  bpp::Vdouble scalingFactorsForward_;
  bpp::VVdouble forwardMatrix_;
  bpp::VVdouble backwardMatrix_;
  bpp::VVdouble posteriorProbMatrix_;
  std::vector< unsigned char > reconstructedHiddenStates_;
  
public:
  HiddenMarkovModelI():
  hiddenStatesInitializationProbabilities_(0),
  scalingFactorsForward_(),
  forwardMatrix_(0, bpp::Vdouble(0)),
  backwardMatrix_(0, bpp::Vdouble(0)),
  posteriorProbMatrix_(0, bpp::Vdouble(0)),
  reconstructedHiddenStates_(0)
  { }

  virtual ~HiddenMarkovModelI() {}

public:

  //NOTE many of these methods are redundant with more efficient HMM implementations (eg zipHMM) and are not really used by iSMC ATM.

  virtual double forwardAlgorithm() = 0;
    
  virtual void backwardAlgorithm() = 0;
 
  virtual void performPosteriorDecoding(size_t genomicStart, size_t genomicEnd) = 0;
  
  virtual void computePosteriorProbabilities() = 0;
       
  virtual const bpp::VVdouble& getForwardMatrix() const = 0;
  virtual bpp::VVdouble& getForwardMatrix() = 0;
  
  virtual const bpp::Vdouble& getScalingFactorsForward() const = 0;
  virtual bpp::Vdouble& getScalingFactorsForward() = 0;
  
  virtual const bpp::VVdouble& getBackwardMatrix() const = 0;
  virtual bpp::VVdouble& getBackwardMatrix() = 0;
    
  virtual const bpp::VVdouble& getPosteriorMatrix() const = 0;
  virtual bpp::VVdouble& getPosteriorMatrix() = 0;
  
  virtual const bpp::Vdouble& getInitializationProbabilities() const = 0;
  virtual bpp::Vdouble& getInitializationProbabilities() = 0;

private:
  virtual void setUpInitializationProbabilities_() = 0;
   
};

#endif
