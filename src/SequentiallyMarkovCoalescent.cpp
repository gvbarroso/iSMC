/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 13/06/2019
 *
 */


#include <boost/bimap.hpp>
#include <boost/algorithm/string.hpp>

#include "OptionsContainer.h"
#include "SequentiallyMarkovCoalescent.h"

using namespace bpp;
using namespace std;


void SequentiallyMarkovCoalescent::setTimeIntervals(unsigned int numIntervals, const string& discMethod, double tMax) {
    
  averageCoalescenceTime_.resize(numIntervals);  
  
  timeIntervals_.resize(numIntervals);
  
  for(size_t i = 0; i < timeIntervals_.size(); ++i) {
      
    if(!lambdaVector_.hasParameter("l" + bpp::TextTools::toString(i))){
      lambdaVector_.addParameter(new bpp::Parameter("l" + bpp::TextTools::toString(i), 1., &bpp::Parameter::R_PLUS_STAR));
    }
  }
  
  if(discMethod == "quantiles") { //quantiles of exp. distribution (Schiffels & Durbin 2014)
    discretizeTimeUsingExpBoundaries_();
  }
  
  else if(discMethod == "log_even") { //even in log-space (Li & Durbin 2011)
    discretizeTimeUsingTmax_(tMax); //new default (changed 28/12/18)  
  }
  
  else {
    throw Exception("iSMC::mis-specified method for time discretisation: " + discMethod);
  }
  
  computeAverageCoalescenceTime();
}

void SequentiallyMarkovCoalescent::setTimeIntervals(unsigned int numIntervals, const string& discMethod, double tMax, const bpp::ParameterList& lambdas) {
    
  averageCoalescenceTime_.resize(numIntervals);  
  timeIntervals_.resize(numIntervals);
  lambdaVector_ = lambdas;
  
  if(discMethod == "quantiles") { //quantiles of exp. distribution (Schiffels & Durbin 2014)
    discretizeTimeUsingExpBoundaries_();
  }
  
  else if(discMethod == "log_even"){ //even in log-space (Li & Durbin 2011)
    discretizeTimeUsingTmax_(tMax); //new default (changed 28/12/18)  
  }
  
  else {
    throw Exception("iSMC::mis-specified method for time discretisation: " + discMethod);
  }
  
  computeAverageCoalescenceTime();
}
  
double SequentiallyMarkovCoalescent::fetchLowerTimeBoundary(double time) {
    
  double lowerBoundary = timeIntervals_[timeIntervals_.size() - 1];
  
  for(size_t i = 0; i < timeIntervals_.size(); ++i) {
      
    if(time < timeIntervals_[i]) {
      lowerBoundary = timeIntervals_[i - 1];
      
      break;
    }
    
    else if (time == timeIntervals_[i]) {
      lowerBoundary = timeIntervals_[i];
      
      break;
    }
  }
  
  return lowerBoundary;
}
  
double SequentiallyMarkovCoalescent::fetchUpperTimeBoundary(double time) {  
    
  double upperBoundary = numeric_limits< double >::max();
  
  for(size_t i = 0; i < timeIntervals_.size(); ++i) {
      
    if(time < timeIntervals_[i]) {
      upperBoundary = timeIntervals_[i];
      
      break; 
    }
    
    else if(time == timeIntervals_[i]) { //else if time is exactly a time boundary,
        
      //but we aren't talking about the LAST time boundary, upperBound is simple the next:
      if(time != timeIntervals_[timeIntervals_.size() - 1]) { 
        upperBoundary = timeIntervals_[i + 1];
      }
      
      break; //if time is exactly the last time boundary, its upper boundary remains infinity
    }
  }
  
  return upperBoundary;
}
    
unsigned int SequentiallyMarkovCoalescent::fetchTimeIndex(double time) {
    
  for(size_t i = 0; i < timeIntervals_.size(); ++i) {
      
    if(time == timeIntervals_[i]) {
      return static_cast< unsigned int >(i);
    }
  }
  
  throw Exception("iSMC::Did not find a match for time " + TextTools::toString(time));
}
  
void SequentiallyMarkovCoalescent::computeAverageCoalescenceTime() {
  
  for(size_t i = 0; i < timeIntervals_.size() - 1; ++i) {
      
    double deltaTime  = timeIntervals_[i + 1] - timeIntervals_[i];
    double lambda = lambdaVector_[i].getValue();
    
    if(lambda < 0.001) { //Dealing with numerical instabilities (see Schiffels & Durbin 2014 Supp. Note)
      averageCoalescenceTime_[i] = (timeIntervals_[i] + timeIntervals_[i + 1]) / 2.;
    }
    else {
      //apply the analytical approximation
      averageCoalescenceTime_[i] = (1. / ((1. - exp(-deltaTime * lambda)) * lambda)) *
                                   (1. + lambda * timeIntervals_[i] - exp(-lambda * deltaTime) *
                                   (1. + lambda * timeIntervals_[i + 1]));
    }
  }
  
  //At the very last time index, the upperBound is "infinity", and so is deltaTime -> approximate 
  size_t index = timeIntervals_.size() - 1;
  double deltaTime = numeric_limits< double >::max();
  double lambda = lambdaVector_[index].getValue();
  
  //Dealing with numerical instabilities (see Schiffels & Durbin 2014 Supp. Note)
  if(lambda < 0.001) { 
    averageCoalescenceTime_[index] = timeIntervals_[index] + (1. / lambda);
  }
  
  else {
    //as in MSMC:                                 
    averageCoalescenceTime_[index] = (1. / ((1. - exp(-deltaTime * lambda)) * lambda)) *
                                     (1. + lambda * timeIntervals_[index]);
  }
  //cout << "last avg. coal. time: " << averageCoalescenceTime_[index] << endl;

}

double SequentiallyMarkovCoalescent::getAverageCoalescenceTime(double time) {
    
  for(size_t i = 0; i < timeIntervals_.size(); ++i) {
      
    if(time == timeIntervals_[i]) {
      return averageCoalescenceTime_[i];
    }
  }
  
  throw Exception("iSMC::Did not find avg coal time for time = " + TextTools::toString(time));
}
    
ParameterList SequentiallyMarkovCoalescent::readParametersFromFile(const string& fileName) { 
  
  ParameterList inputParameters;
  ifstream optimParamsFile;
  optimParamsFile.open(fileName, ios::in);
  
  vector< string > splitLine;
  string line;
  
  if(optimParamsFile.is_open()) {
      
    while(getline(optimParamsFile, line)) {
        
      boost::split(splitLine, line, [](char c){ return c == ' '; });
      
      if(splitLine.size() == 2) { //this is meant to skip the first line 'AIC = X'
          
        string name = splitLine[0];
        double value = atof(splitLine[1].c_str());
        
        inputParameters.addParameter(new Parameter(name, value)); 
      }
    }
    
    optimParamsFile.close();
  }
  
  else {
    throw Exception("iSMC::Unable to open file with optimized parameters!");
  }
  
  return inputParameters;
}
    
void SequentiallyMarkovCoalescent::discretizeTimeUsingExpBoundaries_() {
  
  for(size_t i = 0; i < timeIntervals_.size(); ++i) {
    timeIntervals_[i] = -log(1. - static_cast< double >(i) / static_cast< double >(timeIntervals_.size()));
    //cout << "timeIntervals_[" << i << "] = " << timeIntervals_[i] << endl;
  }
}
  
void SequentiallyMarkovCoalescent::discretizeTimeUsingTmax_(double tMax) {
  
  for(size_t i = 0; i < timeIntervals_.size(); ++i) {
    timeIntervals_[i] = 0.1 * exp(static_cast< double >(i) / static_cast< double >(timeIntervals_.size())
                        * log(1. + 10. * tMax)) - 0.1;
    //cout << "timeIntervals_[" << i << "] = " << timeIntervals_[i] << endl;
  }
}

double SequentiallyMarkovCoalescent::fetchBiHaploidWattersonsTheta_(const vector< unsigned char >& sequence) {
    
  size_t polymorphic = 0;
  size_t callable = 0;
  
  for(size_t i = 0; i < sequence.size(); ++i) {
      
    if(sequence[i] != 2) { //skips missing information
        
      ++callable;
      
      if(sequence[i] == 1) { 
        ++polymorphic;
      }
    }
  }
  
  //cout << "polymorphic = " << polymorphic << "; callable = " << callable << endl;
  return static_cast< double >(polymorphic) / static_cast< double >(callable);
}
  
void SequentiallyMarkovCoalescent::computeMeanThetaAcrossDataSet_(const vector< vector < unsigned char > >& sequences) {
    
  double sumOfThetas = 0.;
  
  for(size_t i = 0; i < sequences.size(); ++i) {
    sumOfThetas += fetchBiHaploidWattersonsTheta_(sequences[i]);
  }
  
  double averageThetaAcrossDataSet = sumOfThetas / static_cast< double >(sequences.size());
  
  //adds a check because an error here is common when testing different ways to parse input seq file
  if(!(averageThetaAcrossDataSet > 0.)) {
    throw Exception("iSMC::Something went wrong when computing theta :(");
  }
  
  setParameterValue("theta", averageThetaAcrossDataSet);
}
