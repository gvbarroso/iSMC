/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 08/01/2019
 *
 */


#include "SmcTransitionProbabilities.h"

using namespace bpp;
using namespace std;


void SmcTransitionProbabilities::setUpExpectedMatrix() { 
  Vdouble rowSumVector(smc_ -> getNumberOfIntervals());
  for(size_t i = 0; i < smc_ -> getNumberOfIntervals(); ++i) {
    double timeBeta = smc_ -> getTimeIntervals()[i];
    double rowSum = 0.;
    for(size_t j = 0; j < smc_ -> getNumberOfIntervals(); ++j) {
     double timeAlpha = smc_ -> getTimeIntervals()[j];
      if(i == j) { //main diagonal
        expectedMatrix_[i][j] = 0.; //temporary value
      }
      else {
        if(timeBeta < timeAlpha) { //checks which transition function to use:
          expectedMatrix_[i][j] = computeTransitionProbabilityIncreasedTime_(timeBeta, timeAlpha, smc_ -> getLambdaVector(), smc_ -> getParameterValue("rho"));
        }
        else if(timeBeta > timeAlpha) { 
          expectedMatrix_[i][j] = computeTransitionProbabilityDecreasedTime_(timeBeta, timeAlpha, smc_ -> getLambdaVector(), smc_ -> getParameterValue("rho"));
        }
      }
      rowSum += expectedMatrix_[i][j];
    }
    rowSumVector[i] = rowSum;
  }
  for(size_t k = 0; k < smc_ -> getNumberOfIntervals(); ++k) { 
    expectedMatrix_[k][k] = 1. - rowSumVector[k]; //main diagonal
  }
}

double SmcTransitionProbabilities::exponentiatedIntegral_(double time1, double time2, const ParameterList& lambdaVec) {
  double lowerBoundT1 = smc_ -> fetchLowerTimeBoundary(time1); //lower time boundary for time1 (a.k.a timeBeta)
  double lowerBoundT2 = smc_ -> fetchLowerTimeBoundary(time2); //lower time boundary for time2 (a.k.a timeAlpha)
  double intensity = -1.; //intensity of the exp. decay
  unsigned int beta = smc_ -> fetchTimeIndex(lowerBoundT1); //time index we are transitioning from
  double lambda = lambdaVec[beta].getValue(); //coal. rate for this time interval
  if(lowerBoundT1 == lowerBoundT2) { //if time1 and time2 belong in the same time interval
    intensity = (time2 - time1) * lambda;
  }
  else { //if time1 and time2 belong in different time intervals
    //intensity of the exp. decay is made of three components: c1, c2 and c3
    double upperBoundT1 = smc_ -> fetchUpperTimeBoundary(time1); //upper time boundary of time1
    double c1 = (upperBoundT1 - time1) * lambda; //first component relates to the "residual" time in the interval where time1 belongs.
    unsigned int alpha = smc_ -> fetchTimeIndex(lowerBoundT2);
    double c2 = 0.; //the second component is the rate at every "whole" time interval up to the one
    for(size_t k = beta + 1; k < alpha; ++k) { //from beta until we find the time interval just before alpha:
      double deltaTime = smc_ -> getTimeIntervals()[k + 1] - smc_ -> getTimeIntervals()[k];
      lambda = lambdaVec[k].getValue();
      c2 += deltaTime * lambda;
    }
    lambda = lambdaVec[alpha].getValue(); //last time interval (alpha).
    double c3 = (time2 - lowerBoundT2) * lambda; //third component is the "residual" time in this last time interval 
    intensity = c1 + c2 + c3; //total intensity of the exp. decay is the sum of the three components
  }
  return exp(- intensity); //probability that the haplotypes are still segregating in the ancestral pop.
}
  
double SmcTransitionProbabilities::exponentiatedIntegralM_(double time1, double time2, const ParameterList& lambdaVec) {
  double lowerBoundT1 = smc_ -> fetchLowerTimeBoundary(time1); //lower time boundary for time1 (a.k.a beta)
  double lowerBoundT2 = smc_ -> fetchLowerTimeBoundary(time2); //lower time boundary for time2 (a.k.a alpha)
  double intensity = -1.; //intensity of the exp. decay
  double lambda = -1.; //helper variable to get the coal. rate for each time interval
  if(lowerBoundT1 == lowerBoundT2) { //if time1 and time2 belong in the same time interval
    unsigned int beta = smc_ -> fetchTimeIndex(lowerBoundT1); //in what time interval are we?
    lambda = lambdaVec[beta].getValue() * 2.; //the coal. rate for this time interval
    intensity = (time2 - time1) * lambda;
  }
  else if(lowerBoundT1 != lowerBoundT2) { //if time1 and time2 belong in different time intervals
    //intensity of the exp. decay is made of three components: c1, c2 and c3
    unsigned int beta = smc_ -> fetchTimeIndex(lowerBoundT1); //in what time interval are we?
    lambda = lambdaVec[beta].getValue() * 2.; 
    double upperBoundT1 = smc_ -> fetchUpperTimeBoundary(time1); 
    double c1 = (upperBoundT1 - time1) * lambda; //first component relates to the "residual" time in the interval where time1 belongs
    unsigned int alpha = smc_ -> fetchTimeIndex(lowerBoundT2);
    //now we move from beta until we find the time interval just before alpha:
    double c2 = 0.; //the second component is the rate at every "whole" time interval up to the one just before where time2 belongs (alpha - 1)
    for(size_t k = beta + 1; k < alpha; ++k) {
      double deltaTime = smc_ -> getTimeIntervals()[k + 1] - smc_ -> getTimeIntervals()[k];
      lambda = lambdaVec[k].getValue() * 2.;
      c2 += deltaTime * lambda;
    }
    lambda = lambdaVec[alpha].getValue() * 2.; //Now we are at the last time interval (alpha).
    double c3 = (time2 - lowerBoundT2) * lambda; //The third component is the "residual" time in this last time interval 
    intensity = c1 + c2 + c3; //the total intensity of the exp. decay is just the some of the three components
  }
  return exp(- intensity); //the probability that 2 haplotypes are still segregating in the ancestral pop.
}
  
double SmcTransitionProbabilities::computeTransitionProbabilityDecreasedTime_(double timeBeta, double timeAlpha, const ParameterList& lambVec, double rho) {
  //We are calculating the transtition prob FROM timeBeta TO timeAlpha and both are time boundaries
  //unsigned int indexBeta = smc_ -> getTimeIndex(timeBeta);
  unsigned int indexAlpha = smc_ -> fetchTimeIndex(timeAlpha);
  double lambda_a = lambVec[indexAlpha].getValue(); //in the PSMC' we have "small lambda" == "big lambda"
  //the coalescence rate for recombining branch m in time interval alpha (including itself):
  double lambda_m_a = 2. * lambVec[indexAlpha].getValue();
  double avgTimeBeta = smc_ -> getAverageCoalescenceTime(timeBeta); //NOTE Schiffels uses timeBeta here
  //double avgTimeAlpha = smc_ -> getAverageCoalescenceTime(timeAlpha); 
  //devides the formula in components to properly fit it to the script:
  //c1 is the term outside parenthesis in Schiffels & Durbin Supp. Note (replacing M by 2 here):
  double c1 = (1. - exp(- 2. * rho * avgTimeBeta)) * (1. / avgTimeBeta) * (1. / 2.) * lambda_a; //NOTE avgtimeBeta
  //double c1 = (1. - exp(- 2. * rho * avgTimeAlpha)) * (1. / avgTimeAlpha) * (1. / 2.) * lambda_a; 
  double deltaTime_a = -1.; //delta time for alpha (ie, how much time there is in that particular time interval)
  if(indexAlpha == smc_ -> getNumberOfIntervals() - 1) { 
    deltaTime_a = numeric_limits< double >::max(); //if we are at the last time interval, time spent in it is "infinite":
  } 
  else { 
    deltaTime_a = smc_ -> getTimeIntervals()[indexAlpha + 1] - timeAlpha;
  }
  double c2 = 0.; //c2 is the prob. of 2 segregating sequences until the relevant time interval
  for(size_t g = 0; g < indexAlpha; ++g) { //for every time interval from 0 up to alpha - 1
    double tempUpperGamma = smc_ -> getTimeIntervals()[g + 1]; 
    double deltaTime_g = smc_ -> getTimeIntervals()[g + 1] - smc_ -> getTimeIntervals()[g];
    double lambda_m_g = 2. * lambVec[g].getValue(); //the coalescence rate in time interval g, after recombination
    double L_m = exponentiatedIntegralM_(tempUpperGamma, timeAlpha, lambVec);
    //further divides the c2 component to properly fit the script
    double c2a = (1. / lambda_m_g) * (1. - exp(- deltaTime_g * lambda_m_g)); 
    double c2b = (1. / lambda_m_a) * (1. - exp(- deltaTime_a * lambda_m_a));
    c2 += 2. * c2a * L_m * c2b; //sum over two possibilities (recombining branch being either i or j), thus the factor of 2:
  }
  //c3 (ie, once we are in alpha; the factor of two again present):
  double c3 = (2. / lambda_m_a) * (deltaTime_a - (1. / lambda_m_a) * (1. - exp(- deltaTime_a * lambda_m_a)));
  //There are numerical issues with c3 when lambda << 1; therefore, we replace c3 by deltaTime_a ^ 2 (its limit as lambda_m_a goes to zero):
  //f(l) := 2/l * (T - (1/l *(1-exp(-T*l)))); f(l) := 2/l * (T - (1/l *(1-exp(-T*l)))); | in mxMaxima
  if(c3 < 0) { 
    c3 = pow(deltaTime_a, 2.);
  }
  return c1 * (c2 + c3);
}
  
double SmcTransitionProbabilities::computeTransitionProbabilityIncreasedTime_(double timeBeta, double timeAlpha, const ParameterList& lambVec, double rho) {
  //We are calculating the transtition prob FROM timeBeta TO timeAlpha and both are time boundaries
  unsigned int indexBeta = smc_ -> fetchTimeIndex(timeBeta);
  unsigned int indexAlpha = smc_ -> fetchTimeIndex(timeAlpha);
  double lambda_a = lambVec[indexAlpha].getValue(); //coalescence rate in time interval alpha:
  double deltaTime_a = -1.; //delta time for alpha (ie, how much time there is in that smcticular time interval)
  if(indexAlpha == smc_ -> getNumberOfIntervals() - 1) { //if we are at the last time interval, time spent in it is "infinite":
    deltaTime_a = numeric_limits< double >::max();
  } 
  else { //if we are not, time spent in the interval is the upper interval - our interval
    deltaTime_a = smc_ -> getTimeIntervals()[indexAlpha + 1] - timeAlpha;
  }
  double avgTimeBeta = smc_ -> getAverageCoalescenceTime(timeBeta);
  //NOTE: avgTimeBeta spans from t = 0.0 until the average coal. time inside timeBeta
  double L_b_a = exponentiatedIntegral_(avgTimeBeta, timeAlpha, lambVec); //probability of segregating pair of lineages from timeBeta to timeAlpha
  //c1 (divided in two to fit the script):
  //NOTE: in the case of the PSMC, the lower case and upper case lambdas are the same and cancel out
  double c1a = (1. - exp(- 2. * rho * avgTimeBeta)) * (1. / avgTimeBeta) * (1. / 2.) * lambda_a * L_b_a;
  double c1b = (1. / lambda_a) * (1. - exp(- deltaTime_a * lambda_a));
  double c1 = c1a * c1b;
  double c2 = 0.;
  //If there are events to be accounted for before we reach time interval beta.
  //This events correspond to the probability of no coalescence until timeBeta for each time interval where recombination could have occurred. 
  if(indexBeta != 0) {
    for(size_t g = 0; g < indexBeta; ++g) {
      double tempUpperGamma = smc_ -> getTimeIntervals()[g + 1]; 
      double deltaTime_g = smc_ -> getTimeIntervals()[g + 1] - smc_ -> getTimeIntervals()[g];
      double lambda_m_g = 2. * lambVec[g].getValue(); //the coalescence rate for recombining branch m (including itself) in time interval g
      //recombining branch m can be either one of the 2 existing branches (hence the factor of 2)
      double L_m = exponentiatedIntegralM_(tempUpperGamma, avgTimeBeta, lambVec); //probability of no coal. since time of recomb. until average coal. time in beta:
      //sum over two possibilities (recombining branch being either i or j), thus the factor of 2
      c2 += 2. * (1. / lambda_m_g) * (1. - exp(- deltaTime_g * lambda_m_g)) * L_m; 
    }
  }
  double lambda_m_b = 2. * lambVec[indexBeta].getValue(); //recombining branch m can be either one of the 2 existing branches (again the factor of 2):
  double c3 = 2. * (1. / lambda_m_b) * (1. - exp(-(avgTimeBeta - timeBeta) * lambda_m_b)); //here also the factor of 2:
  return c1 * (c2 + c3);
}
