/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 02/07/2018
 *
 */


#include <cmath>

#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Parameter.h>

#include "Splines.h"

using namespace std;
using namespace bpp;


void Splines::fireParameterChanged(const ParameterList& parameters) {
  mapLambdasFromSplines(mmsmc_ -> getLambdaVector());
}

void Splines::mapLambdasFromSplines(ParameterList& Lambdas) {
  
  //maps the first lambda here to simplify the loop below:
  Lambdas.setParameterValue("l0", getParameterValue("y0"));
  
  //for all breakpoints and until tmax is reached 
  //fixedKnots_ are placed BETWEEN t0 and tmax
  for(size_t i = 0; i <= fixedKnots_.size(); ++i) { 
    
    double leftBoundary = 0.;
    if(i > 0) {
      leftBoundary = fixedKnots_[i - 1];
    }
    //cout << "leftBoundary: " << leftBoundary << endl;
    double rightBoundary = mmsmc_ -> getTimeIntervals()[mmsmc_ -> getNumberOfIntervals() - 1];
    if(i < fixedKnots_.size()) {
      rightBoundary = fixedKnots_[i];
    }
    
    //cout << "rightBoundary: " << rightBoundary << endl;
    double leftIntercept = getParameterValue("y" + TextTools::toString(i));
    //cout << "leftIntercept: " << leftIntercept << endl;
    double rightIntercept = getParameterValue("y" + TextTools::toString(i + 1));
    //cout << "rightIntercept: " << rightIntercept << endl;
    double leftDerivative = getParameterValue("y" + TextTools::toString(i) + "_prime");
    double rightDerivative = getParameterValue("y" + TextTools::toString(i + 1) + "_prime");
    
    //the curve in the interval defined by breakpoint i and breakpoint i + 1
    double powerOne = getPowerOne_(leftBoundary, rightBoundary, i, leftIntercept, rightIntercept, leftDerivative, rightDerivative);
    double powerTwo = getPowerTwo_(leftBoundary, rightBoundary, i, leftIntercept, rightIntercept, leftDerivative, rightDerivative);
    double powerThree = getPowerThree_(leftBoundary, rightBoundary, i, leftIntercept, rightIntercept, leftDerivative, rightDerivative);
    double constant = getConstant_(leftBoundary, rightBoundary, i, leftIntercept, rightIntercept, leftDerivative, rightDerivative);
    
    //screens all time intervals that fall between leftBoundary and rightBoundary
    for(size_t j = 1; j < mmsmc_ -> getNumberOfIntervals(); ++j) {
      
      double focalTime = mmsmc_ -> getTimeIntervals()[j];
      
      if(focalTime > leftBoundary) {
        
        if(focalTime <= rightBoundary) {
          
          //cout << "focalTime = " << focalTime << " ";
          //the splines function evaluated at point j (focalTime)
          double linear = powerOne * focalTime;
          double quadratic = powerTwo * pow((focalTime), 2.);
          double cubic = powerThree * pow((focalTime), 3.);
          double lambdaMapping = constant + linear + quadratic + cubic;
          //cout << "lambdaMapping = " << lambdaMapping << endl;
          
          if(lambdaMapping <= 1e-2) { lambdaMapping = 1e-2; }
          Lambdas.setParameterValue("l" + TextTools::toString(j), lambdaMapping);
          
        }
        
        else {
          break;
        }
      }
    }
  }    
}   
  
void Splines::includeSplinesParameters_() {
  //from t = 0. until the last knot:
  for(size_t i = 0; i <= numberOfKnots_; ++i) {
    //intercept values
     addParameter_(new Parameter("y" + TextTools::toString(i), 1. + pow(-1., i) * .01, Parameter::R_PLUS_STAR));
    //first derivatives NOTE the interval constraint should be estimated by computing the range where we can NOT get negative lambda mappings
    addParameter_(new Parameter("y" + TextTools::toString(i) + "_prime", 0., shared_ptr<Constraint>(new IntervalConstraint(-1., 1., true, true))));
  }
  //at tmax:
  addParameter_(new Parameter("y" + TextTools::toString(numberOfKnots_ + 1), 1., Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("y" + TextTools::toString(numberOfKnots_ + 1) + "_prime", 0., shared_ptr<Constraint>(new IntervalConstraint(-1., 1., true, true))));
}
  
double Splines::getPowerOne_(double leftKnot, double rightKnot, size_t index,
                             double leftIntercept, double rightIntercept,
                             double leftDerivative, double rightDerivative) 
{
  return leftDerivative - 6. * leftKnot * rightKnot * (rightIntercept - leftIntercept) /
         pow((rightKnot - leftKnot), 3.) + 3. * leftKnot * rightKnot * (leftDerivative + rightDerivative) /
         pow((rightKnot - leftKnot), 2.) - leftKnot * (rightDerivative - leftDerivative) / (rightKnot - leftKnot);
}
  
double Splines::getPowerTwo_(double leftKnot, double rightKnot, size_t index,
                             double leftIntercept, double rightIntercept,
                             double leftDerivative, double rightDerivative) 
{
  return 1. / 2. * (rightDerivative - leftDerivative) / (rightKnot - leftKnot) + 3. * (rightKnot + leftKnot) *
         (rightIntercept - leftIntercept) / pow((rightKnot - leftKnot), 3.) - 3. / 2. * (rightKnot + leftKnot) *
         (rightDerivative + leftDerivative) / pow((rightKnot - leftKnot), 2.);
}
  
double Splines::getPowerThree_(double leftKnot, double rightKnot, size_t index,
                               double leftIntercept, double rightIntercept,
                               double leftDerivative, double rightDerivative) 
{
  return - 2. * (rightIntercept - leftIntercept) / pow((rightKnot - leftKnot), 3.) + (rightDerivative + leftDerivative) / pow((rightKnot - leftKnot), 2.);
}
  
double Splines::getConstant_(double leftKnot, double rightKnot, size_t index,
                             double leftIntercept, double rightIntercept,
                             double leftDerivative, double rightDerivative) 
{
  return leftIntercept - leftKnot * leftDerivative - pow(leftKnot, 2.) * (rightIntercept - leftIntercept) *
         (leftKnot - 3. * rightKnot) / pow((rightKnot - leftKnot), 3.) + 1. / 2. * pow(leftKnot, 2.) *
         (leftDerivative + rightDerivative) * (leftKnot - 3. * rightKnot) / pow((rightKnot - leftKnot), 2.) +
         1. / 2. * pow(leftKnot, 2.) * (rightDerivative - leftDerivative) / (rightKnot - leftKnot);
}

