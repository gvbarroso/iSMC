/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 02/01/2019
 *
 */


#include <thread>

#include "Global.h"
#include "MatrixOptimizationFunction.h"

using namespace std;
using namespace bpp;


void MatrixOptimizationFunction::fireParameterChanged(const ParameterList& parameters) {
  if(parameters.getCommonParametersWith(mmsmc_ -> getParameters()).size() > 0) {
    mmsmc_ -> matchParametersValues(parameters);
    mmsmc_ -> computeAverageCoalescenceTime();
  }
    
  for(size_t i = 0; i < mmsmc_ -> getParameterScalings().size(); ++i) {
    if(parameters.getCommonParametersWith(mmsmc_ -> getParameterScalings()[i] -> getParameters()).size() > 0) {
      mmsmc_ -> getParameterScalings()[i] -> matchParametersValues(parameters);
      mmsmc_ -> getParameterScalings()[i] -> discretize();
    }
  }
  
  for(size_t j = 0; j < mmsmc_ -> getParameterTransitions().size(); ++j) {
    if(parameters.getCommonParametersWith(mmsmc_ -> getParameterTransitions()[j] -> getParameters()).size() > 0) {
      mmsmc_ -> getParameterTransitions()[j] -> matchParametersValues(parameters);
      mmsmc_ -> getParameterTransitions()[j] -> setUpCategoryTransitionMatrix();
    }
  }

  hmi_ -> setUpExpectedMatrix();
  totalDistance_ = computeMedianMatrixDistance_();
}

void MatrixOptimizationFunction::computeAllMatrixDistances_() {

  size_t numberOfTasks = referenceMatrices_.size();

  auto computeMatrixDist = [&] (size_t thread_id) {
    matrixDistancesVector_[thread_id] = computeSingleMatrixDistance_(referenceMatrices_[thread_id]);
  };
   
  thread* threadVector = new thread[numberOfTasks];
  for(size_t i = 0; i < numberOfTasks; ++i) {
    threadVector[i] = thread(computeMatrixDist, i);
  }
   
  for(size_t i = 0; i < numberOfTasks; ++i) {
    threadVector[i].join();
  }
  
  delete [] threadVector; 
}   
    
double MatrixOptimizationFunction::computeSingleMatrixDistance_(const VVdouble& focalMatrix) {

  computeAllMatrixDistances_();

  double singleMatrixDistance = 0.;
  for(size_t i = 0; i < hmi_ -> getExpectedMatrix().size(); ++i) {
    for(size_t j = 0; j < hmi_ -> getExpectedMatrix()[i].size(); ++j) {
      //since we expect these differences to be small, we use the log:
      if(hmi_ -> getExpectedMatrix()[i][j] != 0.) { //dodges 0. entries in TM
        double d1 = log(focalMatrix[i][j]); //BW proposed matrix
        double d2 = log(hmi_ -> getExpectedMatrix()[i][j]); //"parameter-based" matrix
        singleMatrixDistance += pow((d1 - d2), 2.);
      }
    }
  }
  return singleMatrixDistance;
}
  
double MatrixOptimizationFunction::computeSumOfMatrixDistances_() {

  computeAllMatrixDistances_();

  double sumOfDistances = 0.;
  for(size_t i = 0; i < matrixDistancesVector_.size(); ++i) {
    sumOfDistances += matrixDistancesVector_[i];
  }
  return sumOfDistances;
}
  
double MatrixOptimizationFunction::computeMedianMatrixDistance_() {

  computeAllMatrixDistances_();

  double median = -1.;
  if(matrixDistancesVector_.size() == 1) { 
    median = matrixDistancesVector_[0];
  }
  else {
    sort(matrixDistancesVector_.begin(), matrixDistancesVector_.end()); 
    if(matrixDistancesVector_.size() % 2 == 0) {
      size_t medianIndex = matrixDistancesVector_.size() / 2;
      median = (matrixDistancesVector_[medianIndex] + matrixDistancesVector_[medianIndex - 1]) / 2.;
    }
    else if(matrixDistancesVector_.size() % 2 != 0) {
      size_t medianIndex = (matrixDistancesVector_.size() - 1) / 2;
      median = matrixDistancesVector_[medianIndex];
    }
  }
  return median;
}
  
