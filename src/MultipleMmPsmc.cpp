/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 18/03/2019
 *
 */


#include <string>
#include <vector>
#include <thread>
#include <cmath>
#include <algorithm>

#include "Global.h"
#include "MultipleMmPsmc.h"
#include "MmPsmc.h"

using namespace std;
using namespace bpp;


vector< shared_ptr< MmPsmc > > MultipleMmPsmc::fetchTrainingPsmcVector() {

  vector< shared_ptr< MmPsmc > > trainingPsmcVector(0);

  for(size_t i = 0; i < psmcVec_.size(); ++i) {
    if(i != testIndex_) {
      trainingPsmcVector.push_back(psmcVec_[i]);
    }
  }

  return trainingPsmcVector;
}

void MultipleMmPsmc::updatePsmcBackwardMatrices() {
    
  size_t numberOfUsedThreads = min({ NUMBER_OF_AVAILABLE_THREADS, psmcVec_.size() }); 
  size_t numberOfTasks = psmcVec_.size();
  size_t chunkSize = (numberOfTasks + NUMBER_OF_AVAILABLE_THREADS - 1) / NUMBER_OF_AVAILABLE_THREADS; //ceiling
  
  auto backwardAlgorithmParallel = [&] (size_t thread_id) {
      
    size_t from = chunkSize * thread_id;
    size_t to = chunkSize * (thread_id + 1);
    
    for(size_t j = from; j < to; ++j) {
        
      if(j < psmcVec_.size()) {
        psmcVec_[j] -> backwardAlgorithm();
      }
    }
    
  };

  thread* threadVector = new thread[numberOfUsedThreads];
  
  for(size_t i = 0; i < numberOfUsedThreads; ++i) {
    threadVector[i] = thread(backwardAlgorithmParallel, i);
  }    
  
  for(size_t i = 0; i < numberOfUsedThreads; ++i) {
    threadVector[i].join();
  }   
    
  delete [] threadVector;
}
 
double MultipleMmPsmc::fetchCompositeLogLikelihood() {
    
  computeIndividualsLogLikelihood_();
  
  double compositeLogLikelihood = 0.;
  
  for(size_t i = 0; i < psmcVec_.size(); ++i) {
      
    if(i != testIndex_) {
      compositeLogLikelihood += psmcVec_[i] -> getBiHaploidLogLikelihood();
    }
  }
  
  return compositeLogLikelihood;
}
  
double MultipleMmPsmc::fetchMedianLogLikelihood() {
    
  computeIndividualsLogLikelihood_();
  
  Vdouble individualLikelihoods(0);
  
  for(size_t i = 0; i < psmcVec_.size(); ++i) {
      
    if(i != testIndex_) {
      individualLikelihoods.push_back(psmcVec_[i] -> getBiHaploidLogLikelihood());
    }
  }
  
  double median = -1.;
  
  if(individualLikelihoods.size() == 1) { 
    median = individualLikelihoods[0];
  }
  
  else {
      
    sort(individualLikelihoods.begin(), individualLikelihoods.end()); 
    
    if(individualLikelihoods.size() % 2 == 0) {
      size_t medianIndex = individualLikelihoods.size() / 2;
      median = (individualLikelihoods[medianIndex] + individualLikelihoods[medianIndex - 1]) / 2.;
    }
    
    else if(individualLikelihoods.size() % 2 != 0) {
      size_t medianIndex = (individualLikelihoods.size() - 1) / 2;
      median = individualLikelihoods[medianIndex];
    }
  }
  
  return median;
}
  
double MultipleMmPsmc::fetchMinimumLogLikelihood() {
    
  computeIndividualsLogLikelihood_();
  
  Vdouble individualLikelihoods(0);
  
  for(size_t i = 0; i < psmcVec_.size(); ++i) {
    if(i != testIndex_) {
      individualLikelihoods.push_back(psmcVec_[i] -> getBiHaploidLogLikelihood());
    }
  } 
  
  sort(individualLikelihoods.begin(), individualLikelihoods.end()); 
  
  return individualLikelihoods.front();
}

void MultipleMmPsmc::jointForward(const vector< vector< size_t > >& obsVector, const VVdouble& rateTrans, const Vdouble& vpi,
                                  const VVdouble& rateEmiss, Vdouble& forwardScales, VVdouble& rateForwardMatrix) {
     
 // NOTE Should we weight emission probs by TMRCA of each diploid? By 1 - shannon_eq of each diploid?
  
  size_t numCategories = rateTrans.size(); //no. hidden states in HMM layers
  size_t numDiploids = obsVector.size();
  size_t seqLength = obsVector.front().size();
    
  //init
  double scale = 0.;
  for(size_t j = 0; j < numCategories; ++j) {
      
    //cout << "j: " << j << "; ";  
    double jointEmissions = 1.;
    for(size_t l = 0; l < psmcVec_.size(); ++l) { 
        
      size_t obs = obsVector[l][0];  
      jointEmissions *= rateEmiss[j][obs];
    }
    
    //cout << "joint-emiss = " << setprecision(16) << jointEmissions << "; ";
    rateForwardMatrix[0][j] = vpi[j] * jointEmissions;
    //cout << "fwd = " << setprecision(16) << rateForwardMatrix[0][j] << endl;
    scale += rateForwardMatrix[0][j];
  }
  
  //cout << endl << "scale = " << scale << endl << endl;
  forwardScales[0] = scale;
  
  for(size_t j = 0; j < numCategories; ++j) {
    rateForwardMatrix[0][j] /= scale;
  }
  
  //recursion 
  for(size_t i = 1; i < seqLength; ++i) {
      
    scale = 0.; //restarts scale:
    
    for(size_t j = 0; j < numCategories; ++j) {
        
      double a_j = 0.; //a_j(i)
      
      for(size_t k = 0; k < numCategories; ++k) {
        a_j += rateForwardMatrix[i - 1][k] * rateTrans[k][j];
      }
      
      double jointEmissions = 1.; //composite emission prob in all diploids
      
      for(size_t l = 0; l < numDiploids; ++l) {
        size_t obs = obsVector[l][i];
        jointEmissions *= rateEmiss[j][obs];
      }
      
      rateForwardMatrix[i][j] = a_j * jointEmissions;
      scale += rateForwardMatrix[i][j];
    }
    
    for(size_t j = 0; j < numCategories; ++j) {
      rateForwardMatrix[i][j] /= scale;
    }
    
    forwardScales[i] = scale; 
    //cout << scale << " ";
  }
  
  /* TESTING
  double ll = 0.;
  for(vector< double >::iterator it = forwardScales.begin(); it != forwardScales.end(); ++it) {
    ll += log(*it);
  }
  cout << "ll forward = " << setprecision(16) << ll << endl; */

}

void MultipleMmPsmc::jointBackward(const vector< vector< size_t > >& obsVector, const VVdouble& rateTrans,
                                   const VVdouble& rateEmiss, const Vdouble& forwardScales, VVdouble& rateBackwardMatrix) {
    
/* 
 * Should we weight emission probs by TMRCA of each diploid? By 1 - shannon_eq of each diploid?
 */  

  size_t numCategories = rateTrans.size(); //no. hidden states in HMM layers
  size_t numDiploids = obsVector.size();
  size_t seqLength = obsVector.front().size();

  //init
  for(size_t j = 0; j < numCategories; ++j) {
    rateBackwardMatrix[obsVector.front().size() - 1][j] = 1.;
  }

  //recursion 
  for(size_t i = seqLength - 1; i > 0; --i) {
      
    for(size_t j = 0; j < numCategories; ++j) {
        
      double b_j = 0.; //b_j(i-1) 
      
      for(size_t k = 0; k < numCategories; ++k) {
          
        double jointEmissions = 1.0; //composite emission prob in all diploids
        
        for(size_t l = 0; l < numDiploids; ++l) { 
            
          size_t obs = obsVector[l][i];  
          jointEmissions *= rateEmiss[k][obs];
        }
        
        b_j += rateTrans[j][k] * jointEmissions * rateBackwardMatrix[i][k];
      }
      
      rateBackwardMatrix[i - 1][j] = b_j / forwardScales[i]; 
    }
  }
}

void MultipleMmPsmc::computeBatchLogLikelihood_(const vector< shared_ptr< MmPsmc > >& focalPsmcVec, size_t numAvailThreads) {
    
  size_t numberOfTasks = focalPsmcVec.size();
    
  auto likelihoodParallel = [&] (size_t thread_id) {
    focalPsmcVec[thread_id] -> computeBiHaploidLogLikelihood(numAvailThreads); 
  };

  thread* threadVector = new thread[numberOfTasks];
    
  for(size_t i = 0; i < numberOfTasks; ++i) {
    threadVector[i] = thread(likelihoodParallel, i);
  }
  
  for(size_t i = 0; i < numberOfTasks; ++i) {
    threadVector[i].join();
  }
  
  delete [] threadVector;
} 

void MultipleMmPsmc::computeIndividualsLogLikelihood_() {
    
  size_t numAvailThreads = NUMBER_OF_AVAILABLE_THREADS; //resets number of free threads
  
  //create batches 
  size_t numDiploids = psmcVec_.size();
  size_t batchSize = min(numDiploids, numAvailThreads);
  size_t numBatches = (numDiploids + batchSize - 1) / batchSize;
  
  //updates numAvailThreads to be used within each diploid (Psmc)
  numAvailThreads /= batchSize; //will be forwarded down to seq fragments
  
  //computes likelihood per batch
  for(size_t i = 0; i < numBatches; ++i) {
      
    vector< shared_ptr< MmPsmc > > batchPsmcVec(0);
    
    for(size_t j = 0; j < batchSize; ++j) {
      
      size_t index = j + i * batchSize;
            
      //because last batch may not have batchSize PSMC's
      if(index < numDiploids) {
        batchPsmcVec.push_back(psmcVec_[index]); 
      }
      
      else {
        break;
      }
    }
    
    computeBatchLogLikelihood_(batchPsmcVec, numAvailThreads);
  }
}
