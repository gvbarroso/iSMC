/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 04/01/2019
 *
 */


#include <numeric>

#include "BaumWelch.h"

//#define DEBUG_BW

using namespace std;
using namespace bpp;


//designed for HMM layers
void BaumWelch::maximiseRateTransMat(const Vdouble& fwdScales, const VVdouble& fwdMat, const VVdouble& bckMat,
                                     VVdouble& transMat, const VVdouble& emissMat, //TM and EM for rate
                                     const vector< vector< size_t > >& obsVec) {
  
  size_t sequenceLength = fwdMat.size();
  size_t numDiploids = obsVec.size();
  size_t numCategories = transMat.size();
  
  //cout << "dimensions of rate emissMat: " << emissMat.size() << "; x " << emissMat.front().size() << endl;
  //cout << "length of obs. seq. = " << obsVec.front().size() << endl;
  
  VVdouble transitionCounts(numCategories, Vdouble(numCategories)); //stores expected counts j -> k:
  for(size_t i = 0; i < numCategories; ++i) { 
    transitionCounts[i].assign(numCategories, 0.);
  }

  //computes expected transiton counts
  for(size_t i = 1; i < sequenceLength; ++i) { //no transition in the init. of chain
    for(size_t j = 0; j < numCategories; ++j) {
      for(size_t k = 0; k < numCategories; ++k) {
        double compEmissProb = 1.; //composite emission probability (over all diploids)
        for(size_t l = 0; l < numDiploids; ++l) {
          size_t obs = obsVec[l][i];
          compEmissProb *= emissMat[k][obs];
        }  
        double trans = transMat[j][k];
        double fwd = fwdMat[i - 1][j];
        double bck = bckMat[i][k];
        transitionCounts[j][k] += fwd * trans * compEmissProb * bck * (1. / fwdScales[i]); 
      }
    }
  }    

  //normalises
  for(size_t j = 0; j < numCategories; ++j) {
    //from hidden state j to all others
    double tmpSum = accumulate(transitionCounts[j].begin(), transitionCounts[j].end(), 0.);
    for(size_t k = 0; k < numCategories; ++k) {
      transitionCounts[j][k] /= tmpSum;
    }
  }    

  #ifdef DEBUG_BW
  cout << "Computing Baum-Welch Transition Matrix row sums..." << endl;
  for(size_t i = 0; i < transitionCounts.size(); ++i) {
    Vdouble temp = transitionCounts[i];
    double rowSum = accumulate(temp.begin(), temp.end(), 0.);
    if(!(abs(1. - rowSum) < 1e-9)) {
      cout << "Sum of entries for row " << i << " = " << setprecision(20) << rowSum << endl;
    }
  }
  cout << "   done." << endl;
  #endif

  transMat = transitionCounts; //updates TM
}

//designed for HMM layers
void BaumWelch::maximiseRateEmissMat(const VVdouble& fwdMat, const VVdouble& bckMat, VVdouble& rateEmissMat,
                                     const vector< vector< size_t > >& obsVec) {
  size_t sequenceLength = fwdMat.size();
  size_t numDiploids = obsVec.size();
  size_t numCategories = rateEmissMat.size(); 
  size_t numObsStates = rateEmissMat.front().size();
    
  VVdouble emissionCounts(numCategories, Vdouble(numObsStates)); //stores expected counts j -> k:
  for(size_t i = 0; i < numCategories; ++i) { 
    emissionCounts[i].assign(numObsStates, 0.);
  }
    
  for(size_t i = 0; i < sequenceLength; ++i) { 
    for(size_t j = 0; j < numCategories; ++j) {
      for(size_t l = 0; l < numDiploids; ++l) {
        size_t obs = obsVec[l][i]; //observation in diploid l at pos. i
        emissionCounts[j][obs] += fwdMat[i][j] * bckMat[i][j];
      }
    }
  }

  //normalises
  for(size_t j = 0; j < numCategories; ++j) {
    //from hidden state j to all observations (tree transitions or tree emissions)
    double tmpSum = accumulate(emissionCounts[j].begin(), emissionCounts[j].end(), 0.); 
    for(size_t k = 0; k < numObsStates; ++k) { 
      emissionCounts[j][k] /= tmpSum; 
    }
  }

  #ifdef DEBUG_BW
  cout << "Computing Baum-Welch Emission Matrix row sums..." << endl;
  for(size_t i = 0; i < emissionCounts.size(); ++i) {
    Vdouble temp = emissionCounts[i];
    double rowSum = accumulate(temp.begin(), temp.end(), 0.);
    if(!(abs(1. - rowSum) < 1e-9)) { 
      cout << "Sum of entries for row " << i << " = " << setprecision(20) << rowSum << endl;
    }
  }
  cout << "   done." << endl;
  #endif

  rateEmissMat = emissionCounts;
}

//using data from a single diploid
VVdouble BaumWelch::proposeTransMat(const Vdouble& fwdScales, const VVdouble& fwdMat, const VVdouble& bckMat,
                                    const VVdouble& transMat, const VVdouble& emissMat,
                                    const vector< size_t >& seq) {

  size_t sequenceLength = fwdMat.size();
  size_t numHiddenStates = transMat.size();

  VVdouble transitionCounts(numHiddenStates, Vdouble(numHiddenStates)); //stores expected counts j -> k:
  for(size_t i = 0; i < numHiddenStates; ++i) { 
    transitionCounts[i].assign(numHiddenStates, 0.);
  }

  //computes expected transiton counts
  for(size_t i = 1; i < sequenceLength; ++i) { //no transition in the init. of chain
    size_t obs = seq[i];
    for(size_t j = 0; j < numHiddenStates; ++j) {
      for(size_t k = 0; k < numHiddenStates; ++k) {
        double emiss = emissMat[k][obs]; //conditioning on obs at site i
        double trans = transMat[j][k];
        double fwd = fwdMat[i - 1][j];
        double bck = bckMat[i][k];
        transitionCounts[j][k] += fwd * trans * emiss * bck * (1. / fwdScales[i]); 
      }
    }
  }    

  //normalises
  for(size_t j = 0; j < numHiddenStates; ++j) {
    //from hidden state j to all others
    double tmpSum = accumulate(transitionCounts[j].begin(), transitionCounts[j].end(), 0.);
    for(size_t k = 0; k < numHiddenStates; ++k) {
      transitionCounts[j][k] /= tmpSum;
    }
  }

  return transitionCounts;
}

//using data from a single diploid
VVdouble BaumWelch::proposeEmissMat(const VVdouble& fwdMat, const VVdouble& bckMat, const VVdouble& emissMat, const vector< size_t >& seq) {

  size_t sequenceLength = fwdMat.size();
  size_t numHiddenStates = emissMat.size(); 
  size_t numObsStates = emissMat.front().size();

  VVdouble emissionCounts(numHiddenStates, Vdouble(numObsStates)); //stores expected counts j -> k:
  for(size_t i = 0; i < numHiddenStates; ++i) { 
    emissionCounts[i].assign(numObsStates, 0.);
  }
    
  for(size_t i = 0; i < sequenceLength; ++i) {
    size_t obs = seq[i];  
    for(size_t j = 0; j < numHiddenStates; ++j) {
      //only updates counts for obs at site i
      emissionCounts[j][obs] += fwdMat[i][j] * bckMat[i][j];
    }
  }
  
  //normalises
  for(size_t j = 0; j < numHiddenStates; ++j) {
    //from hidden state j to all observations
    double tmpSum = accumulate(emissionCounts[j].begin(), emissionCounts[j].end(), 0.); 
    for(size_t k = 0; k < 2; ++k) { //don't update emissions to missing data
      emissionCounts[j][k] /= tmpSum; 
    }
    emissionCounts[j][2] = 1.; //assigns emissions to missing data
  }

  #ifdef DEBUG_BW
  cout << "Computing Baum-Welch Emission Matrix row sums..." << endl;
  for(size_t i = 0; i < emissionCounts.size(); ++i) {
    Vdouble temp = emissionCounts[i];
    double rowSum = accumulate(temp.begin(), temp.end(), 0.);
    if(!(abs(2. - rowSum) < 1e-9)) { //2 = 1 + 1 from missing data
      cout << "Sum of entries for row " << i << " = " << setprecision(20) << rowSum << endl;
    }
  }
  cout << "done." << endl;
  #endif

  return emissionCounts;
}

//to perform EM on fragments (when "smoothing" demography)
void BaumWelch::incrementTransitionMatrix(const zipHMM::Matrix& fwdMat, const zipHMM::Matrix& bckMat,
                                          VVdouble& newTransMat, const VVdouble& oldTransMat,
                                          const VVdouble& emissMat, const Vdouble& fwdScales,
                                          const vector< size_t >& frag) {

  size_t sequenceLength = frag.size();
  size_t numHiddenStates = oldTransMat.size();

  //computes expected transiton counts for this fragment
  for(size_t i = 1; i < sequenceLength; ++i) { //no transition in the init. of chain WARNING this is a problem when doing it in fragments (use piMmSmc)
    size_t obs = frag[i];
    for(size_t j = 0; j < numHiddenStates; ++j) {
      for(size_t k = 0; k < numHiddenStates; ++k) {
        double emiss = emissMat[k][obs]; //conditioning on obs at site i
        double trans = oldTransMat[j][k];
        double fwd = fwdMat(j, i - 1);
        double bck = bckMat(k, i);
        //increments transition
        newTransMat[j][k] += fwd * trans * emiss * bck * (1. / fwdScales[i]); 
      }
    }
  }
  //normalisation should be done after all fragments have been surveyed
}
