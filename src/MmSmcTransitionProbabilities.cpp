/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 16/01/2019
 *
 */


#include "MmSmcTransitionProbabilities.h"

#include <numeric>

//#define DEBUG_TM


using namespace bpp;
using namespace std;
 
  
void MmSmcTransitionProbabilities::setUpExpectedMatrix() { 

  vector< vector< unsigned char > > hiddenStates = mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates();
  //vector to deal with transitions between same time intervals (genealogy.1 == genealogy.2)
  Vdouble rowSumVector(mmsmc_ -> getNumberOfHiddenStates(), 0.);
  //the sate we are transitioning FROM:
  for(size_t i = 0; i < mmsmc_ -> getNumberOfHiddenStates(); ++i) {
    Vdouble focalRowEntries(mmsmc_ -> getNumberOfHiddenStates()); //TM entries for row i
    unsigned char timeI = hiddenStates[i][0];
    unsigned char thetaI = hiddenStates[i][1]; 
    unsigned char rhoI = hiddenStates[i][2];
    unsigned char neI = hiddenStates[i][3];  
    //the sate we are transitioning TO:
    for(size_t j = 0; j < mmsmc_ -> getNumberOfHiddenStates(); ++j) {
      unsigned char timeJ = hiddenStates[j][0];
      unsigned char thetaJ = hiddenStates[j][1];
      unsigned char rhoJ = hiddenStates[j][2];
      unsigned char neJ = hiddenStates[j][3];    
      //transition probabilities between parameter categories
      double thetaTransProb = mmsmc_ -> getParameterTransitions()[0] -> getCategoryTransitions()[thetaI][thetaJ];
      double rhoTransProb = mmsmc_ -> getParameterTransitions()[1] -> getCategoryTransitions()[rhoI][rhoJ];
      double neTransProb = mmsmc_ -> getParameterTransitions()[2] -> getCategoryTransitions()[neI][neJ];
      vector< bool > stateComparison = compareHiddenStates_(hiddenStates[i], hiddenStates[j]);
      //total transition probability for spatial rates
      double spaceTransProb = thetaTransProb * rhoTransProb * neTransProb;
      //cout << "rho TP: " << rhoTransProb << endl;
      //cout << "total STP: " << spaceTransProb << endl;
      //size_t numChangedParams = static_cast< size_t >(count_if(begin(stateComparison), end(stateComparison), [](size_t sameState){ return sameState == false; }));
      //cout << "HS: " << i << " -> " << j << "; # diff. = " << numChangedParams << endl;
      if(i == j) { //main diagonal
        //cout << "main diagonal" << endl;  
        expectedMatrix_[i][j] = 0.; //we fill in later
      }
      else {
        //if the trees remain the same (NOT main diagonal) there is no equation for the transition probability
        if(timeI == timeJ) { 
          //cout << "same tree" << endl;
          expectedMatrix_[i][j] = 0.; //we fill in later
        }
        else {
          double time1 = mmsmc_ -> getTimeIntervals()[timeI];
          double time2 = mmsmc_ -> getTimeIntervals()[timeJ];
          double focalRho = mmsmc_ -> getParameterValue("rho");
          //if modulated by rho, scales localRho
          if(mmsmc_ -> getParameterTransitions()[1] -> getNumberOfCategories() > 1) {
            //Due to a discrepancy between # categories inside ParamScalings & ParamTransitions, the Gamma+Hotspot model is handled differently  
            if(mmsmc_ -> getParameterTransitions()[1] -> getHeterogeneousRateModel() == "Gamma+Hotspot") {
              if(rhoI == 0) { //hotspot
                focalRho *= mmsmc_ -> getParameterScalings()[1] -> getParameterValue("heat");
              }
              else { //gamma categories (-1 to acess the correct indices inside the gamma dist.)
                focalRho *= mmsmc_ -> getParameterScalings()[1] -> getCategories()[rhoI - 1];
              }
            }
            else { //if gamma OR hotspot dist.
              focalRho *= mmsmc_ -> getParameterScalings()[1] -> getCategories()[rhoI];
            }
          }
          if(timeI < timeJ) { 
            expectedMatrix_[i][j] = spaceTransProb * computeTransitionProbabilityIncreasedTime_(time1, time2, mmsmc_ -> getLambdaVector(), focalRho);
            //cout << " scaled TP[" << i << "][" << j << "]: " << expectedMatrix_[i][j] << endl;
          }
          else if(timeI > timeJ) {
            expectedMatrix_[i][j] = spaceTransProb * computeTransitionProbabilityDecreasedTime_(time1, time2, mmsmc_ -> getLambdaVector(), focalRho);
            //cout << " scaled TP[" << i << "][" << j << "]: " << expectedMatrix_[i][j] << endl;
          }
        }
      }
      //cout << "TM[" << i << "][" << j << "]: " << expectedMatrix_[i][j] << "; ";
      focalRowEntries[j] = expectedMatrix_[i][j];
    }
    //cout << endl;
    //to increase accuracy
    sort(focalRowEntries.begin(), focalRowEntries.end());
    rowSumVector[i] = accumulate(focalRowEntries.begin(), focalRowEntries.end(), 0.); 
  }
  //now we fill entries where the hidden states share the same tree:
  setTransitionProbabilitiesEqualTime_(rowSumVector, hiddenStates);
  
  bool restrictTrans = false;
  if(restrictTrans) { //WARNING ugly as fuck
    Vdouble rowComplementVector(mmsmc_ -> getNumberOfHiddenStates(), 0.);
    //the sate we are transitioning FROM:
    for(size_t i = 0; i < mmsmc_ -> getNumberOfHiddenStates(); ++i) {
      Vdouble focalRowComplement(mmsmc_ -> getNumberOfHiddenStates(), 0.);
      //the sate we are transitioning TO:
      for(size_t j = 0; j < mmsmc_ -> getNumberOfHiddenStates(); ++j) {
        vector< bool > stateComparison = compareHiddenStates_(hiddenStates[i], hiddenStates[j]);
        size_t numChangedParams = static_cast< size_t >(count_if(begin(stateComparison), end(stateComparison), [](size_t sameState){ return sameState == false; }));
        if(numChangedParams > 1) {
          focalRowComplement[j] = expectedMatrix_[i][j];
          expectedMatrix_[i][j] = 0.;
        }
      }
      sort(begin(focalRowComplement), end(focalRowComplement));
      rowComplementVector[i] = accumulate(begin(focalRowComplement), end(focalRowComplement), 0.);
    }
    for(size_t i = 0; i < mmsmc_ -> getNumberOfHiddenStates(); ++i) {
      //the sate we are transitioning TO:
      for(size_t j = 0; j < mmsmc_ -> getNumberOfHiddenStates(); ++j) {  
        vector< bool > stateComparison = compareHiddenStates_(hiddenStates[i], hiddenStates[j]);
        size_t numChangedParams = static_cast< size_t >(count_if(begin(stateComparison), end(stateComparison), [](size_t sameState){ return sameState == false; }));
        if(numChangedParams < 2) {
          expectedMatrix_[i][j] /= (1. - rowComplementVector[i]);
        }
      }
    }
  }
  
  #ifdef DEBUG_TM
  cout << "Computing TM row sums..." << endl;
  for(size_t i = 0; i < expectedMatrix_.size(); ++i) {
    Vdouble temp = expectedMatrix_[i];
    sort(temp.begin(), temp.end());
    double rowSum = VectorTools::sum(temp);
    if(!(abs(1. - rowSum) < 1e-9)) {
      cout << "Sum of entries for row " << i << " = " << setprecision(20) << rowSum << endl;
    }
  }
  cout << "done." << endl;
  #endif
}  

//composition of transition matrices for different values of rho (but not scaled by rho trans probs as in MMHMM)
VVVdouble MmSmcTransitionProbabilities::fetchCompositeTransitionMatrix(bool missingData) {

  size_t numRhoCateg = mmsmc_ -> getParameterTransitions()[1] -> getNumberOfCategories();
  size_t numIntervals = mmsmc_ -> getNumberOfIntervals();

  //composition of transition matrices
  VVVdouble ctm(numRhoCateg, VVdouble(numIntervals, Vdouble(numIntervals)));
  for(size_t x = 0; x < numRhoCateg; ++x) {
    double rho = mmsmc_ -> getParameterValue("rho") * mmsmc_ -> getParameterScalings()[1] -> getCategories()[x]; 
    //gets (x + 1)th SMC transition matrix:                                                                      
    Vdouble subTreeRowSum(numIntervals);
    for(size_t i = 0; i < numIntervals; ++i) {
      double timeBeta = mmsmc_ -> getTimeIntervals()[i];
      double rowSum = 0.;
      for(size_t j = 0; j < numIntervals; ++j) {
        double timeAlpha = mmsmc_ -> getTimeIntervals()[j];
        if(i == j) { //main diagonal of sub-tree
          ctm[x][i][j] = 0.; //temporary value
        }
        // NOTE this SMC transitions incorporate the demography, which we don't necessarily need in the hierarchical HMM
        // However, we must use this transition probabilities because they take into account self-coalescence events
        else { 
          if(timeBeta < timeAlpha) { 
            ctm[x][i][j] = computeTransitionProbabilityIncreasedTime_(timeBeta, timeAlpha, mmsmc_ -> getLambdaVector(), rho);
          }
          else { 
            ctm[x][i][j] = computeTransitionProbabilityDecreasedTime_(timeBeta, timeAlpha, mmsmc_ -> getLambdaVector(), rho);
          }
          //cout << "prob = " << ctm[x][i][j] << endl;
        }
        rowSum += ctm[x][i][j];
      }
      //cout << "rowSum = " << rowSum << endl;
      subTreeRowSum[i] = rowSum;
    }
    for(size_t i = 0; i < numIntervals; ++i) { 
      ctm[x][i][i] = 1. - subTreeRowSum[i]; //main diagonal
    }
  }

  //if there is missing data in the sequences, we add a "missing" TMRCA
  //that has probability 1. of transitioning to and from all other TMRCAs
  if(missingData) { 
    for(size_t x = 0; x < numRhoCateg; ++x) {
      for(size_t i = 0; i < numIntervals; ++i) { //for every row
        ctm[x][i].push_back(1.); //adds a column (transition into "missing" TMRCA)
      }
      //adds row of "missing" TMRCA transitions into all TMRCAs (including self-transiton)
      vector< double > extraRow(numIntervals + 1, 1.); //+1 for self-transition
      ctm[x].push_back(extraRow);
    }
  }

  #ifdef DEBUG_TM
  for(size_t x = 0; x < numRhoCateg; ++x) {
    for(size_t i = 0; i < ctm[x].size(); ++i) {
      for(size_t j = 0; j < ctm[x].size(); ++j) {
        cout << "R: " << x << " F: " << i << " T: " << j << ": " << setprecision(6) << ctm[x][i][j] << endl;
      }
    }
  }
  #endif

  return ctm;
}

//composition of transition matrices for different values of rho (but not scaled by rho trans probs as in MMHMM)
VVVdouble MmSmcTransitionProbabilities::fetchCompositeTransitionMatrix(bool missingData, size_t numRhoCateg, shared_ptr< bpp::DiscreteDistributionInterface > rhoScaling) {

  size_t numIntervals = mmsmc_ -> getNumberOfIntervals();

  //composition of transition matrices
  VVVdouble ctm(numRhoCateg, VVdouble(numIntervals, Vdouble(numIntervals)));
  for(size_t x = 0; x < numRhoCateg; ++x) {
    double rho = mmsmc_ -> getParameterValue("rho") * rhoScaling -> getCategories()[x]; 
    //gets (x + 1)th SMC transition matrix:                                                                      
    Vdouble subTreeRowSum(numIntervals);
    for(size_t i = 0; i < numIntervals; ++i) {
      double timeBeta = mmsmc_ -> getTimeIntervals()[i];
      double rowSum = 0.;
      for(size_t j = 0; j < numIntervals; ++j) {
        double timeAlpha = mmsmc_ -> getTimeIntervals()[j];
        if(i == j) { //main diagonal of sub-tree
          ctm[x][i][j] = 0.; //temporary value
        }
        // NOTE this SMC transitions incorporate the demography, which we don't necessarily need in the hierarchical HMM
        // However, we must use this transition probabilities because they take into account self-coalescence events
        else { 
          if(timeBeta < timeAlpha) { 
            ctm[x][i][j] = computeTransitionProbabilityIncreasedTime_(timeBeta, timeAlpha, mmsmc_ -> getLambdaVector(), rho);
          }
          else { 
            ctm[x][i][j] = computeTransitionProbabilityDecreasedTime_(timeBeta, timeAlpha, mmsmc_ -> getLambdaVector(), rho);
          }
          //cout << "prob = " << ctm[x][i][j] << endl;
        }
        rowSum += ctm[x][i][j];
      }
      //cout << "rowSum = " << rowSum << endl;
      subTreeRowSum[i] = rowSum;
    }
    for(size_t i = 0; i < numIntervals; ++i) { 
      ctm[x][i][i] = 1. - subTreeRowSum[i]; //main diagonal
    }
  }

  //if there is missing data in the sequences, we add a "missing" TMRCA
  //that has probability 1. of transitioning to and from all other TMRCAs
  if(missingData) { 
    for(size_t x = 0; x < numRhoCateg; ++x) {
      for(size_t i = 0; i < numIntervals; ++i) { //for every row
        ctm[x][i].push_back(1.); //adds a column (transition into "missing" TMRCA)
      }
      //adds row of "missing" TMRCA transitions into all TMRCAs (including self-transiton)
      vector< double > extraRow(numIntervals + 1, 1.); //+1 for self-transition
      ctm[x].push_back(extraRow);
    }
  }

  #ifdef DEBUG_TM
  for(size_t x = 0; x < numRhoCateg; ++x) {
    for(size_t i = 0; i < ctm[x].size(); ++i) {
      for(size_t j = 0; j < ctm[x].size(); ++j) {
        cout << "R: " << x << " F: " << i << " T: " << j << ": " << setprecision(6) << ctm[x][i][j] << endl;
      }
    }
  }
  #endif

  return ctm;
}

void MmSmcTransitionProbabilities::setTransitionProbabilitiesEqualTime_(const Vdouble& rowSumVector, const vector< vector < unsigned char > >& hiddenStates) {
  //finding when nothing changes between hidden states(main diagonal) OR only 1 spatial category changes between hidden states
  for(size_t i = 0; i < mmsmc_ -> getNumberOfHiddenStates(); ++i) {
    double rowComplement = 1. - rowSumVector[i];
    //cout << "rowComplement = " << rowComplement << endl;
    unsigned char timeI = hiddenStates[i][0];
    unsigned char thetaI = hiddenStates[i][1]; 
    unsigned char rhoI = hiddenStates[i][2];
    unsigned char neZeroI = hiddenStates[i][3];  
    for(size_t j = 0; j < mmsmc_ -> getNumberOfHiddenStates(); ++j) {
      unsigned char timeJ = hiddenStates[j][0];
      unsigned char thetaJ = hiddenStates[j][1];
      unsigned char rhoJ = hiddenStates[j][2];
      unsigned char neZeroJ = hiddenStates[j][3];  
      double thetaTransProb = mmsmc_ -> getParameterTransitions()[0] -> getCategoryTransitions()[thetaI][thetaJ];
      double rhoTransProb = mmsmc_ -> getParameterTransitions()[1] -> getCategoryTransitions()[rhoI][rhoJ];  
      double neZeroTransProb = mmsmc_ -> getParameterTransitions()[2] -> getCategoryTransitions()[neZeroI][neZeroJ];
      //size_t numChangedParams = static_cast< size_t >(count_if(begin(stateComparison), end(stateComparison), [](size_t sameState){return sameState == false;}));
      if(timeI == timeJ) {
        //cout << "HS from = " << i << " HS to = " << j << endl;
       // cout << "timeI = " << static_cast< size_t > (timeI) << "; timeJ = " << static_cast< size_t > (timeJ) << endl;
        //cout << "thetaI = " << static_cast< size_t > (thetaI) << "; thetaJ = " << static_cast< size_t > (thetaJ) << endl;
        //cout << "rhoI = " << static_cast< size_t > (rhoI) << "; rhoJ = " << static_cast< size_t > (rhoJ) << endl;
        //cout << "neZeroI = " << static_cast< size_t > (neZeroI) << "; neZeroJ = " << static_cast< size_t > (neZeroJ) << endl;
        double stp = thetaTransProb * rhoTransProb * neZeroTransProb;
        expectedMatrix_[i][j] = rowComplement * stp;
      }
    }
  }
}
  
//takes a particular lambdaVector and rho value and returns a transitionMatrix
VVdouble MmSmcTransitionProbabilities::computeTreeTransitionMatrix_(double rho, const ParameterList& lambdaVector) {
  size_t numberOfTimeIntervals = mmsmc_ -> getNumberOfIntervals();  
  VVdouble treeTransitionMatrix(numberOfTimeIntervals, Vdouble(numberOfTimeIntervals));
  for(size_t i = 0; i < numberOfTimeIntervals; ++i) {
    double timeBeta = mmsmc_ -> getTimeIntervals()[i];
    for(size_t j = 0; j < numberOfTimeIntervals; ++j) {
      double timeAlpha = mmsmc_ -> getTimeIntervals()[j];
      //if we are in the main diagonal, just put a 0 and fill it after everything else is computed:
      if(i == j) { 
        treeTransitionMatrix[i][j] = 0.;
       }
      //else apply the functions defined above:
      else {
        //checks which transition function to use:
        if(timeBeta < timeAlpha) { 
          treeTransitionMatrix[i][j] = computeTransitionProbabilityIncreasedTime_(timeBeta, timeAlpha, lambdaVector, rho);
        }
        else if(timeBeta > timeAlpha) { 
          treeTransitionMatrix[i][j] = computeTransitionProbabilityDecreasedTime_(timeBeta, timeAlpha, lambdaVector, rho);
        }
      }
    }
  }
  //now we fill the main diagonal of the matrix:
  for(size_t k = 0; k < numberOfTimeIntervals; ++k) {
    //the sum of all the elements in the row:
    double rowSum = 0.;
    for(size_t l = 0; l < numberOfTimeIntervals; ++l) {
      rowSum += treeTransitionMatrix[k][l];
    }
    treeTransitionMatrix[k][k] = 1. - rowSum;
  }
  return treeTransitionMatrix;
}  
  
vector< bool > MmSmcTransitionProbabilities::compareHiddenStates_(const vector< unsigned char >& hs1, const vector< unsigned char >& hs2) {
  vector< bool > stateComparison = { true, true, true, true };
  for(size_t i = 0; i < hs1.size(); ++i) {
    if(hs1[i] != hs2[i]) { 
      stateComparison[i] = false;
    }
  }
  return stateComparison;
}
