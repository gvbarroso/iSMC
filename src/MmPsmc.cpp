/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 09/09/2019
 *
 */


#include <string>
#include <vector>
#include <iostream>

#ifdef SIMPLEZIPHMM
  #include <SimpleZipHMM/forwarder.hpp>
  #include <SimpleZipHMM/posterior_decoding.hpp>
  #include <SimpleZipHMM/hmm_io.hpp>
#else
  #include <zipHMM/forwarder.hpp>
  #include <zipHMM/posterior_decoding.hpp>
  #include <zipHMM/hmm_io.hpp>
#endif //SIMPLEZIPHMM

#include "MmPsmc.h"

//#define EXTRA_DECODING //to decode & output landscapes other than the posterior avg. rate

using namespace std;
using namespace bpp;


vector< size_t > MmPsmc::fetchTransitionSequence(vector< size_t >& treeSequence, size_t start, size_t end, bool missingData) { //for rho
  
  size_t numIntervals = mmsmc_ -> getNumberOfIntervals();
  size_t lengthTmrcaDist = treeSequence.size();
  vector< size_t > transitionSequence(lengthTmrcaDist); //ret 
    
  if(missingData) {
    maskTreeSequence_(treeSequence, start, end);
    ++numIntervals; //adds the "missing" TMRCA
  }

  for(size_t i = 0; i < lengthTmrcaDist - 1; ++i) { //for all sites except the last
    //index of the observation in the HmmLib vector< pair< size_t, size_t > >
    transitionSequence[i] = treeSequence[i] * numIntervals + treeSequence[i + 1];
    //cout << "ts: " << i << "; from: " << treeSequence[i] << "; to: " << treeSequence[i + 1] << "; code: " << transitionSequence[i] << endl;
  }
  
  //duplicates last tree transition to avoid rec. landscapes being short of one site per fragment
  transitionSequence[lengthTmrcaDist - 1] = transitionSequence[lengthTmrcaDist - 2];
    
  return transitionSequence;
}

vector< size_t > MmPsmc::fetchEmissionSequence(const vector< size_t >& treeSequence, size_t start, size_t end) { //for theta
    
  vector< unsigned char > seq = fetchFragment(start, end); //gets SNP sequence
  
  vector< size_t > emissionSequence(treeSequence.size()); //ret 
  
  size_t numObs = mmep_ -> getNumberOfObservedStates(); //num obs in the original HMM
  
  for(size_t i = 0; i < treeSequence.size(); ++i) { //for all sites
      
    //index of the observation in the HmmLib vector< pair< size_t, size_t > >
    emissionSequence[i] = treeSequence[i] * numObs + seq[i]; //NOTE has built-in masking of layers?
    //cout << "t: " << treeSequence[i] << "; o: " << seq[i] << "; idx: " << emissionSequence[i] << endl;
  }
  
  return emissionSequence;
}

vector< unsigned char > MmPsmc::fetchLocalTreesIndices() {
    
  vector< unsigned char > localTrees(reconstructedHiddenStates_.size());
  
  for(size_t i = 0; i < reconstructedHiddenStates_.size(); ++i) {
      
    unsigned int state = reconstructedHiddenStates_[i];
    localTrees[i] = mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[state][0];
  }
  
  return localTrees;
}

vector< size_t > MmPsmc::fetchLocalTrees(const zipHMM::Matrix& postProbMatrix) { 
    
  vector< size_t > treeSequence(postProbMatrix.get_width()); //ret 
  vector< double >::iterator maxProb;
  
  VVdouble postProbTrees = getTreePosteriorMatrix_(postProbMatrix);
  
  for(size_t i = 0; i < postProbTrees.size(); ++i) { //for all sites
      
    maxProb = max_element(postProbTrees[i].begin(), postProbTrees[i].end());
    treeSequence[i] = distance(postProbTrees[i].begin(), maxProb);
  }
  return treeSequence;
}

//this polymorphic version fills in treeProbs to make computeAverageRateWithinDecodedTrees_ faster
vector< size_t > MmPsmc::fetchLocalTrees(const zipHMM::Matrix& postProbMatrix, Vdouble& treeProbs) { 
    
  vector< size_t > treeSequence(postProbMatrix.get_width()); //ret 
  vector< double >::iterator maxProb;
  
  VVdouble postProbTrees = getTreePosteriorMatrix_(postProbMatrix);
  
  for(size_t i = 0; i < postProbTrees.size(); ++i) { //for all sites
      
    maxProb = max_element(postProbTrees[i].begin(), postProbTrees[i].end());
    treeProbs[i] = *maxProb;
    treeSequence[i] = distance(postProbTrees[i].begin(), maxProb);
  }
  return treeSequence;
}

void MmPsmc::decodeLocalTmrca(const zipHMM::Matrix& postProbMatrix, const string& fileName) { 
    
  ofstream decFile;
  decFile.open(fileName + "_decoded_TMRCA.txt.gz",
               std::ios_base::out | std::ios_base::binary);
  
  boost::iostreams::filtering_ostream tmrcaStream; 
  
  tmrcaStream.push(boost::iostreams::gzip_compressor());
  tmrcaStream.push(decFile);
    
  vector< double >::iterator maxProb;
  
  VVdouble postProbTrees = getTreePosteriorMatrix_(postProbMatrix);
  
  for(size_t i = 0; i < postProbTrees.size(); ++i) { //for all sites
      
    maxProb = max_element(postProbTrees[i].begin(), postProbTrees[i].end());
    
    size_t treeIndex = distance(postProbTrees[i].begin(), maxProb);
    
    tmrcaStream << mmsmc_ -> getTimeIntervals()[treeIndex] << endl;
  }
  
  boost::iostreams::close(tmrcaStream);
}

vector< unsigned char > MmPsmc::fetchLocalRateIndices(const string& rate) {
    
  size_t index = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second; 
  
  vector< unsigned char > localRates(reconstructedHiddenStates_.size());
  
  for(size_t i = 0; i < reconstructedHiddenStates_.size(); ++i) {
      
    unsigned char state = reconstructedHiddenStates_[i];
    
    localRates[i] = mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[state][index];
  }
  
  return localRates;
}

vector< size_t > MmPsmc::fetchLocalRateIndices(const zipHMM::Matrix& postProbMatrix, const string& rate) {
    
  vector< size_t > rateSequence(postProbMatrix.get_width()); //ret  
  vector< double >::iterator maxProb;
  
  VVdouble postProbRate = getRatePosteriorMatrix_(postProbMatrix, rate);
  
  for(size_t i = 0; i < postProbRate.size(); ++i) {
    maxProb = max_element(postProbRate[i].begin(), postProbRate[i].end());
    rateSequence[i] = distance(postProbRate[i].begin(), maxProb);
  }
  
  return rateSequence;
}

void MmPsmc::computeLocalAverageRate(const string& rate) { //OLD version

  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  
  posteriorAverageByPosition_.resize(selectedFragment_.size());
  
  //computes posterior average based on model of rate heterogeneity
  if(mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getHeterogeneousRateModel() == "Gamma") { 
    computeGammaAverageRate_(rate);
  }
  
  else if(mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getHeterogeneousRateModel() == "Gamma+Hotspot") { 
    computeGammaWithHotspotAverageRate_(rate);     
  }
  
  else if(mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getHeterogeneousRateModel() == "Hotspot") { 
    computeHotspotAverageRate_(rate);
  }
}

//wraps up the different prior distributions on rate values
void MmPsmc::computeLocalAverageRate(const zipHMM::Matrix& postProbMatrix, const string& fileName) {

  for(size_t i = 0; i < mmsmc_ -> getParameterTransitions().size(); ++i) {
    
    if(mmsmc_ -> getParameterTransitions()[i] -> getNumberOfCategories() > 1) {
    
      string rate = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().left.find(i) -> second;
  
      string newName = fileName + "_estimated_" + rate + ".txt.gz";
  
      ofstream rateFile; //write it directly to file to save memory
      rateFile.open(newName, std::ios_base::out | std::ios_base::binary);
  
      boost::iostreams::filtering_ostream landStream; 
  
      landStream.push(boost::iostreams::gzip_compressor());
      landStream.push(rateFile);
  
      if(mmsmc_ -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Gamma") { 
        computeGammaAverageRate_(postProbMatrix, rate, landStream); 
      }
  
      else if(mmsmc_ -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Gamma+Hotspot") { 
        computeGammaWithHotspotAverageRate_(postProbMatrix, rate, landStream);
      }
  
      else if(mmsmc_ -> getParameterTransitions()[i] -> getHeterogeneousRateModel() == "Hotspot") { 
        computeHotspotAverageRate_(postProbMatrix, rate, landStream);
      }
      
      boost::iostreams::close(landStream);
    }
  }
}

void MmPsmc::decodeLocalRate(const VVdouble& ratePostProbs, const string& fileName, const string& rate) {

  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  string newName = fileName + "_decoded_" + rate + ".txt.gz";
  
  ofstream decFile;
  decFile.open(newName, std::ios_base::out | std::ios_base::binary);
  
  boost::iostreams::filtering_ostream landStream; 
  
  landStream.push(boost::iostreams::gzip_compressor());
  landStream.push(decFile);
  
  if(mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getHeterogeneousRateModel() == "Gamma") { 
    computeGammaRateMaxProb_(ratePostProbs, rate, landStream);
  }
  
  else if(mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getHeterogeneousRateModel() == "Hotspot") { 
    computeHotspotRateMaxProb_(ratePostProbs, rate, landStream);
  }
  
  else if(mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getHeterogeneousRateModel() == "Gamma+Hotspot") { 
    computeGammaWithHotspotRateMaxProb_(ratePostProbs, rate, landStream);
  }
  
  boost::iostreams::close(landStream);
}

void MmPsmc::decodeLocalRate(const zipHMM::Matrix& ratePostProbs, const string& fileName, const string& rate) {

  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  
  string newName = fileName + "_decoded_" + rate + ".txt.gz";
  
  ofstream decFile;
  decFile.open(newName, std::ios_base::out | std::ios_base::binary);
  
  boost::iostreams::filtering_ostream landStream; 
  
  landStream.push(boost::iostreams::gzip_compressor());
  landStream.push(decFile);
  
  if(mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getHeterogeneousRateModel() == "Gamma") { 
    computeGammaRateMaxProb_(ratePostProbs, rate, landStream);
  }
  
  else if(mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getHeterogeneousRateModel() == "Hotspot") { 
    computeHotspotRateMaxProb_(ratePostProbs, rate, landStream);
  }
  
  else if(mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getHeterogeneousRateModel() == "Gamma+Hotspot") { 
    computeGammaWithHotspotRateMaxProb_(ratePostProbs, rate, landStream);
  }
  
  boost::iostreams::close(landStream);
}

void MmPsmc::computeLocalAverageTmrca() {

  posteriorAverageByPosition_.resize(selectedFragment_.size());

  for(size_t i = 0; i < selectedFragment_.size(); ++i) {
      
    double avgTmrca = 0.;
    
    for(size_t j = 0; j < mmsmc_ -> getNumberOfIntervals(); ++j) {
        
      double focalTreeProb = 0.;
      
      //we search hidden states for a match between their trees and our tree j 
      for(size_t k = 0; k < mmsmc_ -> getNumberOfHiddenStates(); ++k) {
          
        if(mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][0] == j) {
          //when we find a tree match, we increment the focal probability
          focalTreeProb += posteriorProbMatrix_[i][k];
        }
      }
      
      //we weigh the tree height by its probability, and increment
      avgTmrca += mmsmc_ -> getAverageCoalescenceTimeVector()[j] * focalTreeProb;
    }
    
    posteriorAverageByPosition_[i] = avgTmrca;
  }
}

void MmPsmc::computeLocalAverageTmrca(const zipHMM::Matrix& postProbMatrix, const string& fileName) {

  ofstream tmrcaFile; 
  tmrcaFile.open(fileName + "_estimated_TMRCA.txt.gz", 
                 std::ios_base::out | std::ios_base::binary);

  boost::iostreams::filtering_ostream tmrcaStream; 
  
  tmrcaStream.push(boost::iostreams::gzip_compressor());
  tmrcaStream.push(tmrcaFile);
    
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) {
      
    double avgTmrca = 0.;
    
    for(size_t j = 0; j < mmsmc_ -> getNumberOfIntervals(); ++j) {
        
      double focalTreeProb = 0.;
      
      //we search hidden states for a match between their trees and our tree j 
      for(size_t k = 0; k < mmsmc_ -> getNumberOfHiddenStates(); ++k) {
          
        if(mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][0] == j) {
          //when we find a tree match, we increment the focal probability
          focalTreeProb += postProbMatrix(k, i);
        }
      }
      
      //we weight the tree height by its probability, and increment
      double focalTreeHeight = mmsmc_ -> getAverageCoalescenceTimeVector()[j];
      
      avgTmrca += focalTreeHeight * focalTreeProb;
    }
    
    tmrcaStream << avgTmrca << endl;
  }
  
  boost::iostreams::close(tmrcaStream);
}

//memory-efficient version of the function, directly decodes and writes to file
void MmPsmc::posteriorDecodingUsingZipHMM(size_t genomicStart, size_t genomicEnd, const string& fileName,
                                          bool restricted, zipHMM::Matrix& postProbMatrix, bpp::Vdouble& pi) {
  
  //copy fragment (otherwise there is an obvious problem when decoding fragments in parallel)
  vector< unsigned char > copyOfFragment = fetchFragment(genomicStart, genomicEnd);
  
  zipHMM::posterior_decoding(copyOfFragment,
                             pi, //initialisation probabilities updated per fragment
                             mmtp_ -> getExpectedMatrix(),
                             mmep_ -> getExpectedMatrix(),
                             postProbMatrix);
  
  if(restricted) {
    printDecodedStates(postProbMatrix, fileName); //time-restricted modulation
  }
  
  else {
    computeLocalAverageRate(postProbMatrix, fileName);
  }

  #ifdef EXTRA_DECODING
  computeLocalAverageTmrca(postProbMatrix, fileName);
  
  string shannonTotal = fileName + "_total";
  computePerSiteShannonEquitability_(postProbMatrix, shannonTotal);

  string shannonRate = fileName + "_" + rate;
  computePerSiteShannonEquitability_(getRatePosteriorMatrix_(postProbMatrix, rate), shannonRate);
  
  sting shannonTree = fileName + "_tree";
  computePerSiteShannonEquitability_(getTreePosteriorMatrix_(postProbMatrix), shannonTree);
  #endif

} 

void MmPsmc::computePosteriorProbs(size_t genomicStart, size_t genomicEnd,
                                   zipHMM::Matrix& postProbMatrix, Vdouble& pi) {
  
  //copy fragment (otherwise there is an obvious problem when decoding fragments in parallel)
  vector< unsigned char > copyOfFragment = fetchFragment(genomicStart, genomicEnd);
  zipHMM::posterior_decoding(copyOfFragment,
                             pi, //initialisation probabilities updated per fragment
                             mmtp_ -> getExpectedMatrix(),
                             mmep_ -> getExpectedMatrix(),
                             postProbMatrix);
}

//gets posterior probability matrix considering only rate categories as hidden states
VVdouble MmPsmc::getRatePosteriorMatrix(const VVdouble& postProbMatrix, const string& rate) { 
    
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t rateIndexInLibrary = rateIndexInVectors + 1;
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();
  size_t numCategories = mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getNumberOfCategories();

  VVdouble ratePostProbs(postProbMatrix.size(), Vdouble(numCategories)); 
  
  for(size_t i = 0; i < postProbMatrix.size(); ++i) { //for all sites
      
    for(size_t j = 0; j < numCategories; ++j) {
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
          
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) { 
          focalCatProb += postProbMatrix[i][k];
        }
      }
      
      ratePostProbs[i][j] = focalCatProb;
    } 
  }
  
  return ratePostProbs;
}

//gets posterior probability matrix considering only rate categories as hidden states
VVdouble MmPsmc::getRatePosteriorMatrix(const zipHMM::Matrix& postProbMatrix, const string& rate) { 
    
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t rateIndexInLibrary = rateIndexInVectors + 1;
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();
  size_t numCategories = mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getNumberOfCategories();
  
  VVdouble ratePostProbs(postProbMatrix.get_width(), Vdouble(numCategories)); 
  
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) { //for all sites
      
    for(size_t j = 0; j < numCategories; ++j) {
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
          
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) { 
          focalCatProb += postProbMatrix(k, i);  
        }
      }
      
      ratePostProbs[i][j] = focalCatProb;
    }
  }
  
  return ratePostProbs; 
}

//prints a sheet with tree posterior probs per position
void MmPsmc::printTmrcaPosteriorMatrix(const zipHMM::Matrix& postProbMatrix, const string& fileName) {
 
  size_t numIntervals = mmsmc_ -> getNumberOfIntervals();
  
  string newName = fileName + "_tmrca_posterior_probs.txt.gz";
  
  ofstream tmrcaFile; 
  tmrcaFile.open(newName, std::ios_base::out | std::ios_base::binary);
  
  boost::iostreams::filtering_ostream tmrcaStream; 
  
  tmrcaStream.push(boost::iostreams::gzip_compressor());
  tmrcaStream.push(tmrcaFile);
  
  VVdouble tmrcaPostProbs = getTreePosteriorMatrix_(postProbMatrix); //site -> time interval
  
  for(size_t i = 0; i < tmrcaPostProbs.size(); ++i) { //for all sites
      
    for(size_t j = 0; j < numIntervals; ++j) {
        
      tmrcaStream << tmrcaPostProbs[i][j];
      
      if(j < numIntervals - 1) {
        tmrcaStream << "\t";  
      }
    }
    
    tmrcaStream << endl;
  }
  
  boost::iostreams::close(tmrcaStream);
}

//prints a sheet with rates per position
void MmPsmc::printRatePerInterval(const zipHMM::Matrix& postProbMatrix, const string& rate, const string& fileName) {
 
  size_t numIntervals = mmsmc_ -> getNumberOfIntervals();
  
  VVdouble ratePerInterval = computeAverageRatePerInterval_(postProbMatrix, rate);
  
  string newName = fileName + "_rate_per_interval.txt.gz";
  
  ofstream rateFile; 
  rateFile.open(newName, std::ios_base::out | std::ios_base::binary);
  
  boost::iostreams::filtering_ostream rateStream; 
  
  rateStream.push(boost::iostreams::gzip_compressor());
  rateStream.push(rateFile);
    
  for(size_t i = 0; i < ratePerInterval.size(); ++i) { //for all sites
      
    for(size_t j = 0; j < numIntervals; ++j) {
        
      rateStream << ratePerInterval[i][j];
      
      if(j < numIntervals - 1) {
        rateStream << "\t";  
      }
    }
    rateStream << endl;
  }
  boost::iostreams::close(rateStream);
}

//TESTING
void MmPsmc::printDecodedStates(const zipHMM::Matrix& postProbMatrix, const string& fileName) {
  
  //first we decode the tmrca and fill maxTreeProbs
  Vdouble maxTreeProbs(postProbMatrix.get_width());
  vector< size_t > decodedTmrca = fetchLocalTrees(postProbMatrix, maxTreeProbs); //max prob tree for every site
  
  VVdouble ratesPerInterval(0, Vdouble(0)); //
  
  string rateNames = ""; //
  
  //figures out which rates to decode
  for(size_t j = 0; j < mmsmc_ -> getParameterTransitions().size(); ++j) {
      
    if(mmsmc_ -> getParameterTransitions()[j] -> getNumberOfCategories() > 1) {
      //now we get the average rate per time interval  
      string focalRate = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().left.at(j);
      ratesPerInterval.push_back(computeAverageRateWithinDecodedTrees_(postProbMatrix, decodedTmrca, maxTreeProbs, focalRate));
      
      rateNames += (focalRate + ".");
    }
  }
  
  string newName = fileName + "_Tmrca." + rateNames + "txt.gz";
  
  ofstream rateFile; //write it directly to file to save memory
  rateFile.open(newName, std::ios_base::out | std::ios_base::binary);
  
  boost::iostreams::filtering_ostream landStream; 
  
  landStream.push(boost::iostreams::gzip_compressor());
  landStream.push(rateFile);
  
  //prints decoded TMRCA (max prob) as well as average rates within such interval   
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) { //for all sites
    
    landStream << decodedTmrca[i] << "\t";
    
    for(size_t j = 0; j < ratesPerInterval.size(); ++j) { //for heterogeneous rates  
      landStream << ratesPerInterval[j][i] << "\t";
    }
    
    landStream << endl;
  }

  boost::iostreams::close(landStream);
}

void MmPsmc::computeGammaRateMaxProb_(const VVdouble& ratePostProbs, const string& rate, boost::iostreams::filtering_ostream& file) { 
    
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t numCategories = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();
  double rateVal = mmsmc_ -> getParameterValue(rate);
  
  for(size_t i = 0; i < ratePostProbs.size(); ++i) { 
      
    size_t stateIndex = 0;
    double highestProb = ratePostProbs[i][0];
    
    for(size_t j = 1; j < numCategories; ++j) {
        
      if(ratePostProbs[i][j] > highestProb) {
          
        highestProb = ratePostProbs[i][j];
        stateIndex = j;
      }
    }
    
    if(rate != "ne") {
      file << rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex] << endl;
    }
    
    else {
      file << mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex] << endl;
    }
  }
}

void MmPsmc::computeGammaRateMaxProb_(const zipHMM::Matrix& ratePostProbs, const string& rate, boost::iostreams::filtering_ostream& file) { 
    
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t numCategories = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();
  double rateVal = mmsmc_ -> getParameterValue(rate);
  
  for(size_t i = 0; i < ratePostProbs.get_width(); ++i) { 
      
    size_t stateIndex = 0;
    double highestProb = ratePostProbs(0, i);
    
    for(size_t j = 1; j < numCategories; ++j) {
        
      if(ratePostProbs(j, i) > highestProb) {
        highestProb = ratePostProbs(j, i);
        stateIndex = j;
      }
    }
    
    if(rate != "ne") {
      file << rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex] << endl;
    }
    
    else {
      file << mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex] << endl;
    }
  }
}

void MmPsmc::computeHotspotRateMaxProb_(const VVdouble& ratePostProbs, const string& rate,
                                        boost::iostreams::filtering_ostream& file) { 
    
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  double rateVal = mmsmc_ -> getParameterValue(rate);
  
  //Hotspot model has only 2 categories, choosing state with max. prob. is straightforward
  for(size_t i = 0; i < ratePostProbs.size(); ++i) { 
    size_t stateIndex = 0;
    if(ratePostProbs[i][1] > ratePostProbs[i][0]) {
      stateIndex = 1;
    }
    if(rate != "ne") {
      file << rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex] << endl;
    }
    else {
      file << mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex] << endl;
    }
  }
}

void MmPsmc::computeHotspotRateMaxProb_(const zipHMM::Matrix& postProbMatrix, const string& rate,
                                        boost::iostreams::filtering_ostream& file) { 
    
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  double rateVal = mmsmc_ -> getParameterValue(rate);

  //Hotspot model has only 2 categories, choosing state with max. prob. is straightforward
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) { 
      
    size_t stateIndex = 0;
    
    if(postProbMatrix(1, i) > postProbMatrix(0, i)) {
      stateIndex = 1;
    }
    
    if(rate != "ne") {
      file << rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex] << endl;
    }
    
    else {
      file << mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex] << endl;
    }
    
  }
}

void MmPsmc::computeGammaWithHotspotRateMaxProb_(const VVdouble& ratePostProbs, const string& rate,
                                                 boost::iostreams::filtering_ostream& file) { 
    
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t numCategories = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();
  double rateVal = mmsmc_ -> getParameterValue(rate);

  for(size_t i = 0; i < ratePostProbs.size(); ++i) { 
      
    size_t stateIndex = 0;
    double highestProb = ratePostProbs[i][0];
    
    for(size_t j = 1; j < numCategories; ++j) {
        
      if(ratePostProbs[i][j] > highestProb) {
          
        highestProb = ratePostProbs[i][j];
        stateIndex = j;
      }
    }
    
    if(rate != "ne") {
        
      if(stateIndex == 0) { //if hotspot 
        file << rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getParameterValue("heat") << endl;
      }
      
      else {
        file << mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex - 1] << endl;
      }  
    }
    
    else {
      if(stateIndex == 0) { //if hotspot
        file << mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getParameterValue("heat") << endl;
      }
      else {
        file << mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex - 1] << endl; 
      }
    }
  }
}

void MmPsmc::computeGammaWithHotspotRateMaxProb_(const zipHMM::Matrix& ratePostProbs, const string& rate,
                                                 boost::iostreams::filtering_ostream& file) { 

  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t numCategories = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();

  for(size_t i = 0; i < ratePostProbs.get_width(); ++i) { 
    
    size_t stateIndex = 0;
    double highestProb = ratePostProbs(0, i);
    
    for(size_t j = 1; j < numCategories; ++j) {

      if(ratePostProbs(j, i) > highestProb) {
        highestProb = ratePostProbs(j, i);
        stateIndex = j;
      }
    }
    
    if(rate != "ne") {
        
      if(stateIndex == 0) { //if hotspot 
        file << mmsmc_ -> getParameterValue(rate) * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getParameterValue("heat") << endl;
      }
      
      else {
        file << mmsmc_ -> getParameterValue(rate) * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex - 1] << endl;
      }  
    }
    
    else {
        
      if(stateIndex == 0) { //if hotspot
        file << mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getParameterValue("heat") << endl;
      }
      
      else {
        file << mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[stateIndex - 1] << endl; 
      }
    }
  }
}

void MmPsmc::computeHotspotAverageRate_(const string& rate) {

  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;  
  size_t rateIndexInLibrary = rateIndexInVectors + 1;
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();  
  double rateVal = mmsmc_ -> getParameterValue(rate);

  for(size_t i = 0; i < selectedFragment_.size(); ++i) {
      
    double avgRate = 0.;
    
    for(size_t j = 0; j < 2; ++j) { //only two states
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
          
        //we search for a match between the focal category index and the and the parameter index
        //in hidden state k (query rateIndexInLibrary to match the structure of HiddenStatesLibrary):
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) { 
          focalCatProb += posteriorProbMatrix_[i][k]; //increment the total probability of this category:
        }
      }
      
      if(rate != "ne") { //if the rate we are interested in is either theta or rho
        //the reconstructed rate at this position is the genome-average rate * the focal scaling factor, weighted by focalCatProb  
        avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j] * focalCatProb;
      }
      
      else { //if it is ne (spatial, local, which is not an explicit parameter of the model)
        avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j] * focalCatProb;
      }
    }
    
    posteriorAverageByPosition_[i] = avgRate;
  }      
}

void MmPsmc::computeHotspotAverageRate_(const zipHMM::Matrix& postProbMatrix, const string& rate,
                                        boost::iostreams::filtering_ostream& file) {
  
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();  
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t rateIndexInLibrary = rateIndexInVectors + 1;
  double rateVal = mmsmc_ -> getParameterValue(rate);
  
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) {
      
    double avgRate = 0.;
    
    for(size_t j = 0; j < 2; ++j) { //only two states
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
          
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) { 
          focalCatProb += postProbMatrix(k, i);  
        }
      }
      
      if(rate != "ne") {
        avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j] * focalCatProb;
      }
      
      else { //if it is ne (spatial, local, which is not an explicit parameter of the model)
        avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j] * focalCatProb;
      }
    }
    
    file << avgRate << endl;
  }
}

void MmPsmc::computeGammaAverageRate_(const string& rate) {
  
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;  
  size_t rateIndexInLibrary = rateIndexInVectors + 1;
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();
  size_t numCateg = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();
  double rateVal = mmsmc_ -> getParameterValue(rate);
  
  for(size_t i = 0; i < selectedFragment_.size(); ++i) {
      
    double avgRate = 0.;
    
    for(size_t j = 0; j < numCateg; ++j) {
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
        //we search for a match between the focal category index and the and the parameter index
        //in hidden state k (query rateIndexInLibrary to match the structure of HiddenStatesLibrary):
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) { 
          //when we find it, we increment the total probability of this category:
          focalCatProb += posteriorProbMatrix_[i][k];
        }
      }
      
      if(rate != "ne") { //if the rate we are interested in is either theta or rho
        //the reconstructed rate at this position is the genome-average rate * the focal scaling factor, weighted by focalCatProb
        avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j] * focalCatProb;
      }
      
      else { //if it is ne (spatial, local, which is not an explicit parameter of the model)
        avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j] * focalCatProb;
      }
    }
    
    posteriorAverageByPosition_[i] = avgRate;
  }  
}

void MmPsmc::computeGammaAverageRate_(const zipHMM::Matrix& postProbMatrix, const string& rate,
                                      boost::iostreams::filtering_ostream& file) { 
  
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t rateIndexInLibrary = rateIndexInVectors + 1;
  size_t numCateg = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();
  double rateVal = mmsmc_ -> getParameterValue(rate);
  
  //cout << "rate: " << rate << "; lib. index: " << rateIndexInLibrary << "; vec. index: " << rateIndexInVectors << "; # categ.: " << numCateg << endl;
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) { //for all sites
      
    double avgRate = 0.;  
    
    for(size_t j = 0; j < numCateg; ++j) {
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
        //we search for a match between the focal category index and the and the parameter index
        //in hidden state k (query rateIndexInLibrary to match the structure of HiddenStatesLibrary):
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) { 
          //when we find it, we increment the total probability of this category:
          //cout << "increment prob: j = " << j << "; k = " << k << endl;
          focalCatProb += postProbMatrix(k, i);  
        }
      }
      
      if(rate != "ne") {
        avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j] * focalCatProb;
      }
      
      else {
        avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j] * focalCatProb;
      }
    }
    
    file << avgRate << endl;
  }
}


//TESTING
VVdouble MmPsmc::computeAverageRatePerInterval_(const zipHMM::Matrix& postProbMatrix, const string& rate) {
  
  size_t numIntervals = mmsmc_ -> getNumberOfIntervals();
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t numCateg = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();
  
  VVdouble ratePerInterval(postProbMatrix.get_width(), Vdouble(numIntervals)); //site -> time interval 
      
  double rateVal = mmsmc_ -> getParameterValue(rate); //genome-wide average rate

  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) { //for all sites
      
    for(size_t j = 0; j < numIntervals; ++j) {
      
      double avgRate = 0.; //weighted average rate (over all categories) of rate inside interval j
      double sumPostProb = 0.; //sum of rate posterior probabilities (over all categories), inside interval j -> to compute the conditional prob of being in j;
      
      for(size_t k = 0; k < numCateg; ++k) { //first run to compute sum of post prob of rate inside j
        sumPostProb += postProbMatrix(k + (j * numCateg), i); 
      }
      
      for(size_t k = 0; k < numCateg; ++k) {
          
        size_t hsIndex = k + (j * numCateg); //hidden state matching interval j and category k
        double catPostProb = postProbMatrix(hsIndex, i); //posterior probability for category k inside interval j
        
        if(rate != "ne") {
          avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[k] * catPostProb / sumPostProb;
        }
      
        else {
          avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[k] * catPostProb / sumPostProb;
        }
      }
      ratePerInterval[i][j] = avgRate;
    }
  }
  return ratePerInterval;
}

//TESTING this is a more straight-to-the-point version of computeAverageRatePerInterval_ that only cares about the tree with max prob
Vdouble MmPsmc::computeAverageRateWithinDecodedTrees_(const zipHMM::Matrix& postProbMatrix, const vector< size_t >& decodedTmrca,
                                                      const Vdouble& maxTreeProbs, const string& rate) {
  
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t numCateg = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();
  
  Vdouble ratePerInterval(postProbMatrix.get_width()); //one value per site
      
  double rateVal = mmsmc_ -> getParameterValue(rate); //genome-wide average rate

  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) { //for all sites
      
    size_t tree = decodedTmrca[i]; //index of tree with max prob. at position i
    double treePostProb = maxTreeProbs[i];  
    double avgRate = 0.; //weighted average rate (over all categories) inside tree
       
    for(size_t k = 0; k < numCateg; ++k) {
          
      size_t hsIndex = k + (tree * numCateg); //hidden state matching interval j and category k
      double catPostProb = postProbMatrix(hsIndex, i); //posterior probability for category k inside interval j
        
      if(rate != "ne") {
        avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[k] * catPostProb / treePostProb;
      }
     
      else {
        avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[k] * catPostProb / treePostProb;
      }
    }
    ratePerInterval[i] = avgRate;
  }
  return ratePerInterval;
}

void MmPsmc::computeGammaWithHotspotAverageRate_(const string& rate) {
    
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;  
  size_t rateIndexInLibrary = rateIndexInVectors + 1;  
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();
  size_t numCateg = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();
  double rateVal = mmsmc_ -> getParameterValue(rate);
  
  for(size_t i = 0; i < selectedFragment_.size(); ++i) {
      
    double avgRate = 0.;
    
    for(size_t j = 0; j < numCateg; ++j) {
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
          
        //we search for a match between the focal category index and the and the parameter index
        //in hidden state k (query rateIndexInLibrary to match the structure of HiddenStatesLibrary):
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) { 
          //when we find it, we increment the total probability of this category:
          focalCatProb += posteriorProbMatrix_[i][k];
        }
      }
      
      if(rate != "ne") { //if the rate we are interested in is either theta or rho
        //the reconstructed rate at this position is the genome-average rate * the focal scaling factor, weighted by focalCatProb  
          
        if(j == 0) { //if hotspot
          avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getParameterValue("heat") * focalCatProb;
        }
        
        else {
          avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j - 1] * focalCatProb;
        }
      }
      
      else { //if it is ne (spatial, local, which is not an explicit parameter of the model)
          
        if(j == 0) { //if hotspot
          avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getParameterValue("heat") * focalCatProb;  
        }
        
        else {
          avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j - 1] * focalCatProb;
        }
      }
    }
    
    posteriorAverageByPosition_[i] = avgRate;
  }   
}

void MmPsmc::computeGammaWithHotspotAverageRate_(const zipHMM::Matrix& postProbMatrix, const string& rate,
                                                 boost::iostreams::filtering_ostream& file) {
    
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();
  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t rateIndexInLibrary = rateIndexInVectors + 1;
  size_t numCateg = mmsmc_ -> getParameterTransitions()[rateIndexInVectors] -> getNumberOfCategories();
  double rateVal = mmsmc_ -> getParameterValue(rate);
  
  //cout << "rate: " << rate << "; vec. index: " << rateIndexInVectors << "; Gamma+Hotspot w/ " << numCateg << " categ." << endl;
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) { 
      
    double avgRate = 0.;  
    
    for(size_t j = 0; j < numCateg; ++j) {
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
          
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) {  
          //cout << "increment prob: j = " << j << "; k = " << k << "; vectorIndex = " << vectorIndex << endl;    
          focalCatProb += postProbMatrix(k, i);
        }
      }
      
      if(rate != "ne") {
          
        if(j == 0) { //if hotspot
          avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getParameterValue("heat") * focalCatProb; 
        }
        
        else {
          avgRate += rateVal * mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j - 1] * focalCatProb; 
        }
      }
      
      else {
          
        if(j == 0) { //if hotspot
          avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getParameterValue("heat") * focalCatProb;
        }
        
        else {
          avgRate += mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getCategories()[j - 1] * focalCatProb; 
        }
        
      }
    }
    
    file << avgRate << endl;
  }  
}

//gets posterior probability matrix considering only trees as hidden states
VVdouble MmPsmc::getTreePosteriorMatrix_(const VVdouble& postProbMatrix) { 

  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();  
  size_t numIntervals = mmsmc_ -> getNumberOfIntervals();
  
  //using the jump factor is easier for Trees than for rates 
  //this is because trees are always at pos. 0 of HS config.
  //while rates depend on the number of modulated params.
  size_t jumpFactor = numHiddenStates / numIntervals; //multiples by definition

  VVdouble treePostProbs(postProbMatrix.size(), Vdouble(numIntervals)); 
  
  for(size_t i = 0; i < postProbMatrix.size(); ++i) { //for all sites
      
    for(size_t j = 0; j < numIntervals; ++j) {
        
      double focalTreeProb = 0.;
      
      for(size_t k = 0; k < jumpFactor; ++k) { //jumps HS to get the correct trees
          
        size_t hsIndex = k + j * jumpFactor;
        
        focalTreeProb += postProbMatrix[i][hsIndex];
      }
      
      treePostProbs[i][j] = focalTreeProb;
    }
  }
  
  return treePostProbs; 
}

//gets posterior probability matrix considering only trees as hidden states
VVdouble MmPsmc::getTreePosteriorMatrix_(const zipHMM::Matrix& postProbMatrix) { 

  size_t numIntervals = mmsmc_ -> getNumberOfIntervals();
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();

  //using the jump factor is easier for Trees than for rates 
  //this is because trees are always at pos. 0 of HS config.
  //while rates depend on the number of modulated params.
  size_t jumpFactor = numHiddenStates / numIntervals; //multiples by definition
  
  VVdouble treePostProbs(postProbMatrix.get_width(), Vdouble(numIntervals)); 
  
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) { //for all sites
      
    for(size_t j = 0; j < numIntervals; ++j) {
        
      double focalTreeProb = 0.;
      
      for(size_t k = 0; k < jumpFactor; ++k) {
          
        size_t hsIndex = k + j * jumpFactor;
        
        //cout << "index = " << hsIndex << " ";  
        focalTreeProb += postProbMatrix(hsIndex, i);
      }
      treePostProbs[i][j] = focalTreeProb;
    }
  }
  return treePostProbs; 
}

//gets posterior probability matrix considering only rate categories as hidden states
VVdouble MmPsmc::getRatePosteriorMatrix_(const VVdouble& postProbMatrix, const string& rate) {

  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t rateIndexInLibrary = rateIndexInVectors + 1;
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();
  size_t numCategories = mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getNumberOfCategories();

  VVdouble ratePostProbs(postProbMatrix.size(), Vdouble(numCategories)); 
  
  for(size_t i = 0; i < postProbMatrix.size(); ++i) { //for all sites
      
    for(size_t j = 0; j < numCategories; ++j) {
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
          
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) { 
          focalCatProb += postProbMatrix[i][k];
        }
      }
      
      ratePostProbs[i][j] = focalCatProb;
    } 
  }
  
  return ratePostProbs;
}

//gets posterior probability matrix considering only rate categories as hidden states
VVdouble MmPsmc::getRatePosteriorMatrix_(const zipHMM::Matrix& postProbMatrix, const string& rate) { 

  size_t rateIndexInVectors = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second;
  size_t rateIndexInLibrary = rateIndexInVectors + 1;
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();
  size_t numCategories = mmsmc_ -> getParameterScalings()[rateIndexInVectors] -> getNumberOfCategories();

  VVdouble ratePostProbs(postProbMatrix.get_width(), Vdouble(numCategories)); 
  
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) { //for all sites
      
    for(size_t j = 0; j < numCategories; ++j) {
        
      double focalCatProb = 0.;
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
          
        if(j == mmsmc_ -> getHmmStatesLibrary() -> getHiddenStates()[k][rateIndexInLibrary]) { 
          focalCatProb += postProbMatrix(k, i);  
        }
      }
      
      ratePostProbs[i][j] = focalCatProb;
    }
  }
  
  return ratePostProbs; 
}

//https://stackoverflow.com/questions/42871932/how-to-find-all-positions-of-an-element-using-stdfind?noredirect=1&lq=1
void MmPsmc::maskTreeSequence_(vector< size_t >& treeSeq, size_t start, size_t end) { 

  size_t state = mmsmc_ -> getNumberOfIntervals();

  //to avoid re-copying the fragment of interest
  vector< unsigned char >::iterator start_pos = biHaploidSnpCalling_.begin() +
                                                static_cast< vector< unsigned char >::difference_type >(start);
  vector< unsigned char >::iterator end_pos = biHaploidSnpCalling_.begin() +
                                              static_cast< vector< unsigned char >::difference_type >(end); 

  vector< unsigned char >::iterator iter = start_pos;

  while(iter != end_pos) { //looks for missing data in biHaploidSnpCalling_

    iter = find(iter, end_pos, 2u); //the next site with missing data

    if(iter != end_pos) {

      size_t pos = distance(start_pos, iter);
      //mask tree sequence by assigning the "missing" TMRCA state
      treeSeq[pos] = state;

      ++iter;
    }
  }
}
