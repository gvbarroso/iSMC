/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 05/08/2019
 *
 */


#include <thread>
#include <cmath>

#ifdef SIMPLEZIPHMM
  #include <SimpleZipHMM/forwarder.hpp>
  #include <SimpleZipHMM/posterior_decoding.hpp>
  #include <SimpleZipHMM/hmm_io.hpp>
#else
  #include <zipHMM/forwarder.hpp>
  #include <zipHMM/posterior_decoding.hpp>
  #include <zipHMM/hmm_io.hpp>
#endif //SIMPLEZIPHMM

#include "Psmc.h"

using namespace std;
using namespace bpp;
  
  
void Psmc::computeBiHaploidLogLikelihood(size_t numAvailThreads) { 

  //copy seqBreakpoints_ and sort it by length
  Breakpoints sortedBreakpoints = seqBreakpoints_; 
   
  sort(begin(sortedBreakpoints), end(sortedBreakpoints), [ ]( const pair< size_t, size_t >& bp1, const pair< size_t, size_t >& bp2 )
      { return (bp1.second - bp1.first) < (bp2.second - bp2.first); });
    
  //computes the logLikelihood of sequence resetting the chain at specified breakpoints
  size_t numberOfBreakpoints = sortedBreakpoints.size();  
  size_t batchSize = min(numberOfBreakpoints, numAvailThreads);
  size_t numBatches = (numberOfBreakpoints + batchSize - 1) / batchSize;
  
  //cout << "numberOfBreakpoints: " << numberOfBreakpoints << "; batchSize: " << batchSize << "; numBatches: " << numBatches << endl;
  
  double compLikBreakpoints = 0.; //sum of the logLikelihoods of each fragment (~chr)
  
  //computes likelihood per batch
  for(size_t i = 0; i < numBatches; ++i) {
    
    Breakpoints batchBps(0);
    
    for(size_t j = 0; j < batchSize; ++j) {
      
      size_t index = j + i * batchSize;
      
      if(index < numberOfBreakpoints) {
          
        size_t fragStart = sortedBreakpoints[index].first; 
        size_t fragEnd = sortedBreakpoints[index].second;  
      
        //cout << "i: " << i << "; j: " << j << "; index: " << index << "; start: " << fragStart << "; end: " << fragEnd << endl;
       
        batchBps.push_back(make_pair(fragStart, fragEnd)); 
      }
      
      else {
        break;
      }
    }
    
    compLikBreakpoints += computeBatchLogLikelihood_(batchBps);
  }
  
  biHaploidLogLikelihood_ = compLikBreakpoints;
}

vector< unsigned char > Psmc::fetchFragment(size_t genomicStart, size_t genomicEnd) {
    
  if (genomicStart >= biHaploidSnpCalling_.size())
    throw Exception("Psmc::fetchFragment. Starting position outside data range.");
  if (genomicEnd > biHaploidSnpCalling_.size())
    throw Exception("Psmc::fetchFragment. Ending position outside data range.");
  if (genomicStart > genomicEnd)
    throw Exception("Psmc::fetchFragment. Ending position before starting positions.");

  vector< unsigned char >::const_iterator start = biHaploidSnpCalling_.begin() + static_cast< vector< unsigned char >::difference_type >(genomicStart);
  vector< unsigned char >::const_iterator end = biHaploidSnpCalling_.begin() + static_cast< vector< unsigned char >::difference_type >(genomicEnd); 
  
  vector< unsigned char > selectedFragment(start, end);
    
  return selectedFragment;
}

void Psmc::selectFragment(size_t genomicStart, size_t genomicEnd) {
    
  vector< unsigned char >::const_iterator start = biHaploidSnpCalling_.begin() + static_cast<vector< unsigned char>::difference_type >(genomicStart);
  vector< unsigned char >::const_iterator end = biHaploidSnpCalling_.begin() + static_cast< vector<unsigned char>::difference_type >(genomicEnd);
  
  vector< unsigned char > selectedFragment(start, end);
  
  selectedFragment_ = selectedFragment;
  scalingFactorsForward_.resize(selectedFragment_.size());
}

void Psmc::writeDataStructures() {
    
  //creates a data structure for each breakpoint where the Markov Chain is reset in this diploid
  for(size_t i = 0; i < seqBreakpoints_.size(); ++i) {
      
    zipHMM::Forwarder preparator; 
    
    cout << "   creating data structure for fragment #" << i + 1 << " of " << seqBreakpoints_.size() << "." << endl;
    
    auto focalPair = seqBreakpoints_[i];
    
    //these are in computer coordinates (indexed by 0)
    size_t fragStart = get< 0 >(focalPair); 
    size_t fragEnd = get< 1 >(focalPair); 
        
    vector< unsigned char > focalSegment = fetchFragment(fragStart, fragEnd);

    size_t focalAlphabetSize = smcep_ -> fetchNumberOfObservedStates(focalSegment);
    
    //cout << "fragStart = " << fragStart << "; fragEnd = " << fragEnd << "; seglength = " << focalSegment.size() << endl;
    
    if(focalAlphabetSize != smcep_ -> getExpectedMatrix().front().size()) {
        
      cout << "WARNING!!! Number of observed states in fragment " << i + 1;
      cout << " does not match the number of observed states in the entire sequence!" << endl;
      cout << "This usually happens if the fragment is small and lacks either homozygous,";
      cout << "heterozygous or missing sites." << endl;
      cout << "Fragment size: " << focalSegment.size() << endl;
      map<unsigned char, size_t> counts = VectorTools::countValues(focalSegment);
      for(auto& x : counts)
      {
        cout << static_cast<int>(x.first) << "\t" << x.second << endl;
      }
      throw Exception("iSMC::Could not create data structure with zipHMM!");
      
    }
    
    //polymorphic version from SimpleZipHMM
    preparator.read_seq(focalSegment, focalAlphabetSize, smc_ -> getNumberOfHiddenStates(), 1e+4); 
    
    //writes data structure using human coordinates (starting from 1) in folder names
    preparator.write_to_directory("ziphmm_" + biHaploidName_ + "_" + TextTools::toString(fragStart) + "-" + TextTools::toString(fragEnd));
  }  
}

double Psmc::forwardAlgorithm() {
    
  //this forward algorithm is meant for Baum-Welch use
  size_t numHiddenStates = smc_ -> getNumberOfHiddenStates();
  size_t sequenceLength = selectedFragment_.size();
  forwardMatrix_.resize(sequenceLength);
  
  for(size_t i = 0; i < forwardMatrix_.size(); ++i) { 
    forwardMatrix_[i].resize(numHiddenStates);
  }
  
  //Initialization 
  double emissionProb = 0.;
  double scale = 0.;
  
  for(size_t j = 0; j < numHiddenStates; ++j) {
    emissionProb = smcep_ -> getExpectedMatrix()[j][selectedFragment_[0]];
    forwardMatrix_[0][j] = hiddenStatesInitializationProbabilities_[j] * emissionProb;
    scale += forwardMatrix_[0][j];
  }
  
  scalingFactorsForward_[0] = scale;
  
  for(size_t j = 0; j < numHiddenStates; ++j) {
    forwardMatrix_[0][j] /= scale;
  }
  
  //Forward Induction:
  for(size_t t = 1; t < sequenceLength; ++t) {
      
    scale = 0.; //restarts scale:
    
    for(size_t j = 0; j < numHiddenStates; ++j) {
        
      emissionProb = smcep_ -> getExpectedMatrix()[j][selectedFragment_[t]];
      
      double transitionProb = 0.; //(re)starts transitionProb:
      double integratedRate = 0.; //(re)starts a varibale to store the rate of going from state k to j
      
      for(size_t k = 0; k < numHiddenStates; ++k) {
          
        transitionProb = smctp_ -> getExpectedMatrix()[k][j]; //probability of having gone from state k to j:
        integratedRate += forwardMatrix_[t - 1][k] * transitionProb; //the different ways to go from state k at t-1 to state j at t, weighted by the previous steps
      }
      
      forwardMatrix_[t][j] = integratedRate * emissionProb;
      
      scale += forwardMatrix_[t][j];
    }
    
    for(size_t j = 0; j < numHiddenStates; ++j) { //scales the entries of column t in forwardMatrix_:
      forwardMatrix_[t][j] /= scale;
    }
    
    scalingFactorsForward_[t] = scale; 
  }
  
  //computez logLikelihood for BaumWelch purposes
  Vdouble tmp = scalingFactorsForward_;
  
  sort(tmp.begin(), tmp.end());
  
  double logLikelihood = 0.;
  
  for(size_t t = 0; t < sequenceLength; ++t) {
    logLikelihood += log(tmp[t]);
  }
  
  return logLikelihood;
}
    
void Psmc::backwardAlgorithm() { 
  //this forward algorithm is meant for posterior decoding (doesn't return logLikelihood)
  size_t sequenceLength = selectedFragment_.size();
  size_t numHiddenStates = smc_ -> getNumberOfHiddenStates();
  backwardMatrix_.resize(sequenceLength);
  for(size_t i = 0; i < backwardMatrix_.size(); ++i) { 
    backwardMatrix_[i].resize(numHiddenStates);
  }
  //Initialization 
  for(size_t i = 0; i < numHiddenStates; ++i) { 
    backwardMatrix_[sequenceLength - 1][i] = 1.;
  }
  //Backward Induction
  //going backwards for every position t - 1 in the sequence:
  for(size_t t = sequenceLength - 1; t > 0; --t) {
    unsigned int observation = selectedFragment_[t];
    //for each possible (hidden) state j in the "current" step t-1:
    for(size_t j = 0; j < numHiddenStates; ++j) {
      //(re)sets weighted prob. of going from state j to k, and then dropping observation:
      double integratedRate = 0.;
      //for each possible state k in the immediate next step t (ie, in the "future past"):
      for(size_t k = 0; k < numHiddenStates; ++k) {
        //prob. of having gone from state j to k:
        double transitionProb = smctp_ -> getExpectedMatrix()[j][k];
        //prob. of emitting the observation at step t:
        double emissionProb = smcep_ -> getExpectedMatrix()[k][observation];
        integratedRate += transitionProb * emissionProb * backwardMatrix_[t][k];
      }
      backwardMatrix_[t - 1][j] = integratedRate;
      //updates probs. using the same scales used in the forward algorithm, one position ahead:
      backwardMatrix_[t - 1][j] /= scalingFactorsForward_[t];
    }
  }
}  
  
void Psmc::performPosteriorDecoding(size_t genomicStart, size_t genomicEnd) {
    
  //reduces the bihaploid sequence to the fragment of interest
  selectFragment(genomicStart, genomicEnd);
  
  //updates arrays with current parameters
  forwardAlgorithm();
  backwardAlgorithm();
  
  computePosteriorProbabilities();
    
  reconstructedHiddenStates_.resize(selectedFragment_.size());
  
  for(size_t t = 0; t < selectedFragment_.size(); ++t) {
      
    //index of the hidden state with highest prob.:
    size_t stateIndex = 0;
    
    //current highest prob. in posteriorProbabilitiesVector:
    double highestProb = posteriorProbMatrix_[t][0];
    
    for(size_t i = 1; i < smc_ -> getNumberOfHiddenStates(); ++i) {
        
      if(posteriorProbMatrix_[t][i] > highestProb) { 
          
        highestProb = posteriorProbMatrix_[t][i];
        
        stateIndex = i;
      } 
    }
    reconstructedHiddenStates_[t] = static_cast< unsigned char >(stateIndex);
  }
  
  //frees memory after we are done:
  VVdouble().swap(forwardMatrix_);
  VVdouble().swap(backwardMatrix_);
  Vdouble().swap(scalingFactorsForward_);
}
  
void Psmc::posteriorDecodingUsingZipHMM(size_t genomicStart, size_t genomicEnd) {
    
  selectFragment(genomicStart, genomicEnd);
  
  posteriorProbMatrix_.reserve(selectedFragment_.size());
  
  for(size_t i = 0; i < selectedFragment_.size(); ++i) {
    posteriorProbMatrix_[i].reserve(smc_ -> getNumberOfHiddenStates());
  }
    
  zipHMM::Matrix tempMatrix;
  
  zipHMM::posterior_decoding(selectedFragment_,
                             hiddenStatesInitializationProbabilities_,
                             smctp_ -> getExpectedMatrix(),
                             smcep_ -> getExpectedMatrix(),
                             tempMatrix);
  
  for(size_t i = 0; i < selectedFragment_.size(); ++i) {
      
    for(size_t j = 0; j < smc_ -> getNumberOfHiddenStates(); ++j) {
        
      posteriorProbMatrix_[i][j] = tempMatrix(j, i); 
    }
  }
} 
  
void Psmc::computePosteriorProbabilities() {
  
  posteriorProbMatrix_.resize(selectedFragment_.size());
  
  for(size_t i = 0; i < selectedFragment_.size(); ++i) {
    posteriorProbMatrix_[i].resize(smc_ -> getNumberOfHiddenStates());
  }
  for(size_t t = 0; t < selectedFragment_.size(); ++t) {
    for(size_t j = 0; j < smc_ -> getNumberOfHiddenStates(); ++j) {
      posteriorProbMatrix_[t][j] = forwardMatrix_[t][j] * backwardMatrix_[t][j]; 
    }
  }
}


double Psmc::computeBatchLogLikelihood_(const Breakpoints& bpBatch) {
  
  size_t numberOfTasks = bpBatch.size(); 
      
  double compLikBatch = 0.; //sum of the logLikelihoods of each fragment in batch

  auto computeLikelihood = [&] (size_t breakpoint_id) {
      
    zipHMM::Forwarder calculator;   
    
    size_t fragStart = bpBatch[breakpoint_id].first; 
    size_t fragEnd = bpBatch[breakpoint_id].second;  
    
    //cout << "fragStart: " << fragStart << "; fragEnd: " << fragEnd << endl;
        
    calculator.read_from_directory("ziphmm_" + biHaploidName_ + "_" + TextTools::toString(fragStart) + "-" + TextTools::toString(fragEnd)); 
    
    compLikBatch += calculator.forward(hiddenStatesInitializationProbabilities_,
                                       smctp_ -> getExpectedMatrix(),
                                       smcep_ -> getExpectedMatrix());
  };

  thread* threadVector = new thread[numberOfTasks];
  
  for(size_t i = 0; i < numberOfTasks; ++i) {
    threadVector[i] = thread(computeLikelihood, i);
  }
  
  for(size_t i = 0; i < numberOfTasks; ++i) {
    threadVector[i].join();
  }
  
  delete [] threadVector; 
  
  return compLikBatch;
}

vector< pair< size_t, size_t > > Psmc::selectInformativeRegions_(size_t missingThreshold) { 
  
  //first identifies all uninformative blocks (contiguous missing data) in the whole sequence
  //vector with coordinates of UNINFORMATIVE blocks (to be removed from INFORMATIVE segments)
  vector< pair < size_t, size_t > > infoGaps(0);
  
  cout << "   setting extra breakpoints to skip " << missingThreshold << "-bp long uninformative regions." << endl;
  
  size_t blockLength = 0; //length of uninformative block
  for(size_t i = 0; i < biHaploidSnpCalling_.size(); ++i) {
      
    if(biHaploidSnpCalling_[i] == 2u) {
        
      if(i < biHaploidSnpCalling_.size() - 1) {
          
        ++blockLength;
      }
      
      else { //if last position in sequence is uninformative
          
        ++blockLength;
        
        if(blockLength >= missingThreshold) { //and is part of a large enough block
            
          size_t blockStart = i - blockLength + 1; //+1 because we identified a block end while sitting on it
          
          size_t blockEnd = i;
          
          infoGaps.push_back(make_pair(blockStart, blockEnd));
        }
      }
    }
    
    else { //if block of missing data has ended
        
      if(blockLength >= missingThreshold) { 
          
        size_t blockStart = i - blockLength;
        size_t blockEnd = i - 1;
        
        infoGaps.push_back(make_pair(blockStart, blockEnd));
      }
      
      blockLength = 0; //reset block length
    }
  }

  //now merges seqBreakpoints_ and infoGaps, dismantling the pairs 
  vector< size_t > unpairedCoordinates(0);
  
  for(size_t i = 0; i < seqBreakpoints_.size(); ++i) {
    
    auto focalPair = seqBreakpoints_[i];
    
    size_t first = get< 0 >(focalPair);
    unpairedCoordinates.push_back(first);
    
    size_t second = get< 1 >(focalPair);
    unpairedCoordinates.push_back(second);
  }
  
  for(size_t i = 0; i < infoGaps.size(); ++i) {
    
    auto focalPair = infoGaps.at(i);
    
    size_t first = get< 0 >(focalPair);
    unpairedCoordinates.push_back(first);
    
    size_t second = get< 1 >(focalPair);
    unpairedCoordinates.push_back(second);
  }
  
  sort(unpairedCoordinates.begin(), unpairedCoordinates.end());
  
  //assembles pairs of INFORMATIVE blocks
  vector< pair < size_t, size_t > > infoBlocks(0);
  
  for(size_t i = 0; i <= unpairedCoordinates.size() - 2; i+= 2) {
      
    size_t fragStart = unpairedCoordinates[i];
    size_t fragEnd = unpairedCoordinates[i + 1];
    
    //handles short segments of informative sites betwen two uninformative blocks
    if((fragEnd - fragStart) > 1000000) { //NOTE: only includes them if long enough
        
      infoBlocks.push_back(make_pair(fragStart, fragEnd));
    }
  }
  
  return infoBlocks;
}

Vdouble Psmc::fetchPerSiteShannonEquitability_(const VVdouble& postProbMatrix) {
  
  Vdouble shannonPerSite(postProbMatrix.size());
  
  double normFactor = log(postProbMatrix.front().size());
  for(size_t i = 0; i < postProbMatrix.size(); ++i) {
    double eqIndex = 0.;
    for(size_t j = 0; j < postProbMatrix[i].size(); ++j) {
      eqIndex -= postProbMatrix[i][j] * log(postProbMatrix[i][j]);
    }
    shannonPerSite[i] = eqIndex / normFactor;
  }
  
  return shannonPerSite;
}

void Psmc::computePerSiteShannonEquitability_(const VVdouble& postProbMatrix, const string& sequenceFileName) {
  
  ofstream shannonFile;
  shannonFile.open(sequenceFileName + "_shannon_eq.txt.gz", 
                   std::ios_base::out | std::ios_base::binary);
  
  boost::iostreams::filtering_ostream shannonStream; 
  shannonStream.push(boost::iostreams::gzip_compressor());
  shannonStream.push(shannonFile);
  
  double normFactor = log(postProbMatrix.front().size());
  for(size_t i = 0; i < postProbMatrix.size(); ++i) {
    double eqIndex = 0.;
    for(size_t j = 0; j < postProbMatrix[i].size(); ++j) {
      eqIndex -= postProbMatrix[i][j] * log(postProbMatrix[i][j]);
    }
    double shannonEq = eqIndex / normFactor;
    shannonStream << shannonEq << endl;
  }
  
  boost::iostreams::close(shannonStream);
}
  
Vdouble Psmc::fetchPerSiteShannonEquitability_(const zipHMM::Matrix& postProbMatrix) {
  
  Vdouble shannonPerSite(postProbMatrix.get_width());
  double normFactor = log(postProbMatrix.get_width());
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) {
    double eqIndex = 0.;
    for(size_t j = 0; j < postProbMatrix.get_height(); ++j) {
      eqIndex -= postProbMatrix(j, i) * log(postProbMatrix(j, i));
    }
    shannonPerSite[i] = eqIndex / normFactor;
  }
  
  return shannonPerSite;
}

void Psmc::computePerSiteShannonEquitability_(const zipHMM::Matrix& postProbMatrix, const string& sequenceFileName) {
  
  ofstream shannonFile;
  shannonFile.open(sequenceFileName + "_shannon_eq.txt.gz", 
                   std::ios_base::out | std::ios_base::binary);
  
  boost::iostreams::filtering_ostream shannonStream; 
  shannonStream.push(boost::iostreams::gzip_compressor());
  shannonStream.push(shannonFile);

  double normFactor = log(postProbMatrix.get_width());
  for(size_t i = 0; i < postProbMatrix.get_width(); ++i) {
    double eqIndex = 0.;
    for(size_t j = 0; j < postProbMatrix.get_height(); ++j) {
      eqIndex -= postProbMatrix(j, i) * log(postProbMatrix(j, i));
    }
    double shannonEq = eqIndex / normFactor;
    shannonStream << shannonEq << endl;
  }
  
  boost::iostreams::close(shannonStream);
}
