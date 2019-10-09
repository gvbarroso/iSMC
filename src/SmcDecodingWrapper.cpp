/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 09/07/2018
 * Last modified: 09/09/2019
 *
 */

#include <stdlib.h> 
#include <fstream>
#include <thread>
#include <chrono>
#include <utility>
    
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>

#include "ParameterCategoryTransitions.h"
#include "SmcDecodingWrapper.h"
#include "Global.h"

using namespace std;
using namespace bpp;


void SmcDecodingWrapper::writePosteriorCoordinatesToFile() {

  ofstream postDecCoordinates; //coordinates of the fragments

  postDecCoordinates.open(smcOptions_ -> getLabel() + "_decode_coordinates.txt");
  postDecCoordinates << "frag_start" << "\t" << "frag_end" << "\t" << "seq_block" << "\t" << "file_index" << endl;

  size_t fragSize = smcOptions_ -> getFragmentSize();
  size_t numBreakpoints = mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints().size();
    
  for(size_t k = 0; k < numBreakpoints; ++k) {
    
    auto focalPair = mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints()[k];

    size_t breakpointStart = get< 0 >(focalPair);
    size_t breakpointEnd = get< 1 >(focalPair);
    size_t breakpointLength = breakpointEnd - breakpointStart;
    size_t numFragments = (breakpointLength + fragSize - 1) / fragSize;
    
    for(size_t i = 0; i < numFragments; ++i) {

      size_t start = i * fragSize + breakpointStart;
      size_t end = start + fragSize;
      
      //in case breakpointLength is not a multiple of fragSize
      if(end > breakpointEnd) {
        end = breakpointEnd + 1;
      }
      
      postDecCoordinates << start + 1;
      postDecCoordinates << "\t";

      postDecCoordinates << end;
      postDecCoordinates << "\t";
      
      postDecCoordinates << k + 1;
      postDecCoordinates << "\t";
      
      postDecCoordinates << i + 1;
      postDecCoordinates << endl;
    }
  }

  postDecCoordinates.close();
}

void SmcDecodingWrapper::writeDiploidLabelsToFile() {

  ofstream diploidLabels;

  diploidLabels.open(smcOptions_ -> getLabel() + "_diploid_decoding_labels.txt");
  
  diploidLabels << "decoding_label" << "\t" << "diploid_id" << endl;

  for(size_t i = 0; i < mPsmc_ -> getWholePsmcVector().size(); ++i) {

    diploidLabels << smcOptions_ -> getLabel() + "_diploid_" + TextTools::toString(i + 1); //indexed from 1
    diploidLabels << "\t";
    
    diploidLabels << mPsmc_ -> getWholePsmcVector()[i] -> getBiHaploidName();
    diploidLabels << endl;
  }

  diploidLabels.close();
}

//NOTE for a single pair of genomes
void SmcDecodingWrapper::decodeDiploid() { 

  writeDiploidLabelsToFile();
  writePosteriorCoordinatesToFile();
  
  size_t numBreakpoints = mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints().size();
  size_t batchSize = 1; //default
  
  if(smcOptions_ -> decodeBreakpointsParallel()) { //default = false since it requires a lot of RAM
    batchSize = min(numBreakpoints, NUMBER_OF_AVAILABLE_THREADS);
  }
    
  size_t numBatches = (numBreakpoints + batchSize - 1) / batchSize;
    
  //computes likelihood per batch of breakpoints (~chr)
  for(size_t i = 0; i < numBatches; ++i) {
      
    vector< pair< size_t, size_t > > batchBreakpoints(0);
    vector< size_t > indexVector(0); 
    
    for(size_t j = 0; j < batchSize; ++j) {
      
      size_t index = j + i * batchSize;
            
      if(index < numBreakpoints) {
          
        indexVector.push_back(index);
        batchBreakpoints.push_back(mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints()[index]); 
      }
      
      else {
        break;
      }
    }
    
    decodeBatchOfBreakpoints_(batchBreakpoints, indexVector);
  }
}

void SmcDecodingWrapper::computeAverageLandscapes() {
  
  for(size_t i = 0; i < mmsmc_ -> getParameterTransitions().size(); ++i) {
    
    if(mmsmc_ -> getParameterTransitions()[i] -> getNumberOfCategories() > 1) {
    
      string rate = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().left.find(i) -> second;
      computeAverageRateOverDiploids_(rate);
    }
  }
}

Vdouble SmcDecodingWrapper::readLandscapeFromFile(const string& file) { 

  Vdouble rateLandscape;
  
  boost::iostreams::filtering_istream landStream; 
  landStream.push(boost::iostreams::gzip_decompressor());
  
  ifstream landFile;
  landFile.open(file, std::ios_base::in | std::ios_base::binary);
  
  landStream.push(landFile);

  if(landFile.is_open()) {

    string line;
    vector< string > splitLine;
    
    while(getline(landStream, line)) {

      boost::split(splitLine, line, [](char c){ return c == '\t'; });
      size_t ratePos = 0;
      if(smcOptions_ -> restrictedModulation()) { //time-restricted landscapes are formatted per line: TMRCA \t rateVal
        ratePos = 1;
      }
      double rateVal = stod(splitLine[ratePos]); 
      rateLandscape.push_back(rateVal);
    }

    boost::iostreams::close(landStream);
  }

  else {
    throw Exception("iSMC::could not open landscape file: " + file);
  }
  
  return rateLandscape;
}

void SmcDecodingWrapper::writeSingleNucleotideLandscapeToFile(const Vdouble& landscape, const string& file) {

  boost::iostreams::filtering_ostream mapStream; 
  
  ofstream mapFile;
  mapFile.open(file + ".gz", std::ios_base::out | std::ios_base::binary);  
  
  mapStream.push(boost::iostreams::gzip_compressor());
  mapStream.push(mapFile);

  for(auto& value : landscape) {
    mapStream << value << endl;
  }

  boost::iostreams::close(mapStream);
}

//NOTE for multiple pairs of genomes
void SmcDecodingWrapper::decodeDataset(bool missing) {

  writePosteriorCoordinatesToFile();  
  writeDiploidLabelsToFile();

  size_t numAvailThreads = NUMBER_OF_AVAILABLE_THREADS; //copies from global variable
    
  //all diploids have the same BPs
  size_t numBreakpoints = mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints().size(); 
  //creates batches 
  size_t batchSize = min(numBreakpoints, numAvailThreads);
  size_t numBatches = (numBreakpoints + batchSize - 1) / batchSize; 
  
  //if we are using the threads to decode breakpoints in parallel
  if(smcOptions_ -> decodeBreakpointsParallel()) {
    //updates numAvailThreads to be used by each diploid (will be forwarded down)
    numAvailThreads /= batchSize;
  }
  
  //decodes breakpoints in batches
  for(size_t i = 0; i < numBatches; ++i) {
      
    vector< pair< size_t, size_t > > batchBreakpoints(0);
    vector< size_t > indexVector(0); 
    
    for(size_t j = 0; j < batchSize; ++j) {
      
      size_t index = j + i * batchSize;

      if(index < numBreakpoints) {
          
        indexVector.push_back(index);
        //all diploids have the same BPs so we just get the one from the first diploid (index = 0)
        batchBreakpoints.push_back(mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints()[index]); 
      }
      
      else {
        break;
      }
    }
    
    decodeBatchOfBreakpoints_(numAvailThreads, missing, batchBreakpoints, indexVector);
  }
}

void SmcDecodingWrapper::jointlyDecodeUsingBaumWelch(VVdouble& rateEmiss, VVdouble& rateTrans, const string& rate, bool missing) {

  //size_t fragSize = smcOptions_ -> getFragmentSize();
  //size_t numHs = mmsmc_ -> getNumberOfHiddenStates();
  size_t numDiploids = mPsmc_ -> getWholePsmcVector().size();
  size_t seqLength = mPsmc_ -> getWholePsmcVector().front() -> getSequenceLength(); //all diploids: same seqLength
  size_t rateIndex = mmsmc_ -> getHmmStatesLibrary() -> getParameterAlphabet().right.find(rate) -> second; 

  clock_t t1, t2;
  double diff, runTime = -1.; 

  cout << "Preparing Baum-Welch dedicated optimisation of " << rate << " under a joint-model." << endl;
  
  vector< vector < size_t > > wholeSeqObs(numDiploids, vector< size_t >(seqLength));
  
  if(true) { //fake scope to let entireTreeSeq disappear

    t1 = clock(); //decodes TMRCA in all diploids
    vector< vector < size_t > > entireTreeSeq = fetchWholeTmrcaDist_(); 
    t2 = clock();

    diff = static_cast< double >(t2) - static_cast< double >(t1);
    runTime = diff / CLOCKS_PER_SEC;

    cout << "Run-time get full TMRCA dist: " << runTime << " seconds." << endl;
    cout << "Length of TMRCA dist. = " << entireTreeSeq.front().size() << endl;

    t1 = clock();

    for(size_t i = 0; i < numDiploids; ++i) {

      if(rate == "rho") {
        wholeSeqObs[i] = mPsmc_ -> getWholePsmcVector()[i] -> fetchTransitionSequence(entireTreeSeq[i], 0, seqLength, missing);
      }

      else if(rate == "theta") {
        wholeSeqObs[i] = mPsmc_ -> getWholePsmcVector()[i] -> fetchEmissionSequence(entireTreeSeq[i], 0, seqLength);   
      }
    }

    t2 = clock();
    diff = static_cast< double >(t2) - static_cast< double >(t1);
    runTime = diff / CLOCKS_PER_SEC;

    cout << "Run-time get new observations: " << runTime << " seconds." << endl;
    cout << "Length of new obs. seq. = " << wholeSeqObs.front().size() << endl;
  }
  
  // Baum-Welch optimisation of HMM matrices
  // Filtering by Shannon Eq. might make more sense after BW

  //sets matrices (check if they need to be reset at every step)
  size_t numCategories = mmsmc_ -> getParameterTransitions()[rateIndex] -> getNumberOfCategories(); 

  VVdouble rateFwdMat(seqLength, Vdouble(numCategories)); 
  VVdouble rateBckMat(seqLength, Vdouble(numCategories)); 
  
  Vdouble fwdScales(seqLength);
  Vdouble piRate(numCategories);

  for(size_t j = 0; j < numCategories; ++j) {
    piRate[j] = 1. / static_cast< double >(numCategories);
  }
  
  BaumWelch baumWelcher;

  size_t numExpMaxSteps = 20;

  for(size_t i = 0; i < numExpMaxSteps; ++i) {
      
    cout << "Performing EM step #" << i << " of " << numExpMaxSteps - 1 << "..." << endl; //cout.flush();
     
    //forward-backward 
    t1 = clock();
    mPsmc_ -> jointForward(wholeSeqObs, rateTrans, piRate, rateEmiss, fwdScales, rateFwdMat);
    t2 = clock();

    diff = static_cast< double >(t2) - static_cast< double >(t1);
    runTime = diff / CLOCKS_PER_SEC;
    
    cout << "Run-time joint-forward: " << runTime << " seconds." << endl;
    
    t1 = clock();
    mPsmc_ -> jointBackward(wholeSeqObs, rateTrans, rateEmiss, fwdScales, rateBckMat);
    t2 = clock();

    diff = static_cast< double >(t2) - static_cast< double >(t1);
    runTime = diff / CLOCKS_PER_SEC;
    
    cout << "Run-time joint-backward: " << runTime << " seconds." << endl;

    t1 = clock();    
    baumWelcher.maximiseRateTransMat(fwdScales, rateFwdMat, rateBckMat, rateTrans, rateEmiss, wholeSeqObs);
    t2 = clock();

    diff = static_cast< double >(t2) - static_cast< double >(t1);
    runTime = diff / CLOCKS_PER_SEC;
    
    cout << "Run-time Baum-Welch TransMat: " << runTime << " seconds." << endl;
    
    /*t1 = clock();
    baumWelcher.maximiseRateEmissMat(rateFwdMat, rateBckMat, rateEmiss, wholeSeqObs);
    t2 = clock();
    diff = static_cast< double >(t2) - static_cast<double>(t1);
    runTime = diff / CLOCKS_PER_SEC;
    cout << "Run-time Baum-Welch EmissMat: " << runTime << " seconds." << endl;
    */
  }  

  //posterior decoding using HMM matrices optimised w/ Baum-Welch  
  cout << endl << "Baum-Welch is done. Computing posterior average..."; cout.flush();
  
  double rateVal = mmsmc_ -> getParameterValue(rate);
  
  string fileName = smcOptions_ -> getLabel() + "_baum-welch_" + rate + ".txt.gz";
  ofstream decFile;
  decFile.open(fileName, std::ios_base::out | std::ios_base::binary);
  
  boost::iostreams::filtering_ostream bwStream; 
  bwStream.push(boost::iostreams::gzip_compressor());
  bwStream.push(decFile);

  t1 = clock();

  for(size_t i = 0; i < seqLength; ++i) {

    double postAvg = 0.;

    for(size_t j = 0; j < numCategories; ++j) { 
      postAvg += rateFwdMat[i][j] * rateBckMat[i][j] * rateVal * mmsmc_ -> getParameterScalings()[rateIndex] -> getCategories()[j];
    }
    
    bwStream << postAvg << endl;
  }

  boost::iostreams::close(bwStream);
  
  cout << " done." << endl;
  
  t2 = clock();
  diff = static_cast< double >(t2) - static_cast< double >(t1);
  runTime = diff / CLOCKS_PER_SEC;
  
  cout << "Run-time post. avg. Baum-Welch: " << runTime << " seconds." << endl;
}

//returns TMRCA landscape for fragment of interest in all diploids
vector< vector< size_t > > SmcDecodingWrapper::fetchLocalTmrcaDist_(VVdouble& piMmSmc, size_t fragStart, size_t fragEnd, size_t fragId) {
  
  size_t fragLength = fragEnd - fragStart + 1;
  size_t numHiddenStates = mmsmc_ -> getNumberOfHiddenStates();
  size_t numDiploids = mPsmc_ -> getWholePsmcVector().size();

  vector< vector < size_t > > treeObs(numDiploids, vector< size_t >(fragLength)); //ret

  auto extractDiploidTrees = [&] (size_t diploid_id) { 

    zipHMM::Matrix postProbMatrix;

    string diploidFile = smcOptions_ -> getLabel() + "_diploid_" + TextTools::toString(diploid_id + 1) + "_file_index_" + TextTools::toString(fragId + 1);

    mPsmc_ -> getWholePsmcVector()[diploid_id] -> computePosteriorProbs(fragStart, fragEnd, postProbMatrix, piMmSmc[diploid_id]); 
    
    //fills piMmSmc for next fragment
    for(size_t j = 0; j < numHiddenStates; ++j) {
      piMmSmc[diploid_id][j] = postProbMatrix(j, postProbMatrix.get_width() - 1);
    }
  
    //decodes trees
    treeObs[diploid_id] = mPsmc_ -> getWholePsmcVector()[diploid_id] -> fetchLocalTrees(postProbMatrix);

  };

  if(smcOptions_ -> decodeDiploidsParallel()) {
      
    thread* threadVector = new thread[numDiploids];
    
    for(size_t i = 0; i < numDiploids; ++i) {
      threadVector[i] = thread(extractDiploidTrees, i);
    }
    
    for(size_t i = 0; i < numDiploids; ++i) {
      threadVector[i].join();
    }
    
    delete [] threadVector; 
  }
  
  else {
      
    for(size_t i = 0; i < numDiploids; ++i) {
      extractDiploidTrees(i);
    }
  }

  return treeObs;
}

//gets TMRCA landscapes for the entire sequence in all diploids
vector< vector < size_t > > SmcDecodingWrapper::fetchWholeTmrcaDist_() {
  
  size_t fragSize = smcOptions_ -> getFragmentSize();
  size_t numDiploids = mPsmc_ -> getWholePsmcVector().size();
  size_t numHs = mmsmc_ -> getNumberOfHiddenStates();

  //init. probs. for the Markov-modulated SMC for each diploid
  VVdouble piMmSmc(numDiploids, Vdouble(numHs));
  for(size_t i = 0; i < numDiploids; ++i) { 
    for(size_t j = 0; j < numHs; ++j) {
      piMmSmc[i][j] = 1. / static_cast< double >(numHs);
    }
  } 

  vector< vector < size_t > > wholeSeqTmrcas(numDiploids, vector< size_t >(0)); //ret 

  //NOTE assumes all diploids share the same breakpoints
  size_t numBreakpoints = mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints().size();

  for(size_t k = 0; k < numBreakpoints; ++k) {
    
    auto focalPair = mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints()[k];

    size_t breakpointStart = get< 0 >(focalPair);
    size_t breakpointEnd = get< 1 >(focalPair);
    size_t breakpointLength = breakpointEnd - breakpointStart + 1;
    size_t numFragments = (breakpointLength + fragSize - 1) / fragSize;

    for(size_t i = 0; i < numFragments; ++i) {

      size_t start = i * fragSize + breakpointStart;
      size_t end = start + fragSize;
      
      //in case breakpointLength is not a multiple of fragSize
      if(end > breakpointEnd) {
        end = breakpointEnd;
      }

      //decode w/ max. prob. in all diploids
      vector< vector < size_t > > focalTmrca = fetchLocalTmrcaDist_(piMmSmc, start, end, i);

      //appends to obs. for the entire sequence 
      for(size_t j = 0; j < numDiploids; ++j) { 
        move(focalTmrca[j].begin(), focalTmrca[j].end(), back_inserter(wholeSeqTmrcas[j]));
      }
    }
  }

  return wholeSeqTmrcas;
}

void SmcDecodingWrapper::jointlyDecodeFragment_(VVdouble& piMmSmc, size_t fragStart, size_t fragEnd, size_t blockId, size_t fragId, size_t numAvailThreads, bool missing) {
  
  size_t fragLength = fragEnd - fragStart; 
  size_t numDiploids = mPsmc_ -> getWholePsmcVector().size();
  size_t numCatTheta = mmsmc_ -> getParameterTransitions()[0] -> getNumberOfCategories();
  size_t numCatRho = mmsmc_ -> getParameterTransitions()[1] -> getNumberOfCategories();
  bool timeRestricted = smcOptions_ -> restrictedModulation();
  
  //for all diploids, for the fragment of interest, gets new type of obs. states (for each rate with > 1 category)
  vector< vector < size_t > > rhoObsVector(numDiploids, vector< size_t >(fragLength));
  vector< vector < size_t > > thetaObsVector(numDiploids, vector< size_t >(fragLength));

  auto extractInfoDiploid = [&] (size_t diploid_id) {
      
    //gets posterior matrix
    zipHMM::Matrix postProbMatrix;
    
    string diploidFile = smcOptions_ -> getLabel() + "_diploid_" + TextTools::toString(diploid_id + 1) + "_block_" +
                         TextTools::toString(blockId) + "_file_index_" + TextTools::toString(fragId + 1);
    
    mPsmc_ -> getWholePsmcVector()[diploid_id] -> posteriorDecodingUsingZipHMM(fragStart, fragEnd, diploidFile, timeRestricted, postProbMatrix, piMmSmc[diploid_id]); 
    
    vector< size_t > treeSeq = mPsmc_ -> getWholePsmcVector()[diploid_id] -> fetchLocalTrees(postProbMatrix);
    
    //fills piMmSmc for next fragment
    for(size_t j = 0; j < piMmSmc[diploid_id].size(); ++j) { 
      piMmSmc[diploid_id][j] = postProbMatrix(j, postProbMatrix.get_width() - 1);
    }
    
    if(numCatTheta > 1) {
      thetaObsVector[diploid_id] = mPsmc_ -> getWholePsmcVector()[diploid_id] -> fetchEmissionSequence(treeSeq, fragStart, fragEnd); 
    }
    
    if(numCatRho > 1) {
      rhoObsVector[diploid_id] = mPsmc_ -> getWholePsmcVector()[diploid_id] -> fetchTransitionSequence(treeSeq, fragStart, fragEnd, missing);
    }
    
    cout << "   done for diploid #" << diploid_id << "." << endl; 
  };
  
  cout << "Extracting posterior information from all diploids..." << endl;

  
  if(smcOptions_ -> decodeDiploidsParallel()) {
      
    size_t batchSize = min(numDiploids, numAvailThreads);
    size_t numBatches = (numDiploids + batchSize - 1) / batchSize;  
    
    //organises because last batch may have less than batchSize diploids
    vector< size_t > sizeOfBatches(numBatches); 
    for(size_t i = 0; i < numBatches - 1; ++i) {
      sizeOfBatches[i] = batchSize;
    }
    sizeOfBatches[numBatches - 1] = numDiploids - batchSize * (numBatches - 1);
    
    for(size_t i = 0; i < numBatches; ++i) {
      
      size_t numTasks = sizeOfBatches[i];
      thread* threadVector = new thread[numTasks]; 
      
      for(size_t j = 0; j < numTasks; ++j) {

        size_t diploidIndex = j + i * numTasks;
        threadVector[j] = thread(extractInfoDiploid, diploidIndex);
      }
    
      for(size_t j = 0; j < numTasks; ++j) {
        threadVector[j].join();
      }
    
      delete [] threadVector; 
    }
  }
  
  else {
      
    for(size_t i = 0; i < numDiploids; ++i) {
      extractInfoDiploid(i);
    }
  }

  //WARNING the following is ugly
    
  if(numCatRho > 1) {
      
    cout << "Performing the forward-backward rho hierarchical HMM..."; cout.flush();
    
    mmsmc_ -> getHmmStatesLibrary() -> initRhoEmissionsAlphabet(smcOptions_ -> getNumberOfIntervals(), missing); 
    //emission probs from rho's to "hierarchical" observations
    VVdouble rhoEmiss = ep_ -> fetchRhoEmissions(tp_ -> fetchCompositeTransitionMatrix(missing));
 
    VVdouble rateFwdMat(fragLength, Vdouble(numCatRho)); 
    Vdouble forwardScales(fragLength);
    Vdouble piRate(numCatRho);
      
    for(size_t k = 0; k < numCatRho; ++k) {
      piRate[k] = 1. / static_cast< double >(numCatRho);
    }
        
    mPsmc_ -> jointForward(rhoObsVector, mmsmc_ -> getParameterTransitions()[1] -> getCategoryTransitions(),
                           piRate, rhoEmiss, forwardScales, rateFwdMat);

    VVdouble rateBackMat(fragLength, Vdouble(numCatRho)); 
    
    mPsmc_ -> jointBackward(rhoObsVector, mmsmc_ -> getParameterTransitions()[1] -> getCategoryTransitions(),
                            rhoEmiss, forwardScales, rateBackMat);
  
    cout << " done." << endl << "Computing posterior rho for joint model..."; cout.flush();
  
    string fileName = smcOptions_ -> getLabel() + "_joint_rho" + "_block_" + TextTools::toString(blockId) +
                      "_file_index_" + TextTools::toString(fragId + 1) + ".txt.gz";
                    
    ofstream decFile;
    decFile.open(fileName, std::ios_base::out | std::ios_base::binary);
  
    boost::iostreams::filtering_ostream jointStream; 
    jointStream.push(boost::iostreams::gzip_compressor());
    jointStream.push(decFile);

    double rateVal = mmsmc_ -> getParameterValue("rho"); //genome-wide avg

    for(size_t i = 0; i < fragLength; ++i) {
      
      double avgRate = 0.; //avg rate at site i
    
      for(size_t k = 0; k < numCatRho; ++k) {
        
        double focalCatProb = rateFwdMat[i][k] * rateBackMat[i][k];
        avgRate += focalCatProb * rateVal * mmsmc_ -> getParameterScalings()[1] -> getCategories()[k];
      }
    
      jointStream << avgRate << endl;
    }
  
    boost::iostreams::close(jointStream);
  
    cout << " done." << endl;
  }

  if(numCatTheta > 1) {
      
    cout << "Performing the forward-backward theta hierarchical HMM..."; cout.flush();
    
    mmsmc_ -> getHmmStatesLibrary() -> initThetaEmissionsAlphabet(smcOptions_ -> getNumberOfIntervals(), ep_ -> getNumberOfObservedStates()); 
    //emission probs from theta's to "hierarchical" observations
    VVdouble thetaEmiss = ep_ -> fetchThetaEmissions(ep_ -> fetchCompositeEmissionMatrix());
 
    VVdouble rateFwdMat(fragLength, Vdouble(numCatTheta)); 
    Vdouble forwardScales(fragLength);
    Vdouble piRate(numCatTheta);
      
    for(size_t k = 0; k < numCatTheta; ++k) {
      piRate[k] = 1. / static_cast< double >(numCatTheta);
    }
        
    mPsmc_ -> jointForward(thetaObsVector, mmsmc_ -> getParameterTransitions()[0] -> getCategoryTransitions(),
                           piRate, thetaEmiss, forwardScales, rateFwdMat);

    VVdouble rateBackMat(fragLength, Vdouble(numCatTheta)); 
    
    mPsmc_ -> jointBackward(thetaObsVector, mmsmc_ -> getParameterTransitions()[0] -> getCategoryTransitions(),
                            thetaEmiss, forwardScales, rateBackMat);
  
    cout << " done." << endl << "Computing posterior theta for joint model..."; cout.flush();
  
    string fileName = smcOptions_ -> getLabel() + "_joint_theta" + "_block_" + TextTools::toString(blockId) +
                      "_file_index_" + TextTools::toString(fragId + 1) + ".txt.gz";
                    
    ofstream decFile;
    decFile.open(fileName, std::ios_base::out | std::ios_base::binary);
  
    boost::iostreams::filtering_ostream jointStream; 
    jointStream.push(boost::iostreams::gzip_compressor());
    jointStream.push(decFile);

    double rateVal = mmsmc_ -> getParameterValue("theta"); //genome-wide avg

    for(size_t i = 0; i < fragLength; ++i) {
      
      double avgRate = 0.; //avg rate at site i
    
      for(size_t k = 0; k < numCatTheta; ++k) {
        
        double focalCatProb = rateFwdMat[i][k] * rateBackMat[i][k];
        avgRate += focalCatProb * rateVal * mmsmc_ -> getParameterScalings()[0] -> getCategories()[k];
      }
    
      jointStream << avgRate << endl;
    }
  
    boost::iostreams::close(jointStream);
  
    cout << " done." << endl;
  }
}

//NOTE for a single pair of genomes
void SmcDecodingWrapper::decodeBatchOfBreakpoints_(const vector< pair< size_t, size_t > >& batchBps, const vector< size_t >& blockIdVector) {
  
  size_t fragSize = smcOptions_ -> getFragmentSize();
  size_t numHs = mmsmc_ -> getNumberOfHiddenStates();
  
  auto decodeParallel = [&] (size_t thread_id) {
    
    auto focalPair = batchBps[thread_id];
    
    size_t breakpointStart = get< 0 >(focalPair);
    size_t breakpointEnd = get< 1 >(focalPair);
    
    size_t breakpointLength = breakpointEnd - breakpointStart;
    size_t numFragments = (breakpointLength + fragSize - 1) / fragSize; 
    
    size_t blockId = blockIdVector[thread_id] + 1; //indexes from 1
    
    //init probs
    Vdouble pi(numHs);
    
    for(size_t i = 0; i < numHs; ++i) {
      pi[i] = 1. / static_cast< double >(numHs);
    }
    
    size_t fileIndex = 1; //simple counter because breakpoints can heve different lengths
    
    for(size_t i = 0; i < numFragments; ++i) {

      size_t start = i * fragSize + breakpointStart;
      size_t end = start + fragSize;
            
      //in case breakpointLength is not a multiple of fragSize
      if(end > breakpointEnd) {
        end = breakpointEnd;
      }
      
      string fileName = smcOptions_ -> getLabel() + "_diploid_1_block_" + TextTools::toString(blockId) + 
                        "_file_index_" + TextTools::toString(fileIndex);
      
      bool restricted = smcOptions_ -> restrictedModulation(); //time restricted or not
      zipHMM::Matrix postProbMatrix;
      
      mPsmc_ -> getWholePsmcVector()[0] -> posteriorDecodingUsingZipHMM(start, end, fileName, restricted, postProbMatrix, pi);

      //fills in pi for next fragment
      for(size_t l = 0; l < numHs; ++l) {
        pi[l] = postProbMatrix(l, postProbMatrix.get_width() - 1);
      }

      cout << "   done for " << mPsmc_ -> getWholePsmcVector()[0] -> getBiHaploidName();
      cout << " in fragment: " << start + 1 << " to " << end << endl;
      
      ++fileIndex;
    }
  };

  
  size_t numberOfTasks = batchBps.size();
  
  thread* threadVector = new thread[numberOfTasks];
    
  for(size_t i = 0; i < numberOfTasks; ++i) {
    threadVector[i] = thread(decodeParallel, i);
  }
  
  for(size_t i = 0; i < numberOfTasks; ++i) {
    threadVector[i].join();
  }
  
  delete [] threadVector;
 
}

//NOTE for multiple pairs of genomes
void SmcDecodingWrapper::decodeBatchOfBreakpoints_(size_t numAvailThreads, bool missing,
                                                   const vector< pair< size_t, size_t > >& batchBreakpoints,
                                                   const vector< size_t >& blockIdVector) {

  size_t numDiploids = mPsmc_ -> getWholePsmcVector().size();
  size_t numHs = mmsmc_ -> getNumberOfHiddenStates();
  size_t fragSize = smcOptions_ -> getFragmentSize();
    
  size_t numberOfTasks = batchBreakpoints.size();
    
  auto decodeBreakpoint = [&] (size_t thread_id) {
      
    //init. probs. for the Markov-modulated SMC for each diploid
    VVdouble piMmSmc(numDiploids, Vdouble(numHs));
  
    for(size_t i = 0; i < numDiploids; ++i) { 
      
      for(size_t j = 0; j < numHs; ++j) {
        piMmSmc[i][j] = 1. / static_cast< double >(numHs);
      }
    }
      
    auto focalPair = batchBreakpoints[thread_id];
    
    size_t breakpointStart = get< 0 >(focalPair);
    size_t breakpointEnd = get< 1 >(focalPair);
    
    size_t breakpointLength = breakpointEnd - breakpointStart;
    size_t numFragments = (breakpointLength + fragSize - 1) / fragSize; 
    
    size_t blockId = blockIdVector[thread_id] + 1;
    
    for(size_t i = 0; i < numFragments; ++i) {  
      
      size_t startPos = i * fragSize + breakpointStart;
      size_t endPos = startPos + fragSize;
    
      if(endPos > breakpointEnd) { 
        endPos = breakpointEnd;
      }
    
      jointlyDecodeFragment_(piMmSmc, startPos, endPos, blockId, i, numAvailThreads, missing);

      cout << "Joint-decoding done for fragment: " << startPos + 1 << " to " << endPos << endl << endl;
    }
    
  };
  
  //
  
  if(smcOptions_ -> decodeBreakpointsParallel()) {

    thread* threadVector = new thread[numberOfTasks];
    
    for(size_t i = 0; i < numberOfTasks; ++i) {
      threadVector[i] = thread(decodeBreakpoint, i);
    }
  
    for(size_t i = 0; i < numberOfTasks; ++i) {
      threadVector[i].join();
    }
  
    delete [] threadVector;
  }
  
  else {
      
    for(size_t i = 0; i < numberOfTasks; ++i) {
      decodeBreakpoint(i);
    }
  }

}

void SmcDecodingWrapper::computeAverageRateOverDiploids_(const string& rate) {
    
  cout << "Writing files with average " << rate << " over diploids..."; cout.flush();
  
  size_t numDiploids = mPsmc_ -> getWholePsmcVector().size();
  size_t fragSize = smcOptions_ -> getFragmentSize();
  size_t numBreakpoints = mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints().size();
  
  //read files and compute the sum of rate for every position
  for(size_t i = 0; i < numBreakpoints; ++i) {
    
    auto focalPair = mPsmc_ -> getWholePsmcVector()[0] -> getSequenceBreakpoints()[i];

    size_t breakpointStart = get< 0 >(focalPair);
    size_t breakpointEnd = get< 1 >(focalPair);
    
    size_t breakpointLength = breakpointEnd - breakpointStart;
    size_t numFragments = (breakpointLength + fragSize - 1) / fragSize;
    
    for(size_t j = 0; j < numFragments; ++j) {
        
      VVdouble diploidLandscapes(numDiploids, Vdouble(0)); //diploid -> site
      
      for(size_t k = 0; k < numDiploids; ++k) {
          
        string file = smcOptions_ -> getLabel() + "_diploid_" + TextTools::toString(k + 1) + "_block_" +
                      TextTools::toString(i + 1) + "_file_index_" + TextTools::toString(j + 1) + "_estimated_" + rate + ".txt.gz";//+ "_Tmrca." + rate + ".txt.gz"; WARNING for other decoding
        
        if(smcOptions_ -> restrictedModulation()) {
          file = smcOptions_ -> getLabel() + "_diploid_" + TextTools::toString(k + 1) + "_block_" +
                 TextTools::toString(i + 1) + "_file_index_" + TextTools::toString(j + 1) + "_Tmrca." + rate + ".txt.gz";
        }

        diploidLandscapes[k] = readLandscapeFromFile(file);
      }
      
      //for every position, computes average over diploids, for fragment j
      size_t numPositions = diploidLandscapes.front().size();
      Vdouble averageRate(numPositions);
    
      for(size_t l = 0; l < numPositions; ++l) {
        
        Vdouble ratePos(numDiploids);
        
        for(size_t d = 0; d < numDiploids; ++d) {
          ratePos[d] = diploidLandscapes[d][l]; 
        }
        
        averageRate[l] = accumulate(begin(ratePos), end(ratePos), 0.) / static_cast< double >(numDiploids);
      }
      
      //prints to file
      string fileName = smcOptions_ -> getLabel() + "_average_" + rate + "_block_" + TextTools::toString(i + 1) +
                        "_file_index_" + TextTools::toString(j + 1) + ".txt";
    
      writeSingleNucleotideLandscapeToFile(averageRate, fileName);
    }
  }
  
  cout << " done." << endl; 
}

