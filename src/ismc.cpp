/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 16/04/2018
 * Last modified: 26/06/2019
 * Source code for iSMC.
 * Our model is based on the Sequentially Markov Coalescent (SMC) from Mcvean and Cardin (2005).
 * We use the transition probabilities as described for the MSMC (Schiffels and Durbin, 2014).
 * 
 * Tools for perfoming simulations, plotting and converting files are provided as R scripts.
 * 
 */

////////////////////////////////////////////////////////////////////////////////////////////////////


#include <limits>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <string>
#include <cstring>
#include <algorithm> 
#include <cmath>
#include <thread>
#include <chrono>
#include <iterator>
#include <cstdlib>
#include <utility>

#include <boost/algorithm/string.hpp>
#include <boost/bimap.hpp>

#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/NumTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/BppApplication.h>

#include "Global.h"
#include "HmmStatesLibrary.h"
#include "OptionsContainer.h"
#include "PolymorphismData.h"
#include "GammaWithHotspots.h"
#include "ParameterCategoryTransitions.h"
#include "MarkovModulatedSmc.h"
#include "MmSmcEmissionProbabilities.h"
#include "MmSmcTransitionProbabilities.h"
#include "MmPsmc.h"
#include "MultipleMmPsmc.h"
#include "SmcOptimizationWrapper.h"
#include "SmcDecodingWrapper.h"
#include "Vcf.h"


using namespace bpp;
using namespace std;


int main(int argc, char *argv[]) { 
  
  cout << endl;
  cout << "******************************************************************" << endl;
  cout << "*                    iSMC, version 0.0.23                        *" << endl;
  cout << "*                                                                *" << endl;
  cout << "*                                                                *" << endl;
  cout << "*            Recombination                                       *" << endl;
  cout << "*            A mosaic in the genome                              *" << endl;
  cout << "*            Endless ancestors                                   *" << endl;
  cout << "*                                                                *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: G. Barroso                    Last Modif. 16/Dec/2021 *" << endl;
  cout << "*          J. Dutheil                                            *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if(argc == 1) {
    cout << "To use iSMC, please fill in the params file and simply call it from the command line: ismc params=[params_file].bpp" << endl;
    cout << "For more information, please email gvbarroso@gmail.com" << endl;
    return(0);
  }  

  try {

  ///////////////////////////////////////////////////////////////////////////////////////////
  //Reads params file
  BppApplication iSMC(argc, argv, "iSMC");
  iSMC.startTimer();
  map< string, string > params = iSMC.getParams();
  shared_ptr< OptionsContainer > smcOptions = make_shared< OptionsContainer >(params);  
  
  NUMBER_OF_AVAILABLE_THREADS = smcOptions -> getNumberOfThreads();
    
  //Handles input sequence data 
  shared_ptr< PolymorphismData > dataSet = make_shared< PolymorphismData >(smcOptions);  
  
  if(smcOptions -> printSeqs()) {
    dataSet -> printSequencesToFile();
  }
  ///////////////////////////////////////////////////////////////////////////////////////////  

  
  /////////////////////////////////////////////////////////////////////////////////////////// 
  cout << "Creating iSMC objects..."; cout.flush();
  
  //Different indices than in hidden states combinations 
  ParameterAlphabet parameterAlphabet; //defined in HmmStatesLibrary.h
  parameterAlphabet.insert({0, "theta"});
  parameterAlphabet.insert({1, "rho"});
  parameterAlphabet.insert({2, "ne"});
  parameterAlphabet.insert({3, "timeInterval"}); //in HS combinations (HS library), time has index 0
  
  //the spatial parameters: 
  double initalShape = 1.;
  vector< shared_ptr< bpp::DiscreteDistribution > > paramScalings;

  shared_ptr< bpp::DiscreteDistribution > thetaScaling; //theta
  if(smcOptions -> getThetaVarModel() == "Hotspot") {
      
    Vdouble categoryValues = { 1., 100. };
    map< size_t, Vdouble > categoryRanges;
    Vdouble backgroundRange = { 1., 1. }; 
    Vdouble hotspotRange = { 1., 1e+6 }; //do not put std::limits double as the upper limit
    categoryRanges.insert(make_pair(1, backgroundRange)); //1st position in categoryValues
    categoryRanges.insert(make_pair(2, hotspotRange)); //2nd position in categoryValues
    Vdouble categoryProbs = { 0.5, 0.5 }; //not parameters of iSMC
    thetaScaling = make_shared< SimpleDiscreteDistribution >(categoryValues, categoryRanges, categoryProbs);
    thetaScaling -> setNamespace("theta.");
  }
  
  else if(smcOptions -> getThetaVarModel() == "Gamma+Hotspot") {
      
    //lower bound of hotspot intensity is set to highest scaling factor from the discretized gamma
    GammaDiscreteDistribution toyGamma(smcOptions -> getNumberOfThetaCateg(), 0.05, 0.05);
    toyGamma.discretize();
    double lowerBoundHeat = toyGamma.getCategories()[smcOptions -> getNumberOfThetaCateg() - 1];
    double initialHeat = 10. * lowerBoundHeat;
    Vdouble hotspotRange = { lowerBoundHeat, 1e+6 }; //do not put std::limits double as the upper limit
    thetaScaling = make_shared< GammaWithHotspots >(smcOptions -> getNumberOfThetaCateg(), initalShape, initalShape, initialHeat, hotspotRange);
    thetaScaling -> setNamespace("theta.");
    thetaScaling -> aliasParameters("alpha", "beta");
    thetaScaling -> discretize();
  }
  
  else if(smcOptions -> getThetaVarModel() == "Gamma") {
      
    thetaScaling = make_shared< GammaDiscreteDistribution >(smcOptions -> getNumberOfThetaCateg(), initalShape, initalShape);
    thetaScaling -> setNamespace("theta.");
    thetaScaling -> aliasParameters("alpha", "beta");
    thetaScaling -> discretize();
  }
  
  else {
    throw Exception("Mis-specified heterogeneity model for Theta!");
  }
  paramScalings.push_back(thetaScaling); //theta occupies position 0 in the vector 
  
  shared_ptr< bpp::DiscreteDistribution > rhoScaling; //rho 
  if(smcOptions -> getRhoVarModel() == "Hotspot") {
      
    Vdouble categoryValues = { 1., 100. };
    map< size_t, Vdouble > categoryRanges;
    Vdouble backgroundRange = { 1., 1. }; 
    Vdouble hotspotRange = { 1., 1e+6 }; //do not put std::limits double as the upper limit
    categoryRanges.insert(make_pair(1, backgroundRange)); //1st position in categoryValues
    categoryRanges.insert(make_pair(2, hotspotRange)); //2nd position in categoryValues
    Vdouble categoryProbs = { 0.5, 0.5 }; //not parameters of iSMC
    rhoScaling = make_shared< SimpleDiscreteDistribution >(categoryValues, categoryRanges, categoryProbs);
    rhoScaling -> setNamespace("rho.");
  }
  else if(smcOptions -> getRhoVarModel() == "Gamma+Hotspot") {
      
    //lower bound of hotspot intensity is set to highest scaling factor from the discretized gamma
    GammaDiscreteDistribution toyGamma(smcOptions -> getNumberOfRhoCateg(), 0.05, 0.05);
    toyGamma.discretize();
    
    double lowerBoundHeat = toyGamma.getCategories()[smcOptions -> getNumberOfRhoCateg() - 1];
    
    double initialHeat = 10. * lowerBoundHeat;
    
    Vdouble hotspotRange = { lowerBoundHeat, 1e+6 }; //do not put std::limits double as the upper limit
    
    rhoScaling = make_shared< GammaWithHotspots >(smcOptions -> getNumberOfRhoCateg(), initalShape, initalShape, initialHeat, hotspotRange);
    rhoScaling -> setNamespace("rho.");
    rhoScaling -> aliasParameters("alpha", "beta");
    rhoScaling -> discretize();
  }
  
  else if(smcOptions -> getRhoVarModel() == "Gamma") {
    rhoScaling = make_shared< GammaDiscreteDistribution >(smcOptions -> getNumberOfRhoCateg(), initalShape, initalShape);
    rhoScaling -> setNamespace("rho.");
    rhoScaling -> aliasParameters("alpha", "beta");
    rhoScaling -> discretize();
  }
  
  else {
    throw Exception("Mis-specified heterogeneity model for Rho!");
  }
  
  paramScalings.push_back(rhoScaling); //rho occupies position 1 in the vector 
 
  shared_ptr< bpp::DiscreteDistribution > neScaling; //Ne
  if(smcOptions -> getNeVarModel() == "Hotspot") {
      
    Vdouble categoryValues = { 1., 100. };
    map< size_t, Vdouble > categoryRanges;
    Vdouble backgroundRange = { 1., 1. }; 
    Vdouble hotspotRange = { 1., 1e+6 }; //do not put std::limits double as the upper limit
    categoryRanges.insert(make_pair(1, backgroundRange)); //1st position in categoryValues
    categoryRanges.insert(make_pair(2, hotspotRange)); //2nd position in categoryValues
    Vdouble categoryProbs = { 0.5, 0.5 }; //not parameters of iSMC
    neScaling = make_shared< SimpleDiscreteDistribution >(categoryValues, categoryRanges, categoryProbs);
    neScaling -> setNamespace("ne.");
  }
  
  else if(smcOptions -> getNeVarModel() == "Gamma+Hotspot") {
      
    //lower bound of hotspot intensity is set to highest scaling factor from the discretized gamma
    GammaDiscreteDistribution toyGamma(smcOptions -> getNumberOfNeCateg(), 0.05, 0.05);
    toyGamma.discretize();
    double lowerBoundHeat = toyGamma.getCategories()[smcOptions -> getNumberOfNeCateg() - 1];
    double initialHeat = 10. * lowerBoundHeat;
    Vdouble hotspotRange = { lowerBoundHeat, 1e+6 }; //do not put std::limits double as the upper limit
    neScaling = make_shared< GammaWithHotspots >(smcOptions -> getNumberOfNeCateg(), initalShape, initalShape, initialHeat, hotspotRange);
    neScaling -> setNamespace("ne.");
    neScaling -> aliasParameters("alpha", "beta");
    neScaling -> discretize();
  }
  
  else if(smcOptions -> getNeVarModel() == "Gamma") {
      
    neScaling = make_shared< GammaDiscreteDistribution >(smcOptions -> getNumberOfNeCateg(), initalShape, initalShape);  
    neScaling -> setNamespace("ne.");
    neScaling -> aliasParameters("alpha", "beta");
    neScaling -> discretize();
  }
  
  else {
    throw Exception("Mis-specified heterogeneity model for Ne!");
  }
  
  paramScalings.push_back(neScaling); //ne occupies position 2 in the vector
    
  //transitions between parameter categories
  vector< shared_ptr< ParameterCategoryTransitions > > categoryTransitions;

  //theta
  shared_ptr< ParameterCategoryTransitions > thetaTrans = make_shared< ParameterCategoryTransitions >(smcOptions -> getNumberOfThetaCateg(),
                                                                                                      smcOptions -> getThetaVarModel(), "t");
  categoryTransitions.push_back(thetaTrans); //0
  
  //rho
  shared_ptr< ParameterCategoryTransitions > rhoTrans = make_shared< ParameterCategoryTransitions >(smcOptions -> getNumberOfRhoCateg(),
                                                                                                    smcOptions -> getRhoVarModel(), "r");
  categoryTransitions.push_back(rhoTrans); //1
  
  //ne 
  shared_ptr< ParameterCategoryTransitions > neTrans = make_shared< ParameterCategoryTransitions >(smcOptions -> getNumberOfNeCateg(),
                                                                                                   smcOptions -> getNeVarModel(), "ne");
  categoryTransitions.push_back(neTrans); //2
  ///////////////////////////////////////////////////////////////////////////////////////////

  
  ///////////////////////////////////////////////////////////////////////////////////////////
  //Key pointers that are shared among all bi-haploid samples
  shared_ptr< HmmStatesLibrary > hmmLib;
  shared_ptr< MarkovModulatedSmc > model;
  shared_ptr< MmSmcTransitionProbabilities > transitions;
  shared_ptr< MmSmcEmissionProbabilities > emissions;
  shared_ptr< MultipleMmPsmc > multiMmPsmc;
  
  cout << " done." << endl;
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  
  //////////////////////////////////////////////////////////////////////////////////////////
  //Optimisation:
  if(smcOptions -> optimize()) {
    
    cout << "Setting up optimisation." << endl;
    
    unsigned int numIntervals = smcOptions -> getNumberOfIntervals();  
    hmmLib = make_shared< HmmStatesLibrary >(numIntervals, paramScalings, parameterAlphabet);

    //Builds MultiplePsmc object for optimisation:
    auto psmcVectorOptim = vector< shared_ptr< MmPsmc > >();
    psmcVectorOptim.reserve(dataSet -> getNumberOfDiploids());
    
    model = make_shared< MarkovModulatedSmc >(smcOptions, dataSet -> getSnpCalling(), hmmLib, paramScalings, categoryTransitions, parameterAlphabet);
    
    transitions = make_shared< MmSmcTransitionProbabilities >(model);
    emissions = make_shared< MmSmcEmissionProbabilities >(model, dataSet -> getSnpCalling());
        
    for(size_t i = 0; i < dataSet -> getNumberOfDiploids(); ++i) {
      
      shared_ptr< MmPsmc > biHaploid = make_shared< MmPsmc >(dataSet -> getNames()[i],
                                                             dataSet -> getSnpCalling()[i],
                                                             dataSet -> getBreakpoints(),
                                                             smcOptions -> getMissingBlocksLength(),
                                                             model, transitions, emissions);
      psmcVectorOptim.push_back(biHaploid);
    }
    
    multiMmPsmc = make_shared< MultipleMmPsmc >(psmcVectorOptim);

    SmcOptimizationWrapper smcWrapper(model, emissions, transitions, multiMmPsmc, smcOptions);
    
    //decides whether or not to look for backup parameters in the working dir
    if(smcOptions -> resumeOptim()) {
      ParameterList backupParams = MarkovModulatedSmc::readParametersFromFile(smcOptions -> getLabel() + "_backup_params.txt");
      smcWrapper.optimizeParameters(backupParams);
    }
    else {
      smcWrapper.optimizeParameters();
    }
    
    if(smcOptions -> computeCI()) {
        
      cout << endl << "Computing 95% confidence intervals of parameter estimates..." << endl;
      
      SplinesModel* optimizedSplines = smcWrapper.selectBestModel().get();
      ParameterList optimParams = optimizedSplines -> fetchModelParameters();
      vector< string > paramNames = optimParams.getParameterNames();

      ThreePointsNumericalDerivative* tpnd = new ThreePointsNumericalDerivative(optimizedSplines);
      tpnd -> enableFirstOrderDerivatives(true);
      tpnd -> enableSecondOrderDerivatives(true);
      
      if(smcOptions -> computeCovar()) {
        tpnd -> enableSecondOrderCrossDerivatives(true); 
      }
      
      else {
        tpnd -> enableSecondOrderCrossDerivatives(false);
      }
      
      for(size_t i = 0; i < paramNames.size(); ++i) { //loops over paramNames because optimParams potentially changes size 
          
        shared_ptr<Constraint> paramConstraint = optimParams.getParameter(paramNames[i]).getConstraint();
        
        //cf lines 203 / 204 of ThreePointsNumericalDerivative.cpp
        double h = (1. + abs(optimParams.getParameterValue(paramNames[i]))) * tpnd -> getInterval();
        double lowerPointDerivative = optimParams.getParameterValue(paramNames[i]) - h - tpnd -> getInterval() / 2.; //conservative
        double upperPointDerivative = optimParams.getParameterValue(paramNames[i]) + h + tpnd -> getInterval() / 2.; //conservative
        
        if(!paramConstraint -> includes(lowerPointDerivative, upperPointDerivative)) {
          optimParams.deleteParameter(paramNames[i]);
          cout << "   Numerical derivative can't be computed for " << paramNames[i] << "! Estimate is too close to boundary." << endl;
        }
      }
      
      //updates
      paramNames = optimParams.getParameterNames();
      tpnd -> setParametersToDerivate(paramNames);

      RowMatrix< double > varCovarMatrix;
      
      if(smcOptions -> computeCovar()) {
          
        RowMatrix< double >* hessian = NumTools::computeHessianMatrix(*tpnd, optimParams);
        MatrixTools::inv(*hessian, varCovarMatrix);
        //prints covariance matrix
        ofstream outMatrix("varCovarMatrix.txt", ios::out);
        outMatrix << paramNames[0];
        
        for(size_t j = 1; j < paramNames.size(); j++) {
          outMatrix << " " << paramNames[j];
        }
        outMatrix << endl;
        
        for(size_t i = 0; i < paramNames.size(); i++) {
            
          outMatrix << paramNames[i];
          
          for(size_t j = 0; j < paramNames.size(); j++) {
              
            outMatrix << " " << varCovarMatrix(i, j);
          }
          outMatrix << endl;
        }
        
        outMatrix.close();
      }
      
      else { //else the variance-covariance matrix is diagonal (ie only variances are computed)
          
        RowMatrix< double >* hessian = computeDiagonalHessian(*tpnd, optimParams);
        invertDiagonalMatrix(*hessian, varCovarMatrix);
      }
      
      //print 95% CI's
      ofstream outCI("ConfidenceIntervals.txt", ios::out);
      
      for(size_t i = 0; i < optimParams.size(); ++i) {
          
        double upper = optimParams.getParameterValue(paramNames[i]) + 1.96 * pow(varCovarMatrix(i, i), 0.5);
        double lower = optimParams.getParameterValue(paramNames[i]) - 1.96 * pow(varCovarMatrix(i, i), 0.5);
        outCI << paramNames[i] << " " << lower << " - " << upper << endl;
      }
      
      outCI.close(); 
    }
    
    psmcVectorOptim.clear();
  }
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  //Posterior decoding: ATM, using memory efficient version and focusing on spatial rate
  if(smcOptions -> decode()) {
    
    cout << endl << "Initiating Posterior Decoding." << endl;
    
    if(smcOptions -> decodeDiploidsParallel() && smcOptions -> decodeBreakpointsParallel()) {
      cout << "NOTE: parallelising over both diploids and breakpoints." << endl;
      cout << "Messages will be messy from now on!" << endl;
    }
    
    hmmLib = make_shared< HmmStatesLibrary >(smcOptions -> getNumberOfIntervals(), paramScalings, parameterAlphabet);

    //Builds MultiplePsmc object for decoding 
    auto psmcVectorDec = vector< shared_ptr< MmPsmc > >();
    
    psmcVectorDec.reserve(dataSet -> getNumberOfDiploids()); 
    
    model = make_shared< MarkovModulatedSmc >(smcOptions -> getNumberOfIntervals(),
                                              smcOptions -> getTimeDisc(),
                                              smcOptions -> getTmax(),
                                              dataSet -> getSnpCalling(),
                                              hmmLib, paramScalings, categoryTransitions, parameterAlphabet);
    
    transitions = make_shared< MmSmcTransitionProbabilities >(model);
    emissions = make_shared< MmSmcEmissionProbabilities >(model, dataSet -> getSnpCalling());
    
    for(size_t i = 0; i < dataSet -> getNumberOfDiploids(); ++i) {
    
      shared_ptr< MmPsmc > biHaploid = make_shared < MmPsmc >(dataSet -> getNames()[i], dataSet -> getSnpCalling()[i],
                                                              dataSet -> getBreakpoints(),
                                                              model, transitions, emissions);
      
      psmcVectorDec.push_back(biHaploid);
    }
    
    multiMmPsmc = make_shared< MultipleMmPsmc >(psmcVectorDec);
    
    //reads optimized parameters from file and usem them to update all objects 
    ParameterList optimizedParameters = MarkovModulatedSmc::readParametersFromFile(smcOptions -> getLabel() + "_estimates.txt");
    model -> matchParametersValues(optimizedParameters);
    
    //creates a Splines Model, to map lambdas from optimised splines parameters 
    Splines decSplines(model, smcOptions -> getInitNumberOfKnots(), smcOptions -> getSplinesTypeOption());
    decSplines.matchParametersValues(optimizedParameters);
    decSplines.mapLambdasFromSplines(model -> getLambdaVector());
    model -> computeAverageCoalescenceTime();
    
    ParameterList decodingParams(model -> getLambdaVector());
    decodingParams.addParameters(model -> getParameters());
    
    //updates all parametrized objects
    for(size_t i = 0; i < model -> getParameterScalings().size(); ++i) {
        
      if(optimizedParameters.getCommonParametersWith(model -> getParameterScalings()[i] -> getIndependentParameters()).size() > 0) {
        model -> getParameterScalings()[i] -> matchParametersValues(optimizedParameters);
        model -> getParameterScalings()[i] -> discretize();
        decodingParams.addParameters(model -> getParameterScalings()[i] -> getParameters());
      }
    }
    
    for(size_t j = 0; j < model -> getParameterTransitions().size(); ++j) {
        
      if(optimizedParameters.getCommonParametersWith(model -> getParameterTransitions()[j] -> getIndependentParameters()).size() > 0) {
        model -> getParameterTransitions()[j] -> matchParametersValues(optimizedParameters);
        model -> getParameterTransitions()[j] -> setUpCategoryTransitionMatrix();
        decodingParams.addParameters(model -> getParameterTransitions()[j] -> getParameters());
      }
    }
    
    //updates parametrized matrices
    transitions -> setUpExpectedMatrix();
    emissions -> setUpExpectedMatrix();
    
    SmcDecodingWrapper decoder(model, multiMmPsmc, smcOptions, transitions, emissions);
    
    //Individual decoding + average landscape + joint-decoding + Baum-Welch
    if(multiMmPsmc -> getWholePsmcVector().size() > 1) { 
      
      //size_t numObsStates = emissions -> getNumberOfObservedStates(); 
      bool missing = false; //to trigger filtering of missing sites in HMM layers: numObsStates > 2
     
      decoder.decodeDataset(missing);
      decoder.computeAverageLandscapes(); //sample_mean landscape
    }
    
    else { //will do it for all rates with numCategories > 1
      decoder.decodeDiploid();
    }
    
    decoder.writeDiploidLabelsToFile(); //to be used by mapper
  }
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  iSMC.done();  
  
  } catch(exception& e) {
    cout << "iSMC terminated because of an error." << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
