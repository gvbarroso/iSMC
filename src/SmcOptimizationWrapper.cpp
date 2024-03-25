/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 13/04/2018
 * Last modified: 15/12/2022
 *
 */


#include <iostream>

#include "SmcOptimizationWrapper.h"
#include "BackupListenerOv.h"
#include "Global.h"

using namespace std;
using namespace bpp;


shared_ptr< SplinesModel > SmcOptimizationWrapper::selectBestModel() { 
  size_t bestModelIndex = 0;   
  double bestAic = listOfModels_[0] -> getAic();
  for(size_t i = 1; i < listOfModels_.size(); ++i) {
    if(listOfModels_[i] -> getAic() < bestAic) {
      bestAic = listOfModels_[i] -> getAic();
      bestModelIndex = i;
    }
  }
  return listOfModels_[bestModelIndex];
}
  
void SmcOptimizationWrapper::stepwiseExpectationMaximization() { 

  /*
   * NOTE Draft
   *
   * Will be done fragment-wise to save memory
   * Challenge: locus-skipping EM of Song does not help when the number of hidden states is large
   *

  BaumWelch baumWelcher;
  
  //using the full Markov-modulated SMC
  size_t numDiploids = mPsmc_ -> getWholePsmcVector().size();
  size_t numHs = mmsmc_ -> getNumberOfHiddenStates(); 

  VVVdouble indvTransMats //Baum-Welch proposed TM's for each diploid

  //for each diploid i

    VVdouble newTransMat;

    //for the entire sequence, in fragments 

      //select fragment

      //forwar-backward using zipHMM functions
      vector< double > fwdScales(0);
      zipHM::Matrix fwdMatrix;
      zipHM::Matrix bckwdMatrix;

      zipHMM::forward(fwdMatrix, fwdScales ...);
      zipHMM::backward(bckwdMatrix, fwdScales ...);

      baumWelcher.incrementTransitionMatrix(newTransMat ...);
  
  //normalise newTransMat built from whole seq

  indvTransMats[i] = newTransMat
  
  //compute averge of all individual TM's

  //use MatrixOptimisationFunction

  */

}
  
void SmcOptimizationWrapper::optimizeParameters() { 

  ParameterList params(bestParameters_); 
  params.deleteParameter("theta"); //we fix theta to the Wattersons's estimator
  
  //coalescence rates will be optimised with splines:
  for(size_t i = 0; i < mmsmc_ -> getNumberOfIntervals(); ++i) {
    params.deleteParameter("l" + TextTools::toString(i));
  }
  
  createAndFitSplinesModels_(params);
  //updates params
  shared_ptr< SplinesModel > bestSplines = selectBestModel();
    
  bestSplines -> mapLambdasFromSplines(mmsmc_ -> getLambdaVector());

  params.matchParametersValues(mmsmc_ -> getParameters());   
  params.matchParametersValues(mmsmc_ -> getLambdaVector());

  for(size_t i = 0; i < mmsmc_ -> getParameterScalings().size(); ++i) {
    params.matchParametersValues(mmsmc_ -> getParameterScalings()[i] -> getIndependentParameters());
  }
  
  for(size_t i = 0; i < mmsmc_ -> getParameterTransitions().size(); ++i) {
    params.matchParametersValues(mmsmc_ -> getParameterTransitions()[i] -> getParameters());
  }
  
  fireUpdateBestValues_(bestSplines.get(), params);
}

void SmcOptimizationWrapper::optimizeParameters(const ParameterList& backupParams) { 

  cout << endl << "Resuming optimisation." << endl;
   
  ParameterList params;
  params.addParameters(bestParameters_); //non-splines params w/ constraints
  
  params.deleteParameter("theta");
  //coalescence rates will be optimised with splines:
  for(size_t i = 0; i < mmsmc_ -> getNumberOfIntervals(); ++i) 
    params.deleteParameter("l" + TextTools::toString(i));
  
  //adds splines params to match partially optimised values 
  for(size_t i = 0; i < backupParams.size(); ++i) {
    
    string candidateSplineHeight = "y" + TextTools::toString(i);
    string candidateSplineDeriv = "y" + TextTools::toString(i) + "_prime";
    
    if(backupParams.hasParameter(candidateSplineHeight)) { //if knot i existed before
        
      params.addParameter(new Parameter(candidateSplineHeight, 1., Parameter::R_PLUS_STAR));
      params.addParameter(new Parameter(candidateSplineDeriv, 0.,
                                      make_shared<IntervalConstraint>(-1., 1., true, true)));
    }
  }
  params.matchParametersValues(backupParams);
  
  params.printParameters(cout);
  
  
  shared_ptr< SplinesModel > smf(new SplinesModel(mmsmc_, mmsmcep_, mmsmctp_, mPsmc_, params,
                                                  smcOptions_ -> getIgnoredParameters(),
                                                  smcOptions_ -> getInitNumberOfKnots(),
                                                  smcOptions_ -> getSplinesTypeOption()));

  fitModel_(smf);
  listOfModels_.push_back(smf);

  //updates
  shared_ptr< SplinesModel > bestSplines = selectBestModel();
  bestSplines -> mapLambdasFromSplines(mmsmc_ -> getLambdaVector());
  params.matchParametersValues(mmsmc_ -> getParameters());   

  for(size_t i = 0; i < mmsmc_ -> getParameterScalings().size(); ++i) {
    params.matchParametersValues(mmsmc_ -> getParameterScalings()[i] -> getIndependentParameters());
  }
  
  for(size_t i = 0; i < mmsmc_ -> getParameterTransitions().size(); ++i) {
    params.matchParametersValues(mmsmc_ -> getParameterTransitions()[i] -> getParameters());
  }

  fireUpdateBestValues_(bestSplines.get(), params);
}

void SmcOptimizationWrapper::writeEstimatesToFile(const ParameterList& params, double ll) {
    
  //writes params in a simple txt format that is easy to parse
  ofstream parameterEstimates;
  
  parameterEstimates.open(smcOptions_ -> getLabel() + "_estimates.txt");
  parameterEstimates << setprecision(20) << "LogLikelihood = " << ll << endl << endl;
  
  vector< string > parameterNames = params.getParameterNames();
  
  for(size_t i = 0; i < params.size(); ++i) {
      
    parameterEstimates << parameterNames[i];
    parameterEstimates << " ";
    
    parameterEstimates << setprecision(20);
    
    parameterEstimates << params.getParameterValue(parameterNames[i]);
    parameterEstimates << endl;
  }
  
  parameterEstimates.close();    
}

void SmcOptimizationWrapper::writeDemographyToFile() {
    
  Vdouble timeIntervals = mmsmc_ -> getTimeIntervals();
  
  ofstream demoHistory;
  demoHistory.open(smcOptions_ -> getLabel() + "_demography.txt");
  demoHistory << "lower_bound" << "\t" << "upper_bound" << "\t" << "scaled_coalescence_rate" << endl;

  for(size_t i = 0; i < mmsmc_ -> getNumberOfIntervals() - 1; ++i) {

    demoHistory << setprecision(6);

    demoHistory << timeIntervals[i];
    demoHistory << "\t";
    
    demoHistory << timeIntervals[i + 1];
    demoHistory << "\t";
    
    demoHistory << bestParameters_.getParameterValue("l" + TextTools::toString(i));
    demoHistory << endl;
  }
  
  demoHistory.close();    
}

void SmcOptimizationWrapper::fireUpdateBestValues_(SplinesModel* bestSm, const ParameterList& params) {

  if(bestSm -> getAic() < bestAic_) {
      
    bestAic_ = bestSm -> getAic();
    
    bestParameters_.matchParametersValues(params); //non-lambdas
    bestParameters_.matchParametersValues(mmsmc_ -> getLambdaVector()); //lambdas
    
    writeDemographyToFile();
    
    ParameterList estimParams = bestSm -> getParameters();
    estimParams.addParameter(mmsmc_ -> parameter("theta"));
    
    writeEstimatesToFile(estimParams, bestSm -> getLogLikelihood());
  }
  
  else {
    cout << "Warning!! Optimisation did not reduce AIC." << endl;
  }
}  

void SmcOptimizationWrapper::createAndFitSplinesModels_(ParameterList& params) {

  //creates and fit models from minNumKnots to maxNumKnots  
  unsigned int minNumKnots = smcOptions_ -> getInitNumberOfKnots(); 
  unsigned int maxNumKnots = smcOptions_ -> getMaxNumberOfKnots();
  
  string splinesType = smcOptions_ -> getSplinesTypeOption();
  
  for(size_t i = minNumKnots; i <= maxNumKnots; ++i) {

    cout << endl << "Optimising with " << i << " splines knot(s)..." << endl;

    shared_ptr< SplinesModel > smf(new SplinesModel(mmsmc_, mmsmcep_, mmsmctp_, mPsmc_,
                                                    params, i, splinesType));

    fitModel_(smf);
    listOfModels_.push_back(smf);
  }
}
    
void SmcOptimizationWrapper::fitModel_(shared_ptr<SplinesModel> smf) { 
  
  std::cout << "\nOptimizing the following parameters:\n";
  smf -> fetchModelParameters().printParameters(cout);
  //smf potentially has both splines parameters and spatial rates parameters
  auto rfw = make_shared<ReparametrizationFunctionWrapper>(smf, smf -> fetchModelParameters()); //reparametrization of all params

  auto tpnd = make_shared<ThreePointsNumericalDerivative>(smf);
   
  unique_ptr< OptimizerInterface > chosenOptimizer;

  if(smcOptions_ -> getOptimizerOption() == "Powell") {
    chosenOptimizer.reset(new PowellMultiDimensions(rfw));
  }

  else if(smcOptions_ -> getOptimizerOption() == "NewtonRhapson") {
      
    if(!smcOptions_ -> enforceFlatDemo()) {
      tpnd -> setParametersToDerivate(smf -> fetchModelParameters().getParameterNames());
    }
    
    else { //if a flat demography is specified in the params file, we don't derivate splines
      tpnd -> setParametersToDerivate(smf -> fetchNonSplinesParameters().getParameterNames());
    }
    
    tpnd -> enableFirstOrderDerivatives(true);
    tpnd -> enableSecondOrderDerivatives(true);
    tpnd -> enableSecondOrderCrossDerivatives(false); //Pseudo-Newton
    
    chosenOptimizer.reset(new PseudoNewtonOptimizer(tpnd)); 
  }
  
  else {
    throw Exception("iSMC::Mis-specified numerical_optimizer!");
  }
  
  string optimProfile = smcOptions_ -> getLabel() + "_" + TextTools::toString(smf -> getNumberOfKnots()) + "_knots_profile.txt";
  auto profiler = make_shared<StlOutputStream>(make_unique<ofstream>(optimProfile, ios::out));
  chosenOptimizer -> setProfiler(profiler);
  
  string optimMsgs = smcOptions_ -> getLabel() + "_" + TextTools::toString(smf -> getNumberOfKnots()) + "_knots_messages.txt";
  auto messenger = make_shared<StlOutputStream>(make_unique<ofstream>(optimMsgs, ios::out));
  chosenOptimizer -> setMessageHandler(messenger);
    
  if(smcOptions_ -> getOptimizerOption() == "Powell") {
      
    //usual case: we fit the splines models to capture demography
    if(!smcOptions_ -> enforceFlatDemo()) {
      chosenOptimizer -> init(rfw->getParameters());
    }
    
    else { //to test the effect of assuming a flat demography
      ReparametrizationFunctionWrapper tmp(smf, smf -> fetchNonSplinesParameters()); //reparam. of non-splines params
      chosenOptimizer -> init(tmp.getParameters());  
    }
  }
  
  else if(smcOptions_ -> getOptimizerOption() == "NewtonRhapson") {
      
    chosenOptimizer -> setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    
    //usual case: we fit the splines models to capture demography
    if(!smcOptions_ -> enforceFlatDemo()) {
      chosenOptimizer -> init(smf -> fetchModelParameters());
    }
    
    else { //to test the effect of assuming a flat demography
      chosenOptimizer -> init(smf -> fetchNonSplinesParameters());  
    }
    
  }
  
  shared_ptr< FunctionStopCondition > stopCond;
  
  if(smcOptions_ -> relativeStopCond()) {
    stopCond.reset(new RelativeStopCondition(chosenOptimizer.get(), smcOptions_ -> getFunctionTolerance(), smcOptions_ -> getParametersTolerance())); 
  }
  
  else {
    stopCond.reset(new FunctionStopCondition(chosenOptimizer.get(), smcOptions_ -> getFunctionTolerance()));
  }
  
  chosenOptimizer -> setStopCondition(stopCond);
  
  //handling optimisation issues, e.g. Powell: line minimization failing 
  try {
    auto blo = make_shared<BackupListenerOv>(smcOptions_ -> getLabel() + "_backup_params.txt");
    chosenOptimizer -> addOptimizationListener(blo);
    chosenOptimizer -> optimize();
  }
  
  catch(bpp::Exception& e) {
    cout << endl << "Warning!! Error during optimization, convergence might not be reached." << endl;
    cout << "iSMC will proceed, but check log file." << endl;
  }
  
  smf -> computeAic();

  cout << endl << endl << "logLikelihood = " << setprecision(6) << smf -> getLogLikelihood() << endl;

  ParameterList optimisedParams(smf -> getParameters());
  cout << endl << "Optimized iSMC parameters:" << endl;
  optimisedParams.printParameters(cout);
}


