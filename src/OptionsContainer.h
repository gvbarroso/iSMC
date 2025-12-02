/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 09/02/2020
 *
 */


#ifndef _OPTIONSCONTAINER_H_
#define _OPTIONSCONTAINER_H_

#include <string>
#include <vector>
#include <limits>
#include <thread>

#include <Bpp/App/ApplicationTools.h>


class OptionsContainer {
private:
  
  std::string numericalOptimizer_; 
  std::string fullDataSetLabel_;
  std::string fileType_;
  std::string seqCompressionType_;
  std::string maskCompressionType_;

  std::string splinesType_;
  std::string thetaVarModel_;
  std::string rhoVarModel_;
  std::string neVarModel_;

  std::string timeDisc_; //log_even (Li & Durbin 2011) or quantiles (Schiffels & Durbin 2014)
  std::string decRate_;
  std::string sequenceFilePath_; 
  std::string maskFilePath_;
  std::string maskFileType_;

  //NOTE for VCF input format
  std::string tabFilePath_;
  
  //characters in mask file that represent position passing filtering criteria
  std::vector< std::string > callableCode_; //default = 1, P
  std::vector< std::string > ignoredParams_;
  
  unsigned int numberOfThetaCategories_;
  unsigned int numberOfRhoCategories_;
  unsigned int numberOfNeCategories_;
  double maxRhoValue_;
  unsigned int numberOfIntervals_;
  //number of time intervals for posterior decoding; diff. from optim. because of speed
  unsigned int numberOfDecodingIntervals_;
  unsigned int initNumberOfKnots_; 
  unsigned int maxNumberOfKnots_;
  //from which length of missing segments we start cutting them out in the optimisation
  unsigned int skipMissingBlocks_; //default = 0 means we don't cut them
  unsigned int fragmentSize_; //decode in fragments of this size
  
  size_t numberOfThreads_;
  
  //NOTE for FASTA (or "SNP", used in simulations) input format
  std::vector< unsigned int > breakpointStarts_;
  std::vector< unsigned int > breakpointEnds_;
  std::vector< unsigned int > diploidIndices_;
  
  bool printSeqs_;
  bool resumeOptim_;
  bool optimize_;
  //to fit a model without splines parameters, testing the effect of mis-specified coal rates.
  bool enforceFlatDemography_; 
  bool computeCI_;
  bool computeCovar_;
  bool decode_;
  bool decodeDiploidsParallel_;
  bool decodeBreakpointsParallel_;
  bool relativeTolerance_;
  bool restrictedModulation_; 
  
  //to decouple the numerical tolerance between likelihood and parameter values
  double functionTolerance_;
  double parametersTolerance_;
  double tMax_; //maximum time interval when timeDisc_ == "log_even"
  //for VCF files, a threshold to decide whether to consider coordinate as missing data
  double filterQualityVcf_; //see Vcf class
  
public:
  OptionsContainer(std::map< std::string, std::string > parameterOptions):
  numericalOptimizer_(bpp::ApplicationTools::getStringParameter("numerical_optimizer", parameterOptions, "Powell", "", true, 4)), 
  fullDataSetLabel_(bpp::ApplicationTools::getStringParameter("dataset_label", parameterOptions, "my_dataset", "", true, 4)),
  fileType_(bpp::ApplicationTools::getStringParameter("input_file_type", parameterOptions, "")),
  seqCompressionType_(bpp::ApplicationTools::getStringParameter("seq_compression_type", parameterOptions, "none", "", true, 4)),
  maskCompressionType_(bpp::ApplicationTools::getStringParameter("mask_compression_type", parameterOptions, "none", "", true, 4)),
  splinesType_(bpp::ApplicationTools::getStringParameter("splines_type", parameterOptions, "Sigmoidal", "", true, 4)),  
  thetaVarModel_(bpp::ApplicationTools::getStringParameter("theta_var_model", parameterOptions, "Gamma", "", true, 4)),
  rhoVarModel_(bpp::ApplicationTools::getStringParameter("rho_var_model",  parameterOptions, "Gamma", "", true, 4)),  
  neVarModel_(bpp::ApplicationTools::getStringParameter("ne_var_model", parameterOptions, "Gamma", "", true, 4)),
  timeDisc_(bpp::ApplicationTools::getStringParameter("time_disc", parameterOptions, "log_even", "", true, 4)),
  decRate_(bpp::ApplicationTools::getStringParameter("dec_rate", parameterOptions, "rho", "", true, 4)),
  sequenceFilePath_(bpp::ApplicationTools::getAFilePath("sequence_file_path", parameterOptions)),
  maskFilePath_(bpp::ApplicationTools::getAFilePath("mask_file_path", parameterOptions, false, false, "", false, "none", 4)),
  maskFileType_(bpp::ApplicationTools::getStringParameter("mask_file_type", parameterOptions, "none", "", true, 4)),
  tabFilePath_(bpp::ApplicationTools::getAFilePath("tab_file_path", parameterOptions, false, false, "", false, "none", 4)),
  callableCode_(bpp::ApplicationTools::getVectorParameter<std::string>("callable_code", parameterOptions, ',', "1,P", "", true, 4)),
  ignoredParams_(bpp::ApplicationTools::getVectorParameter<std::string>("ignore_params", parameterOptions, ',', "none", "", true, 4)),
  numberOfThetaCategories_(bpp::ApplicationTools::getParameter<unsigned int>("number_theta_categories", parameterOptions, 1, "", true, 4)),
  numberOfRhoCategories_(bpp::ApplicationTools::getParameter<unsigned int>("number_rho_categories", parameterOptions, 1, "", true, 4)),
  numberOfNeCategories_(bpp::ApplicationTools::getParameter<unsigned int>("number_ne_categories", parameterOptions, 1, "", true, 4)),
  maxRhoValue_(bpp::ApplicationTools::getDoubleParameter("max_rho_value", parameterOptions, 1, "", 100, 4)),
  numberOfIntervals_(bpp::ApplicationTools::getParameter<unsigned int>("number_intervals", parameterOptions, 40, "", true, 4) + 1), // +1 because timeInterval[0] = 0, present time
  numberOfDecodingIntervals_(bpp::ApplicationTools::getParameter<unsigned int>("number_intervals_decoding", parameterOptions, numberOfIntervals_, "", true, 4)),
  initNumberOfKnots_(bpp::ApplicationTools::getParameter<unsigned int>("init_number_knots", parameterOptions, 3, "", true, 4)),
  maxNumberOfKnots_(bpp::ApplicationTools::getParameter<unsigned int>("max_number_knots", parameterOptions, 3, "", true, 4)),
  skipMissingBlocks_(bpp::ApplicationTools::getParameter<unsigned int>("skip_missing_blocks", parameterOptions, 0, "", true, 4)),
  fragmentSize_(bpp::ApplicationTools::getParameter<unsigned int>("fragment_size", parameterOptions, 3000000)),
  numberOfThreads_(bpp::ApplicationTools::getParameter<unsigned int>("number_threads", parameterOptions,
                                                                     std::thread::hardware_concurrency(), "", true, 4)),
  breakpointStarts_(bpp::ApplicationTools::getVectorParameter<unsigned int>("breakpoint_starts", parameterOptions, ',', "0", "", true, 4)),
  breakpointEnds_(bpp::ApplicationTools::getVectorParameter<unsigned int>("breakpoint_ends", parameterOptions, ',', "0", "", true, 4)),
  diploidIndices_(bpp::ApplicationTools::getVectorParameter<unsigned int>("diploid_indices", parameterOptions, ',', "(0,1)", "", true, 1)),
  printSeqs_(bpp::ApplicationTools::getBooleanParameter("print_seqs", parameterOptions, false, "", true, 4)),
  resumeOptim_(bpp::ApplicationTools::getBooleanParameter("resume_optim", parameterOptions, false, "", true, 4)),
  optimize_(bpp::ApplicationTools::getBooleanParameter("optimize", parameterOptions, true)),
  enforceFlatDemography_(bpp::ApplicationTools::getBooleanParameter("enforce_flat_demo", parameterOptions, false, "", true, 4)),
  computeCI_(bpp::ApplicationTools::getBooleanParameter("compute_confidence_interval", parameterOptions, false, "", true, 4)),
  computeCovar_(bpp::ApplicationTools::getBooleanParameter("compute_covariance", parameterOptions, false, "", true, 4)),
  decode_(bpp::ApplicationTools::getBooleanParameter("decode", parameterOptions, true)),
  decodeDiploidsParallel_(bpp::ApplicationTools::getBooleanParameter("decode_diploids_parallel", parameterOptions, false, "", true, 4)),
  decodeBreakpointsParallel_(bpp::ApplicationTools::getBooleanParameter("decode_breakpoints_parallel", parameterOptions, false, "", true, 4)),
  relativeTolerance_(bpp::ApplicationTools::getBooleanParameter("relative_tolerance", parameterOptions, false, "", true, 4)),
  restrictedModulation_(bpp::ApplicationTools::getBooleanParameter("restricted_modulation", parameterOptions, false, "", true, 4)),
  functionTolerance_(bpp::ApplicationTools::getDoubleParameter("function_tolerance", parameterOptions, 1e-4)),
  parametersTolerance_(bpp::ApplicationTools::getDoubleParameter("parameters_tolerance", parameterOptions, std::numeric_limits< double >::max(), "", true, 4)),
  tMax_(bpp::ApplicationTools::getDoubleParameter("max_time", parameterOptions, 15., "", true, 4)),
  filterQualityVcf_(bpp::ApplicationTools::getDoubleParameter("vcf_quality_filter", parameterOptions, 0., "", true, 4))
  {
    if((diploidIndices_.size() % 2) != 0) {
      throw bpp::Exception("iSMC::Odd number of diploid indices!"); 
    }
  }
  
public:
  const std::string& getOptimizerOption() const { 
    return numericalOptimizer_;
  }

  const std::string& getLabel() const {
    return fullDataSetLabel_;
  }
  
  const std::string& getFileType() const {
    return fileType_;
  }
  
  const std::string& getSeqCompressionType() const {
    return seqCompressionType_;
  }

  const std::string& getMaskCompressionType() const {
    return maskCompressionType_;
  }
  
  const std::string& getSplinesTypeOption() const {
    return splinesType_;
  }
  
  const std::string& getThetaVarModel() const {
    return thetaVarModel_;
  }

  const std::string& getRhoVarModel() const {
    return rhoVarModel_;
  }

  const std::string& getNeVarModel() const {
    return neVarModel_;
  }
  
  const std::string& getSequenceFileName() const {
    return sequenceFilePath_;
  }
  
  const std::string& getMaskFileName() const {
    return maskFilePath_;
  }
  
  const std::string& getMaskFileType() const {
    return maskFileType_;
  }
  
  const std::string& getTabFileName() const {
    return tabFilePath_;
  }
  
  const std::string& getTimeDisc() const {
    return timeDisc_;
  }
  
  const std::string& getDecRate() const {
    return decRate_;
  }
  
  const std::vector< std::string >& getCallableCode() const {
    return callableCode_;
  }
  
  const std::vector< std::string >& getIgnoredParameters() const {
    return ignoredParams_;
  }
  
  unsigned int getNumberOfThetaCateg() const {
    return numberOfThetaCategories_;
  }
  
  unsigned int getNumberOfRhoCateg() const {
    return numberOfRhoCategories_;
  }
  
  unsigned int getNumberOfNeCateg() const {
    return numberOfNeCategories_;
  }
  
  unsigned int getNumberOfIntervals() const { 
    return numberOfIntervals_;
  }
  
  double getMaxRhoValue() const {
    return maxRhoValue_;
  }
  
  unsigned int getNumberOfDecodingIntervals() const {
    return numberOfDecodingIntervals_;
  }

  unsigned int getInitNumberOfKnots() const {
    return initNumberOfKnots_;
  }
  
  unsigned int getMaxNumberOfKnots() const {
    return maxNumberOfKnots_;
  }
  
  unsigned int getMissingBlocksLength() const {
    return skipMissingBlocks_;
  }

  unsigned int getFragmentSize() const {
    return fragmentSize_;
  }
  
  size_t getNumberOfThreads() const {
    return numberOfThreads_;
  }
  
  const std::vector< unsigned int >& getBreakpointStarts() const {
    return breakpointStarts_;
  }
  
  const std::vector< unsigned int >& getBreakpointEnds() const {
    return breakpointEnds_;
  }
  
  const std::vector< unsigned int >& getDiploidIndices() const {
    return diploidIndices_;
  }
  
  bool printSeqs() {
    return printSeqs_;
  }
  
  bool resumeOptim() {
    return resumeOptim_;
  }
  
  bool optimize() {
    return optimize_;
  }
  
  bool enforceFlatDemo() {
    return enforceFlatDemography_;
  }

  bool computeCI() {
    return computeCI_;
  }
  
  bool computeCovar() {
    return computeCovar_;
  }

  bool decode() {
    return decode_;
  }
  
  bool decodeDiploidsParallel() {
    return decodeDiploidsParallel_;
  }
  
  bool decodeBreakpointsParallel() {
    return decodeBreakpointsParallel_;
  }
  
  bool relativeStopCond() {
    return relativeTolerance_;
  }
  bool restrictedModulation() {
    return restrictedModulation_;
  }

  double getFunctionTolerance() {
    return functionTolerance_;
  }
  
  double getParametersTolerance() {
    return parametersTolerance_;
  }
  
  double getTmax() {
    return tMax_;
  }
  
  double getVcfFilterQuality() {
    return filterQualityVcf_;
  }

};

#endif
