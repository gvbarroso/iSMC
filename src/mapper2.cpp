/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 09/09/2018
 * Last modified: 03/09/2019
 *
 */

////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <utility>
#include <cstring>
#include <numeric>
#include <algorithm>
#include <set>
#include <deque>

#include <boost/algorithm/string.hpp>
//Using boost to read from compressed file
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/bimap.hpp>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include "OptionsContainer.h"

using namespace std;
using namespace bpp;

//average and joint
void readConsensusLandscapes() {
    
}

//for parsing rate-TMRCA tables when decoding within time intervals (see printDecodedStates() inside MmPsmc.cpp)
void parsePairOfLandscapes(Vdouble& rateLandscape, vector< size_t >& tmrcaLandscape, const string& fileName) { 
    
  boost::iostreams::filtering_istream landStream; 
  landStream.push(boost::iostreams::gzip_decompressor());
  
  ifstream landFile;
  landFile.open(fileName, std::ios_base::in | std::ios_base::binary);
  
  landStream.push(landFile);
    
  if(landFile.is_open()) {
      
    string line;
    vector< string > splitLine;
    
    while(getline(landStream, line)) {
      
      boost::split(splitLine, line, [](char c){ return c == '\t'; });
      
      size_t tmrcaIndex = TextTools::to<size_t>(splitLine[0]); //the time index
      double rateVal = TextTools::to<double>(splitLine[1]); //the posterior average rate within that time index
      
      tmrcaLandscape.push_back(tmrcaIndex);
      rateLandscape.push_back(rateVal);
    }
    
    boost::iostreams::close(landStream);
  }
  
  else {
    throw bpp::Exception("Mapper::could not open double-landscape file: " + fileName);
  }
}

vector< vector< size_t > > readCoordinatesFromFile(const string& file) { 
    
  vector< vector< size_t > > table(4, vector< size_t >(0));
  
  ifstream coordFile;
  coordFile.open(file, ios::in);
  
  if(coordFile.is_open()) {
      
    string line;
    vector< string > splitLine;
    
    size_t tmp = 0;
    
    getline(coordFile, line); //skips header
    
    while(getline(coordFile, line)) {
        
      boost::split(splitLine, line, [](char c){ return c == '\t'; });
      
      vector< size_t > frag(4);
      
      sscanf(splitLine[0].c_str(), "%zu", &tmp); //frag_start
      table[0].push_back(tmp);
      
      sscanf(splitLine[1].c_str(), "%zu", &tmp); //frag_end
      table[1].push_back(tmp);
      
      sscanf(splitLine[2].c_str(), "%zu", &tmp); //block
      table[2].push_back(tmp);
      
      sscanf(splitLine[3].c_str(), "%zu", &tmp); //file_index
      table[3].push_back(tmp);
    }
    
    coordFile.close();
  }
  
  else {
    throw bpp::Exception("Mapper::could not open file with decode coordinates: " + file);
  }
  
  return table;
}

vector< vector< string > > readDiploidLabelsFromFile(const string& file) { 
    
  vector< vector< string > > table(0);
  
  ifstream labFile;
  labFile.open(file, ios::in);
  
  if(labFile.is_open()) {
      
    string line;
    vector< string > splitLine;
    
    getline(labFile, line); //skips header: decoding_label diploid_id
    
    while(getline(labFile, line)) {
        
      boost::split(splitLine, line, [](char c){return c == '\t';});
      
      vector< string > frag(2);
      
      frag[0] = splitLine[0]; //decoding_label
      frag[1] = splitLine[1]; //diploid_id
      
      table.push_back(frag); 
    }
    
    labFile.close();
  }
  
  else {
    throw bpp::Exception("Mapper::could not open file with diploid labels!");
  }
  
  return table;
}

void writeBinnedLandscapesToFile(const string& rate, const string& label, size_t binSize, const vector< vector< string > >& tabFile,
                                 const vector< vector< string > >& labTab, const VVVdouble& blockBinLands) {
  //aesthetics
  double scaledBin = static_cast< double> (binSize) / 1000.;
  
  string suffix = "kb";
  
  if(scaledBin >= 1000.) {
    scaledBin /= 1000.;
    suffix = "Mb";
  }
  
  //blockBinLands 3D structure: block -> diploid -> bin
  size_t numBlocks = blockBinLands.size();
  size_t numLandscapes = blockBinLands.front().size();
    
  string fileName = label + "." + rate + "." + TextTools::toString(scaledBin) + suffix + ".bedgraph";

  ofstream binnedLands; 
  binnedLands.open(fileName);
  
  //header
  binnedLands << "chrom" << "\t" << "chromStart" << "\t" << "chromEnd" << "\t";
  
  for(size_t i = 0; i < numLandscapes; ++i) {
      
    if(i < labTab.size()) { 
      binnedLands << labTab[i][1] << "\t";
    }
    
    //if these consensus maps exist
    else if(i == labTab.size()) { 
      binnedLands << "sample_mean";
    }
    
    else if(i == labTab.size() + 1){
      binnedLands << "\t" << "joint_decode";
    }
    
    else {
      binnedLands << "\t" << "baum_welch"; //not used ATM
    }
    
  }
  binnedLands << endl;
  
  for(size_t i = 0; i < numBlocks; ++i) {
    
    size_t numWindows = blockBinLands[i].front().size(); //all diploids have same seq length
    
    //to plot coordinates in the first and second columns of the bedgraph
    size_t lowerCut = TextTools::to<size_t>(tabFile[i][5]);  
    size_t upperCut = TextTools::to<size_t>(tabFile[i][6]);
    size_t alnEnd = TextTools::to<size_t>(tabFile[i][2]) - 1;
    size_t starts = TextTools::to<size_t>(tabFile[i][1]) - 1; // -1 to index from 0
        
    if(lowerCut > starts) {
      starts =  lowerCut;
    }
    
    size_t ends = starts + binSize; //bin ends at

    if(upperCut < alnEnd) {
      alnEnd = upperCut;
    }
      
    for(size_t j = 0; j < numWindows; ++j) { //
      
      string chr = "chr" + TextTools::toString(i + 1);
      binnedLands << chr << "\t" << starts << "\t" << ends;
    
      for(size_t k = 0; k < numLandscapes; ++k) {
        binnedLands << "\t" << setprecision(16) << blockBinLands[i][k][j];
      }
    
      binnedLands << endl;
    
      //moves to next bin
      starts = ends;
      ends += binSize;
    
      if(ends > alnEnd) {
        ends = alnEnd;
      }
    }
  }
  
  binnedLands.close();
}

void writeBinnedInfoToFile(const string& prefix, size_t binSize, const vector< vector< string > >& tabFile,
                           const vector< vector< string > >& labTab, const VVVdouble& binData) {
  //aesthetics
  double scaledBin = static_cast< double> (binSize) / 1000.;
  
  string suffix = "kb";
  
  if(scaledBin >= 1000.) {
    scaledBin /= 1000.;
    suffix = "Mb";
  }
  
  string fileName = prefix + TextTools::toString(scaledBin) + suffix + ".bedgraph";
                    
  ofstream binnedFile; 
  binnedFile.open(fileName);
  
  binnedFile << "chrom" << "\t" << "chromStart" << "\t" << "chromEnd" << "\t";
  for(size_t i = 0; i < labTab.size(); ++i) {
    if(i < labTab.size()) { 
      binnedFile << labTab[i][1] << "\t";
    }
  }
  binnedFile << endl;
  
  
  for(size_t i = 0; i < binData.size(); ++i) { //blocks
    
   //to plot coordinates in the first and second columns of the bedgraph
    size_t lowerCut = TextTools::to<size_t>(tabFile[i][5]);  
    size_t upperCut = TextTools::to<size_t>(tabFile[i][6]);
    size_t alnEnd = TextTools::to<size_t>(tabFile[i][2]) - 1;
    size_t starts = TextTools::to<size_t>(tabFile[i][1]) - 1; // -1 to index from 0
       
    if(lowerCut > starts) {
      starts =  lowerCut;
    }
    
    size_t ends = starts + binSize; //bin ends at

    if(upperCut < alnEnd) {
      alnEnd = upperCut;
    }
    
    for(size_t j = 0; j < binData[i].front().size(); ++j) { //bins
    
      string chr = "chr" + TextTools::toString(i + 1);  
      binnedFile << chr << "\t" << starts << "\t" << ends;
    
      for(size_t k = 0; k < binData[i].size(); ++k) { //individuals
        binnedFile << "\t" << binData[i][k][j];
      }
    
      binnedFile << endl;
    
      //moves to next bin
      starts = ends;
      ends += binSize;
    
      if(ends > alnEnd) {
        ends = alnEnd;
      }
    }
  }
  
  binnedFile.close();
}

vector< vector< char > > readSeqsFromFile(const string& fileName) {
  
  boost::iostreams::filtering_istream seqStream; //from boost iostreams
  seqStream.push(boost::iostreams::gzip_decompressor());
  
  ifstream seqFile;
  seqFile.open(fileName, std::ios_base::in | std::ios_base::binary);
  
  seqStream.push(seqFile);
  
  if(seqFile.is_open()) {
      
    Fasta reader;
    auto alnSeqs = reader.readAlignment(seqStream, AlphabetTools::DEFAULT_ALPHABET);
  
    vector< vector< char > > seqs(alnSeqs -> getNumberOfSequences(), vector< char >(0)); //indv -> site
  
    for(size_t i = 0; i < seqs.size(); ++i) {
      
      seqs[i].resize(alnSeqs -> sequence(i).size());
    
      for(size_t j = 0; j < seqs[i].size(); ++j) {
        seqs[i][j] = alnSeqs -> sequence(i).getChar(j)[0];
      }
    }
    
    return seqs;
  }
  
  else {
    throw bpp::Exception("Mapper::could not open sequence file: " + fileName);
  }
}

void readTmrcaRateLandscapes(VVVdouble& allRateLandscapes, vector< vector < vector< size_t > > >& allTmrcaLandscapes,
                             const string& label, const string& rate, size_t numDiploids,
                             const vector< vector< size_t > >& decoordTable) {
  
  ApplicationTools::displayTask("Reading single-nucleotide TMRCA-" + rate + " landscapes");
  
  size_t numUniqueBlocks = set< size_t >(begin(decoordTable[2]), end(decoordTable[2])).size(); // ~chr
  allRateLandscapes.resize(numUniqueBlocks, VVdouble(0, Vdouble(0))); //block -> diploid -> site
  allTmrcaLandscapes.resize(numUniqueBlocks, vector < vector< size_t > >(0, vector< size_t >(0))); //block -> diploid -> site
  
  for(size_t i = 0; i < numUniqueBlocks; ++i) {
    
    auto focalNumFrags = count(begin(decoordTable[2]), end(decoordTable[2]), i + 1);
    VVdouble blockRate(numDiploids, Vdouble(0)); //indv -> site
    vector< vector< size_t > > blockTmrca(numDiploids, vector< size_t >(0)); //indv -> site
    
    for(size_t j = 0; j < numDiploids; ++j) {
    
      Vdouble diploidRateLand(0); //sites
      vector< size_t > diploidTmrcaLand(0); //sites
      
      for(auto k = 0; k < focalNumFrags; ++k) {
        
        Vdouble tmpRate(0);
        vector< size_t > tmpTmrca(0);
        
        string file = label + "_diploid_" + TextTools::toString(j + 1) + "_block_" + TextTools::toString(i + 1) +
                      "_file_index_" + TextTools::toString(k + 1) + "_Tmrca." + rate + ".txt.gz";
        
        parsePairOfLandscapes(tmpRate, tmpTmrca, file);
        
        move(begin(tmpRate), end(tmpRate), back_inserter(diploidRateLand));
        move(begin(tmpTmrca), end(tmpTmrca), back_inserter(diploidTmrcaLand));
      }
      
      blockRate[j] = diploidRateLand;
      blockTmrca[j] = diploidTmrcaLand;
    }
    
    allRateLandscapes[i].resize(blockRate.size(), Vdouble(0));
    for(size_t r = 0; r < blockRate.size(); ++r) {
      allRateLandscapes[i][r].resize(blockRate[r].size());
    }
    allRateLandscapes[i] = blockRate;
    
    allTmrcaLandscapes[i].resize(blockTmrca.size(), vector< size_t >(0));
    for(size_t r = 0; r < blockTmrca.size(); ++r) {
      allTmrcaLandscapes[i][r].resize(blockTmrca[r].size());
    }
    allTmrcaLandscapes[i] = blockTmrca;
  }

  ApplicationTools::displayTaskDone();
}

vector< vector< string > > readTabFile(const string& file) { 
    
  vector< vector< string > > table(0, vector< string > (7));
  
  ifstream tabFile;
  tabFile.open(file, ios::in);
  
  if(tabFile.is_open()) {
    
    string line;
    vector< string > splitLine(0);
            
    while(getline(tabFile, line)) { //each line is a different chromosome
      
      boost::split(splitLine, line, [](char c){ return c == '\t'; });
      
      vector< string > chrInfo(7);
      
      chrInfo[0] = splitLine[0]; //chr label
      chrInfo[1] = splitLine[1]; //start coord
      chrInfo[2] = splitLine[2]; //end coord
      chrInfo[3] = splitLine[3]; //left bps
      chrInfo[4] = splitLine[4]; //right bps
      chrInfo[5] = splitLine[5]; //bottom cut-off
      chrInfo[6] = splitLine[6]; //top cut-off

      table.push_back(chrInfo);    
    }
    
    tabFile.close();
  }
  
  else {
    throw bpp::Exception("iSMC::could not open Tab file:" + file);
  }
     
  return table;
}
 

int main(int argc, char *argv[]) { 
    
  
  cout << endl;
  cout << "******************************************************************" << endl;
  cout << "*                 iSMC Mapper 2, version 1.0.0                   *" << endl;
  cout << "*                                                                *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: G. Barroso                     Last Modif. 24/03/2025 *" << endl;
  cout << "*          J. Dutheil                                            *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;


  if(argc == 1) {
      
    cout << "To use iSMC Mapper, please fill in a params file with the following options and simply ";
    cout << "call it from the command line: ismc_mapper params=[params_file].bpp" << endl << endl;
    
    cout << "dataset_label = // use the same label provided during optimisation & decoding with iSMC " << endl;
    cout << "bin_sizes = // comma separated integers, default = (1000,10000,100000) " << endl;
    cout << "bin_rate =  // which rate to bin ('rho' or 'theta')" << endl;
    cout << "tmrca_interval =  // constrain the posterior average rate by time interval boundaries (min, max)?" << endl;
    cout << "fasta_seqs =  // are there sequences to bin? (if print_seqs is true in iSMC)" << endl;
    cout << "tab_file_path =  // same as for iSMC" << endl;
    
    cout << "For more information, please email gvbarroso@evolbio.mpg.de " << endl << endl;
    
    return(0);
  }  
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  //Reads params file
  BppApplication mapper(argc, argv, "Mapper");
  
  mapper.startTimer();
  map< string, string > params = mapper.getParams();
  
  string label = ApplicationTools::getStringParameter("dataset_label", params, ""); 
  string rate = ApplicationTools::getStringParameter("bin_rate", params, "rho", "", true, 4);
  string tabFilePath = ApplicationTools::getAFilePath("tab_file_path", params, true, true, "");
  vector< size_t > tmrcaBoundaries = ApplicationTools::getVectorParameter< size_t >("tmrca_interval", params, ',', "(0,40)", "", true, 4);
  
  //are there FASTA sequences to bin?
  bool fastaSeqs = ApplicationTools::getBooleanParameter("fasta_seqs", params, false, "", true, 4);
  
  vector< size_t > binSizes = ApplicationTools::getVectorParameter< size_t >("bin_sizes", params, ',', "(10000,100000,1000000)", "", true, 4);
  
  vector< vector< string > > tabFile = readTabFile(tabFilePath);
  vector< vector< string > > labelsTable = readDiploidLabelsFromFile(label + "_diploid_decoding_labels.txt");
  vector< vector< size_t > > decoordTable = readCoordinatesFromFile(label + "_decode_coordinates.txt");  
  
  size_t numDiploids = labelsTable.size();
  
  VVVdouble allRateLandscapes(0, VVdouble(0, Vdouble(0))); //block -> diploid -> site  
  vector< vector < vector< size_t > > > allTmrcaLandscapes(0, vector < vector< size_t > >(0, vector< size_t >(0))); //block -> diploid -> site  
  readTmrcaRateLandscapes(allRateLandscapes, allTmrcaLandscapes, label, rate, numDiploids, decoordTable);
  
  size_t numBlocks = allRateLandscapes.size();
  size_t numLandscapes = allRateLandscapes.front().size(); //numDiploids (+ 2, if numDiploids > 1)
  
  //block -> indv -> site
  vector< vector< vector< char > > > seqs(0, vector< vector< char > >(0, vector< char >(0)));
  
  if(fastaSeqs) {
      
    ApplicationTools::displayTask("Reading FASTA sequences");
    seqs.resize(numBlocks, vector< vector< char > >(numDiploids, vector< char >(0)));
    
    for(size_t i = 0; i < numBlocks; ++i) {
      seqs[i] = readSeqsFromFile(label + ".block." + TextTools::toString(i + 1) + ".fasta.gz");
    }
    ApplicationTools::displayTaskDone();
  }
  
  //to synchronise landscapes at the nucleotide level
  //eventually with an external map (eg DECODE)
  ApplicationTools::displayTask("Trimming landscape(s) according to tab file");

  for(size_t i = 0; i < numBlocks; ++i) { //block ~chr
        
    //this assumes each block is a chr or scaffold
    //in the original seq file and has not been further chopped
    size_t lowerCut = TextTools::to<size_t>(tabFile[i][5]);  
    size_t upperCut = TextTools::to<size_t>(tabFile[i][6]);
    
    size_t alnStart = TextTools::to<size_t>(tabFile[i][1]) - 1; // -1 to index from 0
    size_t alnEnd = TextTools::to<size_t>(tabFile[i][2]) - 1; // -1 to index from 0
    size_t lowerStart = 0;
    
    if(lowerCut >= alnStart) {
      lowerStart = lowerCut - alnStart;
    }
    else {
      cout << endl << "WARNING!! Requested lowerCut precedes alnStart! Using first position instead." << endl;
      lowerCut = alnStart; //to print in the bedgraph
    }
          
    size_t lowerEnd = alnEnd - alnStart;
        
    if(upperCut <= alnEnd) {
      lowerEnd = upperCut - alnStart;
    }
    else {
      cout << "WARNING!! Requested upperCut succeeds alnEnd! Using last position instead." << endl;
      upperCut = alnEnd; //to print in the bedgraph
    }
        
    for(size_t j = 0; j < numLandscapes; ++j) {       
      
      Vdouble::const_iterator leftEnd = begin(allRateLandscapes[i][j]) + static_cast< Vdouble::difference_type >(lowerStart);
      Vdouble::const_iterator rightEnd = begin(allRateLandscapes[i][j]) + static_cast< Vdouble::difference_type >(lowerEnd);
          
      Vdouble tmpRate(leftEnd, rightEnd);
      
      allRateLandscapes[i][j].resize(tmpRate.size());
      allRateLandscapes[i][j] = tmpRate;
      
      vector< size_t >::const_iterator leftFin = begin(allTmrcaLandscapes[i][j]) + static_cast< vector< size_t >::difference_type >(lowerStart);
      vector< size_t >::const_iterator rightFin = begin(allTmrcaLandscapes[i][j]) + static_cast< vector< size_t >::difference_type >(lowerEnd);
    
      vector< size_t > tmpTmrca(leftFin, rightFin);
      
      allTmrcaLandscapes[i][j].resize(tmpTmrca.size());
      allTmrcaLandscapes[i][j] = tmpTmrca;
      
      //NOTE sitting out to test double-landscapes
      if(fastaSeqs && j < numDiploids) { //condition on j since there's no consensus sequence
          
        vector< char >::const_iterator leftSeqs = begin(seqs[i][j]) + static_cast< Vdouble::difference_type >(lowerStart);
        vector< char >::const_iterator rightSeqs = begin(seqs[i][j]) + static_cast< Vdouble::difference_type >(lowerEnd);
        
        vector< char > tmpSeqs(leftSeqs, rightSeqs);
       
        seqs[i][j].resize(tmpSeqs.size());
        seqs[i][j] = tmpSeqs;
      }
    }
  }
  ApplicationTools::displayTaskDone();
  
  
  //Finally binning
  for(size_t i = 0; i < binSizes.size(); ++i) { 
      
    size_t binSize = binSizes[i];
    
    //block -> diploid -> bin
    VVVdouble binMaskPerBlock(numBlocks, VVdouble(numLandscapes, Vdouble(0))); 
    
    //block -> diploid -> bin
    VVVdouble binPiPerBlock(numBlocks, VVdouble(numLandscapes, Vdouble(0))); 
      
    //block -> diploid -> bin
    VVVdouble binRatePerBlock(numBlocks, VVdouble(numLandscapes, Vdouble(0))); 
    
    ApplicationTools::displayTask("Binning maps of size " + TextTools::toString(binSize));
    
    for(size_t j = 0; j < numBlocks; ++j) {
    
      size_t numSites = allRateLandscapes[j].front().size();  
      size_t numWindows = (numSites + binSize - 1) / binSize;
    
      //proportion of missing data in each bin for each individual, for block j
      VVdouble binMask(numDiploids, Vdouble(numWindows)); 
    
      //nucleotide diversity in each bin for each individual, for block j
      VVdouble binPi(numDiploids, Vdouble(numWindows)); 
      
      //all recombination maps for bin size i, for block j
      VVdouble binRate(numLandscapes, Vdouble(numWindows)); 
    
      for(unsigned int k = 0; k < numLandscapes; ++k) { //diploids + avg + joint
        
        for(unsigned int l = 0; l < numWindows; ++l) {
       
          size_t startPos = l * binSize;
          size_t endPos = (l + 1) * binSize;
        
          if(endPos >= numSites) { 
            endPos = numSites;
          }
          
          double binLen = 0.;
          double sumRate = 0.;
          
          for(size_t m = startPos; m < endPos; ++m) {
            //here we condition on the TMRCA index:  
            if(allTmrcaLandscapes[j][k][m] >= tmrcaBoundaries[0] && allTmrcaLandscapes[j][k][m] < tmrcaBoundaries[1]) { 
              sumRate += allRateLandscapes[j][k][m];
              ++binLen;
            }
          }
          
          binRate[k][l] = sumRate / binLen;
          
          cout << rate << " = " << sumRate / binLen << " computed within window " << l + 1 << " using " << binLen << " sites." << endl;
          
          if(k < numDiploids) { //if diploid landscape (not avg or joint-decode)
            
            if(fastaSeqs) {
                
              double missing = static_cast< double >(count(begin(seqs[j][k]) + static_cast<Vdouble::difference_type>(startPos), begin(seqs[j][k]) + static_cast<Vdouble::difference_type>(endPos), '2'));          
              binMask[k][l] =  missing / binLen;
                        
              double snp = static_cast< double >(count(begin(seqs[j][k]) + static_cast<Vdouble::difference_type>(startPos), begin(seqs[j][k]) + static_cast<Vdouble::difference_type>(endPos), '1'));                
              binPi[k][l] = snp / (binLen - missing);
            }
            
          }
        }
      }
      
      binRatePerBlock[j] = binRate;
      binMaskPerBlock[j] = binMask;
      binPiPerBlock[j] = binPi;
    }
    
    //rate landscapes
    writeBinnedLandscapesToFile(rate, label, binSize, tabFile, labelsTable, binRatePerBlock);
      
    if(fastaSeqs) {
          
      //missing data
      string prefix = label + ".missing.prop.";
      writeBinnedInfoToFile(prefix, binSize, tabFile, labelsTable, binMaskPerBlock);
    
      //diversity (pi)
      prefix = label + ".diversity.";
      writeBinnedInfoToFile(prefix, binSize, tabFile, labelsTable, binPiPerBlock);
    }
    
    ApplicationTools::displayTaskDone();
  }
  
  mapper.done();
}
