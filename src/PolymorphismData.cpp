/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 27/03/2025
 *
 */


#include <boost/algorithm/string.hpp>

#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Text/TextTools.h>

#include "PolymorphismData.h"

using namespace bpp;
using namespace std;
using namespace boost::iostreams;
  
void PolymorphismData::processInputSequences() { 
    
  //input sequence file
  filtering_istream seqStream;
  
  if(opt_ -> getSeqCompressionType() == "gzip") {
    seqStream.push(gzip_decompressor());
  }
  
  else if(opt_ -> getSeqCompressionType() == "zip") {
    seqStream.push(zlib_decompressor());
  }
  
  else if(opt_ -> getSeqCompressionType() == "bzip2") {
    seqStream.push(bzip2_decompressor());
  }
  
  else if(opt_ -> getSeqCompressionType() != "none") {
    throw Exception("iSMC::Mis-specified sequence compression type.");
  }  
  
  ApplicationTools::displayResult("Sequence file:", opt_ -> getSequenceFileName());
  ApplicationTools::displayTask("Opening sequence file");
  cout.flush();
  ifstream seqFile(opt_ -> getSequenceFileName(), std::ios_base::in | std::ios_base::binary);  
  seqStream.push(seqFile);
  ApplicationTools::displayTaskDone();
  
  if(!seqFile.is_open()) {
    throw bpp::Exception("iSMC::could not open seq. file: " + opt_ -> getSequenceFileName());    
  }
  
  //mask file (required only for VCF input)
  filtering_istream maskStream;
  
  if(opt_ -> getMaskFileName() != "none") { //if mask file name is provided
    if(opt_ -> getMaskCompressionType() == "gzip") {
      maskStream.push(gzip_decompressor());
    }
    
    else if(opt_ -> getMaskCompressionType() == "zip") {
      maskStream.push(zlib_decompressor());
    }
    
    else if(opt_ -> getMaskCompressionType() == "bzip2") {
      maskStream.push(bzip2_decompressor());        
    }
    
    else if(opt_ -> getMaskCompressionType() != "none") {
      throw bpp::Exception("Mis-specified mask compression type!");
    }
  }
 
  ApplicationTools::displayResult("Mask file:", opt_ -> getMaskFileName());
  ApplicationTools::displayTask("Opening mask file"); 
  ifstream maskFile(opt_ -> getMaskFileName(), std::ios_base::in | std::ios_base::binary); 
  maskStream.push(maskFile);

  if(!maskFile.is_open() && opt_ -> getMaskFileName() != "none") {
    throw bpp::Exception("iSMC::could not open mask file: " + opt_ -> getMaskFileName());    
  }
  ApplicationTools::displayTaskDone();

  //converts DNA sequences into vectors of { 0, 1, 2 } 
  if(opt_ -> getFileType() == "FASTA") {  
    callSnpsFromFasta(seqStream);
  }
    
  else if(opt_ -> getFileType() == "VCF" || opt_ -> getFileType() == "gVCF") {
    if(opt_ -> getMaskFileName() != "none") {
      callSnpsFromVcf(seqStream, maskStream);
    }
    else {
      callSnpsFromVcf(seqStream);
    }
  }
  
  else if(opt_ -> getFileType() == "SNP") {
    callSnpsFromSnpFile(seqStream);
  }
  
  else {
    throw bpp::Exception("iSMC::Mis-specified input file type.");
  }
}

void PolymorphismData::printSequencesToFile() {
    
  ApplicationTools::displayTask("Writing sequences to file(s)");
 
  //size_t snpCounter = 0;
  for(size_t i = 0; i < seqBreakpoints_.size(); ++i) { //for every block
    
    auto focalBlock = seqBreakpoints_[i];
      
    size_t left = focalBlock.first;
    size_t right = focalBlock.second;
        
    //for focal block, writes FASTA file with one sequence per diploid 
    //(hence they have the same length, which is convenient for mapper)
    ofstream seqFile;  
    seqFile.open(opt_ -> getLabel() + ".block." + TextTools::toString(i + 1) + ".fasta.gz",
                 std::ios_base::out | std::ios_base::binary);
    
    filtering_ostream seqStream;
  
    seqStream.push(boost::iostreams::gzip_compressor());
    seqStream.push(seqFile);
    
    for(size_t j = 0; j < indvNames_.size(); ++j) { //for every diploid
        
      seqStream << ">" + indvNames_[j] + ".block." + TextTools::toString(i + 1) << endl;
      
      vector< unsigned char >::const_iterator leftEnd = begin(indvSeqs_[j]) +
                                                        static_cast< vector< unsigned char >::difference_type >(left);
                                                        
      vector< unsigned char >::const_iterator rightEnd = begin(indvSeqs_[j]) +
                                                         static_cast< vector< unsigned char >::difference_type >(right);
    
      vector< unsigned char > diploidFragment(leftEnd, rightEnd);
      
      for(auto& site : diploidFragment) {
        //if(site == 1u) {
        //  ++snpCounter;
        //}
        seqStream << static_cast< size_t >(site);
      }
      seqStream << endl;
    }
    
    boost::iostreams::close(seqStream);
  }
  ApplicationTools::displayTaskDone();
  //cout << snpCounter << " SNPs post-masking." << endl;
}

void PolymorphismData::callSnpsFromSnpFile(filtering_istream& seqInput) {
  
  ApplicationTools::displayTask("Parsing SNP file"); 
  
  //to get breakpoints
  vector< vector< string > > chrTable = readTabFile_(opt_ -> getTabFileName());
    
  vector< size_t > startCoords(chrTable.size());
  vector< size_t > endCoords(chrTable.size());

  for(size_t i = 0; i < chrTable.size(); ++i) {
    
    startCoords[i] = TextTools::to<size_t>(chrTable[i][1]);
    endCoords[i] = TextTools::to<size_t>(chrTable[i][2]);
  }  
  
  setBreakpoints_(startCoords, endCoords);

  vector< unsigned int > containerIndices(opt_ -> getDiploidIndices());    
  indvSeqs_.resize(containerIndices.size() / 2);
  
  for(size_t i = 0; i < indvSeqs_.size(); ++i)
    indvNames_.push_back("diploid_" + TextTools::toString(i + 1));
  
  string line;
  getline(seqInput, line); // skips header
  
  while(getline(seqInput, line)) {
    
    line = TextTools::removeWhiteSpaces(line);
        
    unsigned int diploidIndex = 0;
      
    for(size_t i = 0; i <= containerIndices.size() - 2; i += 2) {
          
      unsigned int first_hap = containerIndices[i];
      unsigned int second_hap = containerIndices[i + 1];
        
      if(first_hap == second_hap)
        throw Exception("iSMC::Comparing a haplotype with itself.");
        
      if((line[first_hap] == '2') || (line[second_hap] == '2'))
        indvSeqs_[diploidIndex].push_back(2u); // missing data

      else if(line[first_hap] == line[second_hap])
        indvSeqs_[diploidIndex].push_back(0u);
        
      else
        indvSeqs_[diploidIndex].push_back(1u);
        
      ++diploidIndex;
    }
  }
  ApplicationTools::displayTaskDone();
}
  
void PolymorphismData::callSnpsFromFasta(filtering_istream& seqInput) { 
  
  ApplicationTools::displayTask("Parsing FASTA file"); 
  
  //to get breakpoints
  vector< vector< string > > chrTable = readTabFile_(opt_ -> getTabFileName());
    
  vector< size_t > startCoords(chrTable.size());
  vector< size_t > endCoords(chrTable.size());

  for(size_t i = 0; i < chrTable.size(); ++i) {
      
    startCoords[i] = TextTools::to<size_t>(chrTable[i][1]);
    endCoords[i] = TextTools::to<size_t>(chrTable[i][2]);
  }  
  
  setBreakpoints_(startCoords, endCoords);
    
  //parses Fasta
  Fasta reader;
  auto alignedSequences = reader.readAlignment(seqInput, AlphabetTools::DNA_ALPHABET);
  
  vector< unsigned int > containerIndices(opt_ -> getDiploidIndices());
  
  for(size_t i = 0; i <= containerIndices.size() - 2; i += 2) {
      
    VectorSiteContainer diploidSequence(AlphabetTools::DNA_ALPHABET);
    
    unsigned int first_hap = containerIndices[i];
    unsigned int second_hap = containerIndices[i + 1];
        
    if(first_hap == second_hap)
      throw Exception("iSMC::Comparing one haplotype with itself.");
    
    auto seq1 = unique_ptr<Sequence>(alignedSequences -> sequence(first_hap).clone());
    diploidSequence.addSequence(seq1->getName(), seq1);
    auto seq2 = unique_ptr<Sequence>(alignedSequences -> sequence(second_hap).clone());
    diploidSequence.addSequence(seq2->getName(), seq2);
      
    vector< unsigned char > diploidSnps(diploidSequence.getNumberOfSites());
    
    for(size_t j = 0; j < diploidSequence.getNumberOfSites(); ++j) {
        
      if(SiteTools::hasGap(diploidSequence.site(j)) || SiteTools::hasUnknown(diploidSequence.site(j))) {
          
        diploidSnps[j] = 2u;
      }
      
      else {
          
        bool snp = !SiteTools::isConstant(diploidSequence.site(j));
        
        if(snp) {
          diploidSnps[j] = 1u;
        }
        
        else { 
          diploidSnps[j] = 0u;
        }
      }
    }
    
    indvSeqs_.push_back(diploidSnps);
    
    indvNames_.push_back(alignedSequences -> sequence(containerIndices[i]).getName() + "-" +
                         alignedSequences -> sequence(containerIndices[i + 1]).getName());
  }
  ApplicationTools::displayTaskDone();
}
   
void PolymorphismData::callSnpsFromVcf(filtering_istream& seqInput, filtering_istream& maskFile) { 
  
  vector< vector< string > > chrTable = readTabFile_(opt_ -> getTabFileName());
    
  vector< size_t > startCoords(chrTable.size());
  vector< size_t > endCoords(chrTable.size());

  for(size_t i = 0; i < chrTable.size(); ++i) {
      
    startCoords[i] = TextTools::to<size_t>(chrTable[i][1]);
    endCoords[i] = TextTools::to<size_t>(chrTable[i][2]);
  }
  
  setBreakpoints_(startCoords, endCoords);
  
  Vcf vcfReader(startCoords, endCoords,
                opt_ -> getFileType(),
                opt_ -> getCallableCode());
    
  vcfReader.readSequences(seqInput);
    
  if(opt_ -> getMaskFileType() == "FASTA") {
    Fasta fastaReader;
    shared_ptr<SequenceContainerInterface> callableMask = fastaReader.readSequences(maskFile, AlphabetTools::DEFAULT_ALPHABET);
    vcfReader.maskSequences(callableMask, seqBreakpoints_);
  }
  
  else if(opt_ -> getMaskFileType() == "BED") {
    vcfReader.maskSequences(maskFile, chrTable, seqBreakpoints_);
  }
  
  else if(opt_ -> getMaskFileType() != "none") {
    throw bpp::Exception("Mis-specified mask file type: " + opt_ -> getMaskFileType());
  }
  
  indvNames_ = vcfReader.getNames();
  indvSeqs_ = vcfReader.getSeqs();
}

void PolymorphismData::callSnpsFromVcf(filtering_istream& seqInput) { 
  
  vector< vector< string > > chrTable = readTabFile_(opt_ -> getTabFileName());
  
  vector< size_t > startCoords(chrTable.size());
  vector< size_t > endCoords(chrTable.size());

  for(size_t i = 0; i < chrTable.size(); ++i) {
      
    startCoords[i] = TextTools::to<size_t>(chrTable[i][1]);
    endCoords[i] = TextTools::to<size_t>(chrTable[i][2]);
  }
  
  setBreakpoints_(startCoords, endCoords);

  Vcf vcfReader(startCoords, endCoords,
                opt_ -> getFileType(),
                opt_ -> getCallableCode());
    
  vcfReader.readSequences(seqInput);
        
  indvNames_ = vcfReader.getNames();
  indvSeqs_ = vcfReader.getSeqs();
}

vector< vector< string > > PolymorphismData::readTabFile_(const string& file) { 
  
  ApplicationTools::displayTask("Parsing Tab file");
  
  vector< vector< string > > table(0, vector< string >(3));
  
  ifstream tabFile;
  tabFile.open(file, ios::in);
  
  if(tabFile.is_open()) {
    
    string line;
    vector< string > splitLine(0);
            
    while(getline(tabFile, line)) { //each line is a different block (chromosome, scaffold)
      
      if(!bpp::TextTools::isEmpty(line)) {
          
        boost::split(splitLine, line, [](char c){ return c == '\t'; });
      
        vector< string > chrInfo(3);
        //removing white spaces as a small safety move
        chrInfo[0] = TextTools::removeWhiteSpaces(splitLine[0]); //chr label
        chrInfo[1] = TextTools::removeWhiteSpaces(splitLine[1]); //start coord
        chrInfo[2] = TextTools::removeWhiteSpaces(splitLine[2]); //end coord
        table.push_back(chrInfo);    
      }
    }
    
    tabFile.close();
  }
  
  else {
    throw bpp::Exception("iSMC::could not open Tab file:" + file);
  }
  
  ApplicationTools::displayTaskDone();
  
  return table;
}

void PolymorphismData::setBreakpoints_(const vector< size_t >& chrStarts,
                                       const vector< size_t >& chrEnds) {
  
  size_t numChr = chrStarts.size();
  size_t left = 0;
  size_t right = 0;
  
  //determines block breakpoints in the contigous WGS (indvSeqs_ vector)
  for(size_t i = 0; i < numChr; ++i) {
    
    right = chrEnds[i] - chrStarts[i] + left;
    seqBreakpoints_.push_back(make_pair(left, right)); //starts at left = 0
    
    left = right + 1; //NOTE: this implies that all sites in the alignment must be used!
  }
}
