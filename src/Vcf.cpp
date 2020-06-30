/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 07/07/2018
 * Last modified: 05/09/2019
 *
 */


#include <stdlib.h> 
#include <algorithm>

#include "Vcf.h"
#include "OptionsContainer.h"

using namespace std;
using namespace bpp;
using namespace boost::iostreams;


void Vcf::fetchBasicInfo(filtering_istream& seqInput) {
  
  string line;  
  vector< string > splitLine;

  while(getline(seqInput, line)) {
    
    if(line[0] == '#' && line[1] != '#') { //header
      
      boost::split(splitLine, line, [](char c){ return c == '\t'; }); 
      
      //checks VCF format to update number of columns before genotypes      
      if(any_of(splitLine.begin(), splitLine.end(), [](const string& s){ return s == "FORMAT"; })) {
        numInfoCols_ = 9; 
      }
      
      //assumes there are no cols right of genotypes:
      numDiploids_ = splitLine.size() - numInfoCols_; 
      snpCallings_.resize(numDiploids_);
      indvNames_.resize(numDiploids_);
      
      for(size_t i = 0; i < numDiploids_; ++i) {
        indvNames_[i] = splitLine[numInfoCols_ + i];
      }
      
      break;
    }
  }
}

void Vcf::readSequences(filtering_istream& seqInput) {
  
  cout << "Parsing VCF file..."; cout.flush();  

  fetchBasicInfo(seqInput);
  
  string line;  
  vector< string > splitLine(0);
  vector< string > prevLine(0);

  size_t snpCounter = 0;
  //parses the whole vcf as if it contained a single chromosome
  while(getline(seqInput, line)) {
    
    if(line.at(0) != '#') { //gets lines with data  
        
      boost::split(splitLine, line, [](char c) { return c == '\t'; } );      
      string filter = splitLine[6];

      //at the first position of the VCF (no previous line) NOTE clean this
      if(snpCallings_.front().size() == 0) {
          
        if((filter == "PASS") || (filter == ".") || (atof(filter.c_str()) > filterThresh_)) {
        
          for(size_t i = 0; i < numDiploids_; ++i) {
              
            //calls genotype
            if(splitLine[numInfoCols_ + i][0] == '.') { 
              snpCallings_[i].push_back(2u); //missing data
            } 
    
            else if(splitLine[numInfoCols_ + i][0] == splitLine[numInfoCols_ + i][2]) {
              snpCallings_[i].push_back(0u); //homozygote
            }
            
            //NOTE could add && splitLine[numInfoCols_ + i][2] != '.' here:
            else if(splitLine[numInfoCols_ + i][0] != splitLine[numInfoCols_ + i][2]) { 
              snpCallings_[i].push_back(1u); //heterozygote
              ++snpCounter;
            }
          }
        }
        
        else {
         
          for(size_t i = 0; i < numDiploids_; ++i) {
            snpCallings_[i].push_back(2u);
          }
        }
      }
      
      else {
          
        if((filter == "PASS") || (filter == ".") || (atof(filter.c_str()) > filterThresh_)) {
          parseHighQualitySite_(splitLine, prevLine, snpCounter);
        }
        
        else {
          parseLowQualitySite_(splitLine, prevLine);
        }
      }
      prevLine = splitLine;
    }
  }
  
  cout << " done." << endl;
  //cout << snpCounter << " SNPs pre-masking." << endl;
}

//FASTA mask
void Vcf::maskSequences(shared_ptr< VectorSequenceContainer > callableMask,
                        const vector< pair< size_t, size_t > >& chrBreaks) { 
    
  cout << "Masking sequences..."; cout.flush(); 
  
  size_t numChr = startCoordinates_.size();
  for(size_t i = 0; i < numChr; ++i) {
    
    //start of chr i within the WGS (snpCallings_)
    size_t chrStart = chrBreaks[i].first; 
    //start and end of coordinates (mapped to ref. genome) of chr i within the WGS (snpCallings_)
    size_t startCoord = startCoordinates_.at(i);
    size_t endCoord = endCoordinates_.at(i);
    
    for(size_t j = startCoord; j < endCoord; ++j) {
        
      //NOTE assumes mask is sorted containing only chromosomes in the VCF
      //-1 becase snpCallings_ is indexed from 0 and coordinates in VCFs aren't  
      string state = callableMask -> getSequence(i).getChar(j - 1);
      
      //position of site j in chr i within the contiguous WGS (snpCallings_, indexed from 0)
      //adjusts coordinate since mask file has entire chromosome
      size_t seqPos = j - startCoord + chrStart;
      
      //if the character at position i in maskfile IS NOT PRESENT in callableCode_, we mask it
      if(find(callableCode_.begin(), callableCode_.end(), state) == callableCode_.end()) {
        maskSite_(seqPos); 
      }
    }
  }
  cout << " done." << endl;
}

//BED mask
void Vcf::maskSequences(filtering_istream& mask, vector< vector< string > >& chrTable,
                        const vector< pair< size_t, size_t > >& chrBreaks) { 
    
  cout << "Masking sequences..."; cout.flush(); 
  
  string maskLine;
  vector< string > splitMaskLine(0);
  size_t maskLineCounter = 1;
  
  size_t chrCounter = 0;
  vector< string > splitTabLine = chrTable.at(chrCounter);
  //start coordinate of chr in VCF
  size_t startChrVcf = TextTools::to<size_t>(splitTabLine[1]);
  size_t endChrVcf = TextTools::to<size_t>(splitTabLine[2]);
   
  while(getline(mask, maskLine)) { 
    
    boost::split(splitMaskLine, maskLine, [](char c){ return c == '\t'; });
    
    if(maskLine.empty()) {
      throw bpp::Exception("iSMC::parsing empty line from mask file! Line: " + TextTools::toString(maskLineCounter));
    }
    
    if(splitMaskLine[0] == splitTabLine[0]) { //checks if same chromosome
        
      //start of chr within the WGS (snpCallings_)
      size_t chrStart = chrBreaks[chrCounter].first; 
    
      //start and end of coordinates (mapped to ref. genome) of segment in BED file (indexed from ZERO)
      size_t startBedInterval = TextTools::to<size_t>(splitMaskLine[1]);
      size_t endBedInterval = TextTools::to<size_t>(splitMaskLine[2]);
      
      //only masks site if it's actually present in the VCF
      if(startBedInterval >= startChrVcf && endBedInterval <= endChrVcf) { 
          
        for(size_t j = startBedInterval; j < endBedInterval; ++j) {
          //masks sites present in BED file
          size_t seqPos = j - startChrVcf + chrStart;
          maskSite_(seqPos); 
        }
      }
    }
    
    else { //changes chromosome
        
      ++chrCounter;
            
      splitTabLine = chrTable.at(chrCounter);
      startChrVcf = TextTools::to<size_t>(splitTabLine[1]);
    }
    
    ++maskLineCounter;
  }
  
  cout << " done." << endl;
}
   
void Vcf::maskSite_(size_t pos) {

  for(size_t i = 0; i < numDiploids_; ++i) {
    snpCallings_[i][pos] = 2u; 
  }
}

void Vcf::parseLowQualitySite_(const vector< string >& splitLine,
                               const vector< string >& prevLine) { 
  
  string chr = splitLine[0];
  size_t siteCoord = static_cast< size_t >(stol(splitLine[1]));
  
  string prevChr = prevLine[0];
  size_t prevSiteCoord = static_cast< size_t >(stol(prevLine[1])); 
  
  for(size_t i = 0; i < numDiploids_; ++i) {
    
    if(chr == prevChr) {
      //how many sites in-between positions (avoid -1 in weird cases of repeated coordinate)
      size_t runOfState = max(siteCoord - prevSiteCoord - 1, static_cast< size_t >(0)); 
      snpCallings_[i].insert(snpCallings_[i].end(), runOfState, stdState_); 
    }
    
    snpCallings_[i].push_back(2u);
  }
}
  
void Vcf::parseHighQualitySite_(const vector< string >& splitLine,
                                const vector< string >& prevLine, size_t& snpCounter) {
  
  string chr = splitLine[0];
  int siteCoord = static_cast< int >(stoi(splitLine[1]));
  
  string prevChr = prevLine[0];
  int prevSiteCoord = static_cast< int >(stoi(prevLine[1])); 
  
  for(size_t i = 0; i < numDiploids_; ++i) {
   
    if(chr == prevChr) {
      //how many sites in-between positions (max to avoid -1 in weird cases of repeated coordinate)
      int runOfState = siteCoord - prevSiteCoord - 1;//max(siteCoord - prevSiteCoord - 1, static_cast< size_t >(0)); 
            
      if(runOfState < 0) {
        throw bpp::Exception("iSMC::VCF position is likely duplicated = " + TextTools::toString(siteCoord));
      }
      
      snpCallings_[i].insert(snpCallings_[i].end(), static_cast<size_t>(runOfState), stdState_); 
    }
    
    //calls genotype
    if(splitLine[numInfoCols_ + i][0] == '.') { 
      snpCallings_[i].push_back(2u); //missing data
    }
    
    else if(splitLine[numInfoCols_ + i][0] == splitLine[numInfoCols_ + i][2]) {
      snpCallings_[i].push_back(0u); //homozygote
    }
     
    else if(splitLine[numInfoCols_ + i][0] != splitLine[numInfoCols_ + i][2]) {
      snpCallings_[i].push_back(1u); //heterozygote
      ++snpCounter;
    }
      
    else {
      throw bpp::Exception("iSMC::Could not parse VCF position = " + TextTools::toString(siteCoord));
    }
  }
}
