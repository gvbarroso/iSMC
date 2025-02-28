/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 07/07/2018
 * Last modified: 13/05/2019
 *
 */


#ifndef _VCF_
#define _VCF_

#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/bimap.hpp>

//for parsing mask files in FASTA format
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/SiteTools.h>

//Using boost to read from compressed file
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>


class Vcf {
private:
  bool isComplete_; //complete VCF contains all callable positions (all others must be masked)
  unsigned char stdState_; //the state that fills the observed sequence for positions absent from VCF (0 for std VCF, 2 for complete VCF)
  size_t numInfoCols_; //number of columns before genotypes
  size_t numDiploids_; //number of diploid individuals (each with n chromosomes)
  double filterThresh_;
  std::vector< size_t > startCoordinates_; //start positions for each chromosome
  std::vector< size_t > endCoordinates_; //end positions for each chromosome
  std::vector< std::string > callableCode_; //for masking, characters that represent that a position passed filters
  std::vector< std::string > indvNames_; //indv names
  std::vector< std::vector< unsigned char > > snpCallings_; //indv -> site (contiguous WGS)
  
public:
  Vcf(std::vector< size_t > starts,
      std::vector< size_t > ends,
      const std::string& vcfType,
      const std::vector< std::string >& callableCode):
  isComplete_(false),
  stdState_(0u),
  numInfoCols_(8),
  numDiploids_(0),
  filterThresh_(0.),
  startCoordinates_(starts),
  endCoordinates_(ends),
  callableCode_(callableCode),
  indvNames_(0),
  snpCallings_(0, std::vector< unsigned char >(0))
  {
    if(vcfType == "cVCF") { //"complete" VCF
      isComplete_ = true; 
      stdState_ = 2u; //if complete VCF, absent positions are considered missing
    }
  }

public:
  bool isVcfComplete() {
    return isComplete_;
  }

  unsigned char getStandardState() {
    return stdState_;
  }

  size_t getNumberOfInfoColumns() {
    return numInfoCols_;
  }
  
  size_t getNumberOfDiploids() {
    return numDiploids_;
  }
  
  const std::vector< std::string >& getCallableCode() const {
    return callableCode_;
  }
  
  const std::vector< std::string >& getNames() const {
    return indvNames_;
  }

  const std::vector< std::vector< unsigned char > >& getSeqs() const {
    return snpCallings_;
  }

  //gets basic information from the VCF: startCoodtinate_, numDiploids_ and indvNames_	
  //this method exists on its own to make readSequences straightforward (clean) and maybe a tiny bit faster
  void fetchBasicInfo(boost::iostreams::filtering_istream& seqInput);

  void readSequences(boost::iostreams::filtering_istream& seqInput);
  
  //mask in FASTA format
  void maskSequences(std::shared_ptr< bpp::SequenceContainerInterface > callableMask,
                     const std::vector< std::pair< size_t, size_t > >& chrBreaks);
  
  //mask in BED format
  void maskSequences(boost::iostreams::filtering_istream& mask, std::vector< std::vector< std::string > >& chrTable,
                     const std::vector< std::pair< size_t, size_t > >& chrBreaks);

private:      
  void parseLowQualitySite_(const std::vector< std::string >& splitLine,
                            const std::vector< std::string >& prevLine);
  
  void maskSite_(size_t pos);
  
  void parseHighQualitySite_(const std::vector< std::string >& splitLine,
                             const std::vector< std::string >& prevLine,
			     size_t& snpCounter);  

 };

#endif
