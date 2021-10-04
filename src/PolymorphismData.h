/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 04/10/2021
 *
 */

#ifndef _POLYMORPHISMDATA_H_
#define _POLYMORPHISMDATA_H_

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

//Using boost to read from compressed file
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/bimap.hpp>

#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/SiteTools.h>

#include "OptionsContainer.h"
#include "Vcf.h"


typedef std::vector< std::pair< size_t, size_t > > Breakpoints;


class PolymorphismData { 
private:
  std::vector< std::vector< unsigned char > > indvSeqs_; //diploid -> site (contiguous WGS)
  std::vector< std::string > indvNames_;
  Breakpoints seqBreakpoints_; //chr breakpoints are shared among all diploids
  std::shared_ptr< OptionsContainer > opt_;

public:
  PolymorphismData(std::shared_ptr< OptionsContainer > smcOptions):
  indvSeqs_(0, std::vector< unsigned char >(0)),
  indvNames_(0),
  seqBreakpoints_(),
  opt_(smcOptions)
  { 
    processInputSequences();
    // makes sure that all independent blocks have at least 1 masked site
    // thus all share the observed states and ziphmm will work
    maskFirstSites_();
  }
  
public: 
  const std::vector< std::vector< unsigned char > >& getSnpCalling() const { 
    return indvSeqs_;
  }
  
  const std::vector< std::string >& getNames() const {
    return indvNames_;
  }
  
  size_t getNumberOfDiploids() const {
    return indvNames_.size();
  }
  
  const Breakpoints& getBreakpoints() const {
    return seqBreakpoints_;
  }

  void processInputSequences();
  
  void printSequencesToFile();
  
  void printMasksToFile();
  
  //from simulations (file w/ mutations mapped into haploids)
  void callSnpsFromSnpFile(boost::iostreams::filtering_istream& seqFile); 
  
  void callSnpsFromFasta(boost::iostreams::filtering_istream& seqInput);
   
  void callSnpsFromVcf(boost::iostreams::filtering_istream& seqInput, boost::iostreams::filtering_istream& maskFile);

  void callSnpsFromVcf(boost::iostreams::filtering_istream& seqInput);

private:  
  void setBreakpoints_(const std::vector< size_t >& left,
                       const std::vector< size_t >& right);

  std::vector< std::vector< std::string > > readTabFile_(const std::string& file);
  
  // masks first position of each chr
  void maskFirstSites_() {
      
    for(size_t i = 0; i < indvSeqs_.size(); ++i)
    {
      for(size_t j = 0; j < seqBreakpoints_.size(); ++j)
      {
        size_t start = seqBreakpoints_[j].first;
        indvSeqs_[i][start] = 2u;
      }
    }
  }
    
};

#endif
