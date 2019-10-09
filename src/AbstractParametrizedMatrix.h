/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 01/01/2019
 *
 */


#ifndef _ABSTRACTPARAMETRIZEDMATRIX_H_
#define _ABSTRACTPARAMETRIZEDMATRIX_H_

#include <string>
#include <vector>

#include "ParametrizedMatrixI.h"


class AbstractParametrizedMatrix: 
  public ParametrizedMatrixI {
public:
  AbstractParametrizedMatrix():
  ParametrizedMatrixI()
  { }
  
public:
  bpp::VVdouble& getExpectedMatrix() {
    return expectedMatrix_;
  }   
  const bpp::VVdouble& getExpectedMatrix() const {
    return expectedMatrix_;
  }
  
  //to be used in baum-welch
  void setMatrix(const bpp::VVdouble& matrix) {
    expectedMatrix_ = matrix;
  }
      
  void scaleMatrix(bpp::VVdouble& matrix, double scale);

};

#endif
