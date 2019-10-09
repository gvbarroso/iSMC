/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 01/01/2019
 *
 */


#ifndef _PARAMETRIZEDMATRIX_H_
#define _PARAMETRIZEDMATRIX_H_

#include <string>
#include <vector>
#include <memory>

#include <Bpp/Numeric/Parametrizable.h>
#include <Bpp/Numeric/VectorTools.h>

class ParametrizedMatrixI {
protected:
  bpp::VVdouble expectedMatrix_;
  
public:
  ParametrizedMatrixI():
  expectedMatrix_(0, bpp::Vdouble(0))
  { }

  virtual ~ParametrizedMatrixI() {}
  
public:
  virtual void setUpExpectedMatrix() = 0;
  virtual void scaleMatrix(bpp::VVdouble& matrix, double scale) = 0;
  virtual void setMatrix(const bpp::VVdouble& bwMatrix) = 0;
  virtual bpp::VVdouble& getExpectedMatrix() = 0;
  
};

#endif
