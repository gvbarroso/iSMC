/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 17/04/2018
 * Last modified: 18/04/2018
 *
 */


#include <thread>
#include <iostream>

#include "Global.h"
#include "OptionsContainer.h"

size_t NUMBER_OF_AVAILABLE_THREADS = std::thread::hardware_concurrency(); //updated in ismc.cpp

using namespace std;
using namespace bpp;


RowMatrix< double >* computeDiagonalHessian(DerivableSecondOrder& function, const ParameterList& parameters) {

  size_t n = parameters.size();

  //parameters.printParameters(cout);
  vector< string > variables = parameters.getParameterNames();
  RowMatrix< double >* hessian = new RowMatrix< double >(n, n);

  for(unsigned int i = 0; i < n; ++i) {

    for(unsigned int j = 0; j < n; ++j) {

      if(j == i) {
        (*hessian)(i,j) = function.d2f(variables[i], parameters);
      }

      else {
        (*hessian)(i,j) = 0.;
      }
    }
  }

  return hessian;
}

void invertDiagonalMatrix(const RowMatrix< double >& A, RowMatrix< double >& O) {
  O = A;
  for(unsigned int i = 0; i < A.getNumberOfRows(); ++i) {
    O(i,i) = 1. / A(i,i);
  }
}