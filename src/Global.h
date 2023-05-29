/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 16/04/2018
 * Last modified: 16/04/2018
 *
 */


#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <cstddef>

#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/ParameterList.h>

extern size_t NUMBER_OF_AVAILABLE_THREADS;

bpp::RowMatrix< double >* computeDiagonalHessian(bpp::SecondOrderDerivable& function, const bpp::ParameterList& parameters);

void invertDiagonalMatrix(const bpp::RowMatrix< double >& A, bpp::RowMatrix< double >& O);

//NOTE add help function
   
#endif
