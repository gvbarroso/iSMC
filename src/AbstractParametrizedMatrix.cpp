/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 21/11/2018
 *
 */

#include <string>
#include <vector>

#include "AbstractParametrizedMatrix.h"

using namespace bpp;
using namespace std;
  
void AbstractParametrizedMatrix::scaleMatrix(VVdouble& matrix, double scale) {
  for(size_t i = 0; i < matrix.size(); ++i) {
    for(size_t j = 0; j < matrix.size(); ++j) {
      matrix[i][j] *= scale;
    }
  }
}