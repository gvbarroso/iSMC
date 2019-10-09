/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 12/04/2018
 * Last modified: 20/12/2018
 *
 */


#include "MarkovModulatedSmc.h"

using namespace std;
using namespace bpp;

void MarkovModulatedSmc::scaleCoalescenceRates(size_t NeCategory) {
  //NOTE: this version does not (yet) model heterogeneous Ne in particular epochs
  for(size_t i = 0; i < timeIntervals_.size(); ++i) {
    double scaledLambda = lambdaVector_.getParameterValue("l" +  TextTools::toString(i)) / paramScalings_[2] -> getCategories()[NeCategory]; //ne index in vectors
    lambdaVector_.setParameterValue("l" + TextTools::toString(i), scaledLambda);
  }
}