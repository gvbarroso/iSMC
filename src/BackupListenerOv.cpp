/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 02/07/2019
 * Last modified: 02/07/2019
 *
 */


#include "BackupListenerOv.h"

using namespace std;
using namespace bpp;

//prints list of current best parameter values (in original space) to backup file
void BackupListenerOv::optimizationStepPerformed(const OptimizationEvent& event) {
    
  const auto reparamFun = event.getOptimizer() -> getFunction();
  ParameterList pl;
  try {
    pl = dynamic_cast<const ReparametrizationFunctionWrapper&>(*reparamFun).function().getParameters();
  } catch(bad_cast& e) {
    const auto reparamFun2 = dynamic_cast<const ThreePointsNumericalDerivative&>(*reparamFun).getFunction();
    pl = dynamic_cast<const ReparametrizationFunctionWrapper&>(*reparamFun2).function().getParameters();
  }

  ofstream bck(backupFile_.c_str(), ios::out);
  double AIC = 2. * static_cast< double >(pl.size()) + 2. * event.getOptimizer()->function().getValue();
  bck << "AIC = " << setprecision(20) << AIC << endl << endl;
  
  for(size_t i = 0; i < pl.size(); ++i) {
    bck << pl[i].getName() << " " <<  setprecision(20) << pl[i].getValue() << endl;
  }
  
  bck.close();
}
