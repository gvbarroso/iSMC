/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 02/07/2019
 * Last modified: 02/07/2019
 *
 */

#ifndef _BACKUPLISTENEROV_H_
#define _BACKUPLISTENEROV_H_

#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/Functions.h>


class BackupListenerOv:
  public bpp::OptimizationListener {

private:
  std::string backupFile_;

public:
  BackupListenerOv(const std::string& backupFile):
    backupFile_(backupFile) {}

  virtual ~BackupListenerOv() {}

public:
  void optimizationInitializationPerformed(const bpp::OptimizationEvent& event) {}
  
  void optimizationStepPerformed(const bpp::OptimizationEvent& event);

  bool listenerModifiesParameters() const { return false; }
};

#endif
