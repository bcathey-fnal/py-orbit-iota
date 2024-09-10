#ifndef WRAP_SC_SPACECHARGE_ANALYTIC_GAUSSIAN_H
#define WRAP_SC_SPACECHARGE_ANALYTIC_GAUSSIAN_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_spacecharge{
    void initSpaceChargeCalcAnalyticGaussian(PyObject* module);
  }
	
#ifdef __cplusplus
}
#endif // __cplusplus

#endif 
