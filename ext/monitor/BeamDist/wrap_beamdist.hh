#ifndef WRAP_BEAMDIST_H
#define WRAP_BEAMDIST_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_monitor{
    void initBeamDist(PyObject* pymodule);
  }

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // WRAP_BEAMDIST_H
