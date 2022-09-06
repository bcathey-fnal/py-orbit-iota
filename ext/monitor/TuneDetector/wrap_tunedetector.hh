#ifndef WRAP_TUNEDETECTOR_H
#define WRAP_TUNEDETECTOR_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_monitor{
    void initTuneDetector(PyObject* pymodule);
  }

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // WRAP_TUNEDETECTOR_H
