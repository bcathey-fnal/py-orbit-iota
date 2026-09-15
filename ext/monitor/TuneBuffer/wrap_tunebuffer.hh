#ifndef WRAP_TUNEBUFFER_H
#define WRAP_TUNEBUFFER_H

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_monitor{
    void initTuneBuffer(PyObject* pymodule);
  }

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // WRAP_TUNEBUFFER_H
