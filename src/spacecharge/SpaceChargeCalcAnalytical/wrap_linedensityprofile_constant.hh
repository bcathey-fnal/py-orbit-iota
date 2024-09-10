#ifndef WRAP_CONSTANT_LINEDENSITY_PROFILE_H
#define WRAP_CONSTANT_LINEDENSITY_PROFILE_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_spacecharge
  {
    void initConstantLineDensityProfile(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
