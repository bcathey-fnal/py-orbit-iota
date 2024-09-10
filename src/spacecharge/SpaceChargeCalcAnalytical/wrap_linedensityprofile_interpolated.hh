#ifndef WRAP_INTERPOLATED_LINEDENSITY_PROFILE_H
#define WRAP_INTERPOLATED_LINEDENSITY_PROFILE_H

#include "Python.h"

#ifdef __cplusplus
extern "C"
{
#endif

  namespace wrap_spacecharge
  {
    void initInterpolatedLineDensityProfile(PyObject* module);
  }

#ifdef __cplusplus
}
#endif

#endif
