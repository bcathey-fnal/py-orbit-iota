#ifndef WRAP_ELECTRONBEAM_H
#define WRAP_ELECTRONBEAM_H

#include "Python.h"

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#include "ElectronBeam.hh"
#include "wrap_scattering.hh"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_scattering{
    void initElectronBeam(PyObject* module);
  }
	
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // WRAP_ELECTRONBEAM_H
