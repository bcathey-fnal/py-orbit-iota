#ifndef WRAP_ECOOLER_H
#define WRAP_ECOOLER_H

#include "Python.h"
#include <vector>

#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#include "wrap_scattering.hh"
#include "wrap_bunch.hh"

#include "ECooler.hh"
#include "ElectronBeam.hh"

#ifdef __cplusplus
extern "C" {
#endif

  namespace wrap_scattering{
    void initECooler(PyObject* module);
  }
	
#ifdef __cplusplus
}
#endif // __cplusplus

// Helper functions
int import_python_list(PyObject *pylist, std::vector<double> &vec);
PyObject* export_python_list(std::vector<double> &vec);

#endif // WRAP_ECOOLER_H
