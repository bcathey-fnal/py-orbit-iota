#ifndef WRAP_SCATTERING_H
#define WRAP_SCATTERING_H

#include "Python.h"

#include "wrap_ecooler.hh"
#include "wrap_electronbeam.hh"

namespace wrap_scattering{

#ifdef __cplusplus
extern "C" {
#endif
    // Register module so that python knows about it.
    void initscattering(void);
    // Convenience function to obtain the python type from class name
	PyObject* getScatteringType(const char* name);
	
#ifdef __cplusplus
}
#endif  // __cplusplus

} //end of namespace wrap_scattering

#endif // WRAP_SCATTERING_H
