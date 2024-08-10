#ifndef WRAP_MONITOR_H
#define WRAP_MONITOR_H

#include "Python.h"
#include <vector>

#include "wrap_beamdist.hh"
#include "wrap_tunedetector.hh"
// #include "wrap_schottky.hh"
// #include "wrap_bpm.hh"

namespace wrap_monitor{

#ifdef __cplusplus
extern "C" {
#endif
    // Register module so that python knows about it.
    void initmonitor(void);
	// Convenience function to obtain the python type from class name
    PyObject* getMonitorType(const char* name);

#ifdef __cplusplus
}
#endif  // __cplusplus

} //end of namespace wrap_monitor

#endif // WRAP_MONITOR_H
