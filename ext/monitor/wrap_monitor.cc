#include <iostream>
#include "wrap_monitor.hh"

static PyMethodDef monitorMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

    void initmonitor()
    {
        //create new module
        PyObject* pymodule = Py_InitModule("monitor", monitorMethods);
        // add the other classes init
        wrap_monitor::initBeamDist(pymodule);
        wrap_monitor::initTuneDetector(pymodule);
        // wrap_monitor::initSchottky(pymodule);
        // wrap_monitor::initBPM(pymodule);
    }
	
    PyObject* getMonitorType(const char* name)
    {
        PyObject* mod = PyImport_ImportModule("monitor");
        PyObject* pyType = PyObject_GetAttrString(mod,name);
        Py_DECREF(mod);
        Py_DECREF(pyType);
        return pyType;
    }
		
#ifdef __cplusplus
}
#endif
