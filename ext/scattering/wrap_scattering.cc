#include "wrap_scattering.hh"

static PyMethodDef scatteringMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

    void initscattering()
    {
        //create new module
        PyObject* module = Py_InitModule("scattering", scatteringMethods);
        //add the other classes init
        wrap_scattering::initECooler(module);
        wrap_scattering::initElectronBeam(module);
    }
	
    PyObject* getScatteringType(const char* name)
    {
        PyObject* mod = PyImport_ImportModule("scattering");
        PyObject* pyType = PyObject_GetAttrString(mod,name);
        Py_DECREF(mod);
        Py_DECREF(pyType);
        return pyType;
    }
		
#ifdef __cplusplus
}
#endif
