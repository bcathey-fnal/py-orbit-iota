#include "wrap_electronbeam.hh"

#include <iostream>

//using namespace OrbitUtils;

namespace wrap_scattering{

#ifdef __cplusplus
extern "C" {
#endif

    //--------------------------------------//
    // Python ElectronBeam class definition //
    // ------------------------------------ //

    //------------------------------------------------------------------------
    // Member functions for python class wrapping an instance of ElectronBeam
    
    // Constructor for python class. It never will be called directly.
    static PyObject* ElectronBeam_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        // std::cout << "Creating new python electronbeam object." << std::endl;
        return (PyObject *) self;
    }

    // Initializer for python ecooler class - Creates the ElectronBeam instance.
    // This is implementation of the __init__ method
    // ElectronBeam(double beta_e, double Tpara_e, double Tperp_e, double radius_e, double current_e)
    static int ElectronBeam_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
    {
        double beta, Vjitter, Tperp, radius, current, pipe_radius=0.0; // Declare the constructor arguments
        // Declare keywords
        static char *kwlist[] = {(char*)"beta", (char*)"Vjitter", (char*)"Tperp",
            (char*)"radius", (char*)"current", (char*)"pipe_radius", NULL};
        if(!PyArg_ParseTupleAndKeywords(args, kwds, "ddddd|d:__init__", kwlist,
                    &beta, &Vjitter, &Tperp, &radius, &current, &pipe_radius))
            ORBIT_MPI_Finalize("Constructor is of the form scattering.electronbeam("
                    "beta, Vjitter, Tperp, radius, current, pipe_radius = radius)");
        self->cpp_obj = new ElectronBeam(beta, Vjitter, Tperp, radius, current, pipe_radius); // Allocate the ElectronBeam object
        ((ElectronBeam*) self->cpp_obj)->setPyWrapper((PyObject*) self); // Point to Python object
        return 0;
    }

    // Destructor for python ecooler class - Deletes the ElectronBeam instance.
    static void ElectronBeam_del(pyORBIT_Object* self)
    {
        // Get the internal C++ object
        ElectronBeam* cpp_ElectronBeam = (ElectronBeam*) self->cpp_obj;
        // Delete the ElectronBeam instance if present
        if(cpp_ElectronBeam != NULL)
            delete cpp_ElectronBeam;
        // Finally destroy the python object
        self->ob_type->tp_free((PyObject*)self);
        //std::cout << "Destroyed ElectronBeam instance." << std::endl << std::flush;
    }

    // Interface to ElectronBeam::getBeamPropertiesatRadius(double x, double y, double B, double &nebeam, double *vebeam)
    static PyObject* ElectronBeam_getBeamPropertiesatRadius(PyObject *self, PyObject *args)
    {
        int nVars = PyTuple_Size(args); // Get the number of arguments passed to python
        pyORBIT_Object* pyElectronBeam = (pyORBIT_Object*) self; // Get the python ElectronBeam object
        ElectronBeam* cpp_ElectronBeam = (ElectronBeam*) pyElectronBeam->cpp_obj; // Get the internal ElectronBeam instance
        // Radial position, number density, beam velocity, radial field, voltage depression
        double x, y, radius, B, nebeam, vebeam[3], Er, Vdep;

        // Validate arguments
        if(!PyArg_ParseTuple(args,"ddd:getBeamPropertiesatRadius",&x, &y, &B)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("scattering.electronbeam: getBeamPropertiesatRadius(x, y, B) takes 3 arguments.");
        cpp_ElectronBeam->getBeamPropertiesatRadius(x, y, B, nebeam, vebeam); // Get the density and velocity values
        cpp_ElectronBeam->getStaticField(sqrt(x*x+y*y), Er, Vdep); // Get the field and potential
        return Py_BuildValue("dddddd", nebeam, vebeam[0], vebeam[1], vebeam[2], Er, Vdep);
    }

    // Interface to ElectronBeam::getTpara()
    static PyObject* ElectronBeam_getTpara(PyObject *self, PyObject *args)
    {
        int nVars = PyTuple_Size(args); // Get the number of arguments passed to python
        if(nVars > 0)
            ORBIT_MPI_Finalize("scattering.electronbeam: getTpara() takes no argument.");
        pyORBIT_Object* pyElectronBeam = (pyORBIT_Object*) self; // Get the python ElectronBeam object
        ElectronBeam* cpp_ElectronBeam = (ElectronBeam*) pyElectronBeam->cpp_obj; // Get the internal ElectronBeam instance
        return Py_BuildValue("d", cpp_ElectronBeam->getTpara());
    }

    // Interface to ElectronBeam::getTperp()
    static PyObject* ElectronBeam_getTperp(PyObject *self, PyObject *args)
    {
        int nVars = PyTuple_Size(args); // Get the number of arguments passed to python
        if(nVars > 0)
            ORBIT_MPI_Finalize("scattering.electronbeam: getTperp() takes no argument.");
        pyORBIT_Object* pyElectronBeam = (pyORBIT_Object*) self; // Get the python ElectronBeam object
        ElectronBeam* cpp_ElectronBeam = (ElectronBeam*) pyElectronBeam->cpp_obj; // Get the internal ElectronBeam instance
        return Py_BuildValue("d", cpp_ElectronBeam->getTperp());
    }

    // Interface to ElectronBeam::getOmegap(double r)
    static PyObject* ElectronBeam_getOmegap(PyObject *self, PyObject *args)
    {
        double radius;
        // Validate arguments
        if(!PyArg_ParseTuple(args,"d:getOmegap",&radius)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("scattering.electronbeam: getOmegap(radius) takes 1 argument.");
        pyORBIT_Object* pyElectronBeam = (pyORBIT_Object*) self; // Get the python ElectronBeam object
        ElectronBeam* cpp_ElectronBeam = (ElectronBeam*) pyElectronBeam->cpp_obj; // Get the internal ElectronBeam instance
        return Py_BuildValue("d", cpp_ElectronBeam->getOmegap(radius));
    }

    // Interface to ElectronBeam::getMeanSpacing(double r)
    static PyObject* ElectronBeam_getMeanSpacing(PyObject *self, PyObject *args)
    {
        double radius;
        // Validate arguments
        if(!PyArg_ParseTuple(args,"d:getMeanSpacing",&radius)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("scattering.electronbeam: getMeanSpacing(radius) takes 1 argument.");
        pyORBIT_Object* pyElectronBeam = (pyORBIT_Object*) self; // Get the python ElectronBeam object
        ElectronBeam* cpp_ElectronBeam = (ElectronBeam*) pyElectronBeam->cpp_obj; // Get the internal ElectronBeam instance
        return Py_BuildValue("d", cpp_ElectronBeam->getMeanSpacing(radius));
    }

    
    //------------------------------------------------------------------------
    // scattering.ecooler class manifest - lists methods, members and other
    // properties of the python object
    // Declaration of methods
    static PyMethodDef ElectronBeamClassMethods[] = {
        { "getBeamPropertiesatRadius", ElectronBeam_getBeamPropertiesatRadius, METH_VARARGS,
            "Get the beam density (1/m^3) and velocity (m/s) - getBeamPropertiesatRadius(x, y, B)"},
        { "getTpara", ElectronBeam_getTpara, METH_VARARGS,
            "Get the beam longitudinal temperature (K). - getTpara()"},
        { "getTperp", ElectronBeam_getTperp, METH_VARARGS,
            "Get the beam transverse temperature (K). - getTperp()"},
        { "getOmegap", ElectronBeam_getOmegap, METH_VARARGS,
            "Get the local plasma frequency (1/s) . - getOmegap(radius)"},
        { "getMeanSpacing", ElectronBeam_getMeanSpacing, METH_VARARGS,
            "Get the mean distance between electrons (m) . - getMeanSpacing(radius)"},
        {NULL}
        //{NULL, NULL}
    };

    // Declaration of members
    static PyMemberDef ElectronBeamClassMembers [] = {
        {NULL}
    };

    // Definition of PyTypeObject object
    static PyTypeObject pyORBIT_ElectronBeam_Type = {
        PyObject_HEAD_INIT(NULL)
        0, /*ob_size*/
        "electronbeam", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) ElectronBeam_del , /*tp_dealloc*/
        0, /*tp_print*/
        0, /*tp_getattr*/
        0, /*tp_setattr*/
        0, /*tp_compare*/
        0, /*tp_repr*/
        0, /*tp_as_number*/
        0, /*tp_as_sequence*/
        0, /*tp_as_mapping*/
        0, /*tp_hash */
        0, /*tp_call*/
        0, /*tp_str*/
        0, /*tp_getattro*/
        0, /*tp_setattro*/
        0, /*tp_as_buffer*/
        Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
        "The ElectronBeam python wrapper", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        ElectronBeamClassMethods, /* tp_methods */
        ElectronBeamClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) ElectronBeam_init, /* tp_init */
        0, /* tp_alloc */
        ElectronBeam_new, /* tp_new */
    };

    // Initialization function of the ElectronBeam class
    // This makes it known to python
    void initElectronBeam(PyObject* module)
    {
        if (PyType_Ready(&pyORBIT_ElectronBeam_Type) < 0) return;
        Py_INCREF(&pyORBIT_ElectronBeam_Type);
        PyModule_AddObject(module, "electronbeam", (PyObject *)&pyORBIT_ElectronBeam_Type);
    }

#ifdef __cplusplus
}
#endif

} //end of namespace wrap_scattering
