#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#include "wrap_monitor.hh"
#include "wrap_beamdist.hh"
#include "wrap_bunch.hh"

#include "BeamDist.hh"

#include <iostream>

//using namespace OrbitUtils;

namespace wrap_monitor{

#ifdef __cplusplus
extern "C" {
#endif

    //---------------------------------//
    // Python BeamDist class definition //
    // ------------------------------- //

    //------------------------------------------------------------------------
    // Member functions for python class wrapping an instance of BeamDist

    // Constructor for python class. It never will be called directly.
    static PyObject* BeamDist_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        // std::cout << "Creating new python beamdist object." << std::endl;
        return (PyObject *) self;
    }

    // Initializer for python ecooler class - Creates the BeamDist instance.
    // This is implementation of the __init__ method beamdist(npart, trig_mode)
    static int BeamDist_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
    {
        int npart; TRIG_TYPE trig_mode; // Declare constructor arguments
        // Declare keywords
        static char *kwlist[] = {(char*)"npart", (char*)"trig_mode", NULL};
        // Try parsing the input parameters
        if(PyArg_ParseTupleAndKeywords(args, kwds, "ii:__init__", kwlist, &npart, &trig_mode))
        {
            self->cpp_obj = new BeamDist(npart, trig_mode);
        }
        else return 1; // Parse failed!

        ((BeamDist*) self->cpp_obj)->setPyWrapper((PyObject*) self); // Point to Python object
        // std::cout << "Created new BeamDist instance." << std::endl;
        return 0;
    }

    // Destructor for python beamdist class - Deletes the BeamDist instance.
    static void BeamDist_del(pyORBIT_Object* self)
    {
        // Get the internal C++ object
        BeamDist* cpp_BeamDist = (BeamDist*) self->cpp_obj;
        // Delete the BeamDist instance if present
        if(cpp_BeamDist != NULL)
            delete cpp_BeamDist;
        // Finally destroy the python object
        self->ob_type->tp_free((PyObject*)self);
        //std::cout << "Destroyed BeamDist instance." << std::endl << std::flush;
    }

    // Interface to BeamDist::trigger()
    static PyObject* BeamDist_trigger(PyObject *self, PyObject *args)
    {
        pyORBIT_Object* pyBeamDist = (pyORBIT_Object*) self; // Get the python beamdist object
        BeamDist* cpp_BeamDist = (BeamDist*) pyBeamDist->cpp_obj; // Get the internal BeamDist instance
        cpp_BeamDist->trigger(); // Trigger acquisition!
        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Interface to BeamDist::trackBunch(Bunch* bunch)
    static PyObject* BeamDist_trackBunch(PyObject *self, PyObject *args)
    {
        pyORBIT_Object* pyBeamDist = (pyORBIT_Object*) self; // Get the python beamdist object
        BeamDist* cpp_BeamDist = (BeamDist*) pyBeamDist->cpp_obj; // Get the internal BeamDist instance
        PyObject* pyBunch; // Declare a pointer to the python bunch object

        // Validate arguments
        if(!PyArg_ParseTuple(args,"O:trackBunch",&pyBunch)) // Try to obtain the argument
            ORBIT_MPI_Finalize("monitor.beamdist: trackBunch(bunch) takes 1 argument.");
        // Obtain the python bunch type
        PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
        // Check whether the python object is of the bunch type
        if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type))
            ORBIT_MPI_Finalize("monitor.beamdist: trackBunch(bunch) - The argument is not a bunch.");
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
        cpp_BeamDist->trackBunch(cpp_bunch); // Finally call the tracker!

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Copy the distribution data stored in BeamDist
    static PyObject* BeamDist_copydata(PyObject *self, PyObject *args)
    {
        int i, nrows, nparticles; // Number of rows in the data buffer
        pyORBIT_Object* pyBeamDist = (pyORBIT_Object*) self; // Get the python beamdist object
        BeamDist* cpp_BeamDist = (BeamDist*) pyBeamDist->cpp_obj; // Get the internal BeamDist instance
        PyObject *pyDataBuff, *pyDataSingleParticle; // Declare a pointer to the numpy array memory object

        // Validate arguments
        if(PyArg_ParseTuple(args,"O!:copydata", &PyList_Type, &pyDataBuff)) // Try to obtain the argument
        {
            nparticles = cpp_BeamDist->getx().size();
            nrows = PyList_Size(pyDataBuff); // Obtain the number of rows, which is the number of particles
            if(nrows < nparticles)
                ORBIT_MPI_Finalize("monitor.beamdist: copydata(databuff). databuff contains less rows than number of particles.");
        }
        else ORBIT_MPI_Finalize("monitor.beamdist: copydata(databuff) takes 1 argument which is a python list.");

        // Extract local references to the data vectors
        std::vector<double> &x = cpp_BeamDist->getx();
        std::vector<double> &y = cpp_BeamDist->gety();
        std::vector<double> &z = cpp_BeamDist->getz();
        std::vector<double> &xp = cpp_BeamDist->getxp();
        std::vector<double> &yp = cpp_BeamDist->getyp();
        std::vector<double> &dE = cpp_BeamDist->getdE();

        // Copy all the data to the list object
        for(i=0; i<nparticles; i++)
        {
            pyDataSingleParticle = PyList_GetItem(pyDataBuff, i); // Get the coords for a single particles
            // Assign the coord data directly!
            PyFloat_AS_DOUBLE(PyList_GetItem(pyDataSingleParticle, 0)) = x[i];
            PyFloat_AS_DOUBLE(PyList_GetItem(pyDataSingleParticle, 1)) = y[i];
            PyFloat_AS_DOUBLE(PyList_GetItem(pyDataSingleParticle, 2)) = z[i];
            PyFloat_AS_DOUBLE(PyList_GetItem(pyDataSingleParticle, 3)) = xp[i];
            PyFloat_AS_DOUBLE(PyList_GetItem(pyDataSingleParticle, 4)) = yp[i];
            PyFloat_AS_DOUBLE(PyList_GetItem(pyDataSingleParticle, 5)) = dE[i];
        }

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    } 

    //------------------------------------------------------------------------
    // monitor.beamdist class manifest - lists methods, members and other
    // properties of the python object
    // Declaration of methods
    static PyMethodDef BeamDistClassMethods[] = {
		{"trigger", BeamDist_trigger, METH_VARARGS, "Trigger data acquisition. - trigger()"},
        { "trackBunch", BeamDist_trackBunch, METH_VARARGS, "Track the bunch. - trackBunch(bunch)"},
		{"copydata", BeamDist_copydata, METH_VARARGS, "Copy the acquired data. - copydata()"},
        {NULL}
        //{NULL, NULL}
    };

    // Declaration of members
    static PyMemberDef BeamDistClassMembers [] = {
        {NULL}
    };

    // Definition of PyTypeObject object
    static PyTypeObject pyORBIT_BeamDist_Type = {
        PyObject_HEAD_INIT(NULL)
        0, /*ob_size*/
        "beamdist", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) BeamDist_del , /*tp_dealloc*/
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
        "The BeamDist python wrapper", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        BeamDistClassMethods, /* tp_methods */
        BeamDistClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) BeamDist_init, /* tp_init */
        0, /* tp_alloc */
        BeamDist_new, /* tp_new */
    };

    // Initialization function of the BeamDist class
    // This makes it known to python
    void initBeamDist(PyObject* pymodule)
    {
        if (PyType_Ready(&pyORBIT_BeamDist_Type) < 0) return;
        Py_INCREF(&pyORBIT_BeamDist_Type);
        PyModule_AddObject(pymodule, "beamdist", (PyObject *)&pyORBIT_BeamDist_Type);
    }


#ifdef __cplusplus
}
#endif

} //end of namespace wrap_monitor
