#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#include "wrap_monitor.hh"
#include "wrap_tunedetector.hh"
#include "wrap_bunch.hh"

#include "TuneDetector.hh"

#include <iostream>

//using namespace OrbitUtils;

namespace wrap_monitor{

#ifdef __cplusplus
extern "C" {
#endif

    //--------------------------------------//
    // Python TuneDetector class definition //
    // ------------------------------------ //

    //------------------------------------------------------------------------
    // Member functions for python class wrapping an instance of TuneDetector

    // Constructor for python class. It never will be called directly.
    static PyObject* TuneDetector_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        return (PyObject *) self;
    }

    // Initializer for python tunedetector class - Creates the TuneDetector instance.
    // This is implementation of the __init__ method TuneDetector(int npart, double T[][4])
    static int TuneDetector_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
    {
        int i,j, npart; // Array indices and number of particles
        // Try parsing the input parameters
        if(!PyArg_ParseTuple(args, "i:__init__", &npart))
            ORBIT_MPI_Finalize("monitor.tunedetector(npart): npart is an integer.");
        self->cpp_obj = new TuneDetector(npart); // Create a new instance of the C++ object
        ((TuneDetector*) self->cpp_obj)->setPyWrapper((PyObject*) self); // Point to Python object
        return 0;
    }

    // Destructor for python tunedetector class - Deletes the TuneDetector instance.
    static void TuneDetector_del(pyORBIT_Object* self)
    {
        // Get the internal C++ object
        TuneDetector* cpp_TuneDetector = (TuneDetector*) self->cpp_obj;
        // Delete the TuneDetector instance if present
        if(cpp_TuneDetector != NULL)
            delete cpp_TuneDetector;
        // Finally destroy the python object
        self->ob_type->tp_free((PyObject*)self);
    }

    // Interface to TuneDetector::reset(int npart)
    static PyObject* TuneDetector_reset(PyObject *self, PyObject *args)
    {
        pyORBIT_Object* pyTuneDetector = (pyORBIT_Object*) self; // Get the python tunedetector object
        TuneDetector* cpp_TuneDetector = (TuneDetector*) pyTuneDetector->cpp_obj; // Get the internal TuneDetector instance
        // Validate argument
        int npart;
        if(!PyArg_ParseTuple(args,"i:reset",&npart)) // Try to obtain the argument
            ORBIT_MPI_Finalize("monitor.tunedetector: reset(npart) takes 1 argument.");
        cpp_TuneDetector->reset(npart); // Reset the tune detector
        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Interface to TuneDetector::trackBunch(Bunch* bunch)
    static PyObject* TuneDetector_trackBunch(PyObject *self, PyObject *args)
    {
        int i,j;
        double Tfloquet[4][4]; // Floquet matrix
        pyORBIT_Object* pyTuneDetector = (pyORBIT_Object*) self; // Get the python tunedetector object
        TuneDetector* cpp_TuneDetector = (TuneDetector*) pyTuneDetector->cpp_obj; // Get the internal TuneDetector instance
        // Declare the python objects: floquet matrix, single row of matrix, bunch
        PyObject *pyTfloquet, *pyTfloquetrow, *pyBunch; int isfirst;
        
        // Validate arguments
        if(!PyArg_ParseTuple(args,"OO!i:trackBunch",&pyBunch, &PyList_Type, &pyTfloquet, &isfirst))
            ORBIT_MPI_Finalize("monitor.beamdist: trackBunch(bunch, Tfloquet, isfirst) takes 3 arguments. Tfloquet is a 4x4 python list object.");
        // Obtain the python bunch type
        PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
        // Check whether the pyBunch is of the bunch type
        if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type))
            ORBIT_MPI_Finalize("monitor.beamdist: trackBunch(bunch, Tfloquet, isfirst) - The first argument is not a bunch.");
        
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj; // Get the C++ instance of the bunch object
        if(PyList_Size(pyTfloquet) == 4) // The number of rows should be 4
        {   
            //std::cout << WARNING_PREFIX << "Tfloquet = [";
            for(i=0;i<4;i++)
            {   
                //std::cout << "[";
                pyTfloquetrow = PyList_GetItem(pyTfloquet, i);
                for(j=0;j<4;j++)
                {
                    Tfloquet[i][j] = PyFloat_AS_DOUBLE(PyList_GetItem(pyTfloquetrow, j));
                    //std::cout << Tfloquet[i][j] << ", ";
                }
                //std::cout << "], ";
            }
            //std::cout << "]" << std::endl;
        }
        else
        {
            std::cout << WARNING_PREFIX << "Setting Tfloquet to the identity matrix." << std::endl;
            for(i=0;i<4;i++)
                for(j=0;j<4;j++)
                    if(i==j)Tfloquet[i][j] = 1.0;
                    else Tfloquet[i][j] = 0.0;
        } 
        
        cpp_TuneDetector->trackBunch(cpp_bunch, Tfloquet, isfirst); // Finally call the tracker!

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Copy the tune data stored in TuneDetector
    static PyObject* TuneDetector_gettunes(PyObject *self, PyObject *args)
    {
        // iterator, number of turns, number of rows in the data buffer, number of particles
        int i, nturns, nrows, nparticles;
        pyORBIT_Object* pyTuneDetector = (pyORBIT_Object*) self; // Get the python tunedetector object
        TuneDetector* cpp_TuneDetector = (TuneDetector*) pyTuneDetector->cpp_obj; // Get the internal TuneDetector instance
        PyObject *pyDataBuff, *pyDataSingleParticle; // Declare a pointer to the numpy array memory object

        // Validate arguments
        if(PyArg_ParseTuple(args,"iO!:gettunes", &nturns, &PyList_Type, &pyDataBuff)) // Try to obtain the argument
        {
            nparticles = cpp_TuneDetector->getParticleCount(); // Get the number of particles in the tune detector
            nrows = PyList_Size(pyDataBuff); // Obtain the number of rows, which should be the number of particles
            if(nrows < nparticles) // Complain!
                ORBIT_MPI_Finalize("monitor.tunedetector: gettunes(nturns, databuff). databuff contains less rows than number of particles.");
        }
        else ORBIT_MPI_Finalize("monitor.tunedetector: gettunes(nturns, databuff) takes 2 arguments. databuff is a python list");

        // Tune data vectors
        std::vector<double> mux, muy, muz;
        // Extract the tune data
        cpp_TuneDetector->getTuneX(mux, nturns);
        cpp_TuneDetector->getTuneY(muy, nturns);
        cpp_TuneDetector->getTuneZ(muz, nturns);

        // Copy all the data to the list object
        for(i=0; i<nparticles; i++)
        {
            pyDataSingleParticle = PyList_GetItem(pyDataBuff, i); // Get the coords for a single particles
            // Assign the coord data directly! BEWARE!!! This assumes a very specific list structure.
            // Should probably do this differently!
            PyFloat_AS_DOUBLE(PyList_GetItem(pyDataSingleParticle, 6)) = mux[i];
            PyFloat_AS_DOUBLE(PyList_GetItem(pyDataSingleParticle, 7)) = muy[i];
            PyFloat_AS_DOUBLE(PyList_GetItem(pyDataSingleParticle, 8)) = muz[i];
        }

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    } 

    //------------------------------------------------------------------------
    // monitor.tunedetector class manifest - lists methods, members and other
    // properties of the python object
    // Declaration of methods
    static PyMethodDef TuneDetectorClassMethods[] = {
		{"reset", TuneDetector_reset, METH_VARARGS, "Reset the tune counters. - reset(npart)"},
        { "trackBunch", TuneDetector_trackBunch, METH_VARARGS, "Track the bunch. - trackBunch(bunch, Tfloquet)"},
		{"gettunes", TuneDetector_gettunes, METH_VARARGS, "Obtain the tunes. - gettunes()"},
        {NULL}
        //{NULL, NULL}
    };

    // Declaration of members
    static PyMemberDef TuneDetectorClassMembers [] = {
        {NULL}
    };

    // Definition of PyTypeObject object
    static PyTypeObject pyORBIT_TuneDetector_Type = {
        PyObject_HEAD_INIT(NULL)
        0, /*ob_size*/
        "tunedetector", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) TuneDetector_del , /*tp_dealloc*/
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
        "The TuneDetector python wrapper.", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        TuneDetectorClassMethods, /* tp_methods */
        TuneDetectorClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) TuneDetector_init, /* tp_init */
        0, /* tp_alloc */
        TuneDetector_new, /* tp_new */
    };

    // Initialization function of the TuneDetector class
    // This makes it known to python
    void initTuneDetector(PyObject* pymodule)
    {
        if (PyType_Ready(&pyORBIT_TuneDetector_Type) < 0) return;
        Py_INCREF(&pyORBIT_TuneDetector_Type);
        PyModule_AddObject(pymodule, "tunedetector", (PyObject *)&pyORBIT_TuneDetector_Type);
    }


#ifdef __cplusplus
}
#endif

} //end of namespace wrap_monitor
