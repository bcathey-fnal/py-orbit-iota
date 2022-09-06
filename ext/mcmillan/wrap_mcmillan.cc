// Wrapper code for the C++ McMillan class making it usable from python
// PyORBIT includes
#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#include "wrap_bunch.hh"
// McMillan class specific includes
#include "wrap_mcmillan.hh"
#include "McMillan.hh"
// Standard includes
#include <iostream>

//using namespace OrbitUtils;

// McMillan module doesn't have any methods
static PyMethodDef mcmillanMethods[] = { {NULL,NULL} };

#ifdef __cplusplus
extern "C" {
#endif

    // Initialization function of the mcmillan module.
    // Python looks for this function on seeing "import mcmillan"
    void initmcmillan()
    {
        //create new module
        PyObject* module = Py_InitModule("mcmillan", mcmillanMethods);
        // Register the McMillan class
        initclassMcMillan(module);
    }

    //----------------------------------//
    // Python McMillan class definition //
    //----------------------------------//

    //------------------------------------------------------------------------
    // Member functions for python class wrapping an instance of McMillan

    // Constructor for python class. It never will be called directly.
    static PyObject* McMillan_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        // std::cout << "Creating new python McMillan object." << std::endl;
        return (PyObject *) self;
    }

    // Initializer for python McMillan class - Creates the McMillan instance.
    // This is implementation of the __init__ method.
    // Interface to McMillan::McMillan(double K, double A, double N)
    static int McMillan_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
    {
        // Declare the constructor argument
        double K, A; int O, N;
        // Declare keywords
        static char *kwlist[] = {(char*)"K", (char*)"A", (char*)"O", (char*)"N", NULL};
        // Try parsing the input parameters
        if(PyArg_ParseTupleAndKeywords(args, kwds, "ddii:__init__", kwlist,&K,&A,&O,&N))
            self->cpp_obj = new McMillan(K, A, O, N); // Allocate the McMillan object
        else return 1; // Fail
        ((McMillan*) self->cpp_obj)->setPyWrapper((PyObject*) self); // Point to Python object
        // std::cout << "Created new McMillan instance." << std::endl;
        return 0;
    }

    // Destructor for python McMillan class - Deletes the McMillan instance.
    static void McMillan_del(pyORBIT_Object* self)
    {
        // Get the internal C++ object
        McMillan* cpp_McMillan = (McMillan*) self->cpp_obj;
        // Delete the McMillan instance if present
        if(cpp_McMillan != NULL)
            delete cpp_McMillan;
        // Finally destroy the python object
        self->ob_type->tp_free((PyObject*)self);
        //std::cout << "Destroyed McMillan instance." << std::endl << std::flush;
    }

    // Interface to double McMillan::getStrengthOverLength(), double McMillan::getWidth()
    // and double McMillan::getNslices()
    static PyObject* McMillan_getParams(PyObject *self, PyObject *args)
    {
        pyORBIT_Object* pyMcMillan = (pyORBIT_Object*) self; // Get the python McMillan object
        McMillan* cpp_McMillan = (McMillan*) pyMcMillan->cpp_obj; // Get the internal McMillan instance
        return Py_BuildValue("ddi", cpp_McMillan->getStrengthOverLength(),
                cpp_McMillan->getWidth(), cpp_McMillan->getNslices()); // Return the parameters
    }

    // Interface to void McMillan::setNslices(int N)
    static PyObject* McMillan_setNslices(PyObject *self, PyObject *args)
    {
        int Nslices; // Parameter to obtain and set
        pyORBIT_Object* pyMcMillan = (pyORBIT_Object*) self; // Get the python McMillan object
        McMillan* cpp_McMillan = (McMillan*) pyMcMillan->cpp_obj; // Get the internal McMillan instance
        if(!PyArg_ParseTuple(args,"i:setNslices",&Nslices)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("mcmillan.McMillan: setNslices(N) takes 1 argument.");
        cpp_McMillan->setNslices(Nslices);

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Interface to McMillan::trackBunch(Bunch* bunch, double length)
    static PyObject* McMillan_trackBunch(PyObject *self, PyObject *args)
    {
        pyORBIT_Object* pyMcMillan = (pyORBIT_Object*) self; // Get the python McMillan object
        McMillan* cpp_McMillan = (McMillan*) pyMcMillan->cpp_obj; // Get the internal McMillan instance
        PyObject* pyBunch; // Declare a pointer to the python bunch object

        double length, ks; // We will need the length and solenoid field through which the particles are tracked.

        // Validate arguments
        if(!PyArg_ParseTuple(args,"Odd:trackBunch",&pyBunch,&length,&ks)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("mcmillan.McMillan: trackBunch(bunch,length,ks) takes 3 arguments.");
        // Obtain the python bunch type
        PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
        // Check whether the python object is of the bunch type
        if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type))
            ORBIT_MPI_Finalize("mcmillan.McMillan: trackBunch(bunch,length,ks) - The first argument is not a bunch.");
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
        cpp_McMillan->trackBunch(cpp_bunch,length,ks); // Finally call the tracker!

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }


    //------------------------------------------------------------------------
    // mcmillan.McMillan class manifest - lists methods, members and other
    // properties of the python object
    // Declaration of methods
    static PyMethodDef McMillanClassMethods[] = {
        {"getParams", McMillan_getParams, METH_VARARGS, "Get all parameters of the McMillan lens. K, A, N = getParams()"},
        {"setNslices", McMillan_setNslices, METH_VARARGS, "Set the number of slices - setNslices(N)"},
        {"trackBunch", McMillan_trackBunch, METH_VARARGS, "Track the bunch - trackBunch(bunch,length,ks)"},
        {NULL}
        //{NULL, NULL}
    };

    // Declaration of members
    static PyMemberDef McMillanClassMembers [] = {
        {NULL}
    };

    // Definition of PyTypeObject object
    // Contains all the information relevant to python
    static PyTypeObject pyORBIT_McMillan_Type = {
        PyObject_HEAD_INIT(NULL)
        0, /*ob_size*/
        "McMillan", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) McMillan_del , /*tp_dealloc*/
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
        "The McMillan class python wrapper", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        McMillanClassMethods, /* tp_methods */
        McMillanClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) McMillan_init, /* tp_init */
        0, /* tp_alloc */
        McMillan_new, /* tp_new */
    };


    // Initialization function of the McMillan class
    // This makes it known to python
    void initclassMcMillan(PyObject* module)
    {
        if (PyType_Ready(&pyORBIT_McMillan_Type) < 0) return;
        Py_INCREF(&pyORBIT_McMillan_Type);
        PyModule_AddObject(module, "McMillan", (PyObject *)&pyORBIT_McMillan_Type);
    }

#ifdef __cplusplus
}
#endif

