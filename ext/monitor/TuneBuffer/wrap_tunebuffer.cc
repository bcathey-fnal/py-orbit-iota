#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#include "wrap_monitor.hh"
#include "wrap_tunebuffer.hh"
#include "wrap_bunch.hh"

#include "TuneBuffer.hh"

#include <iostream>

namespace wrap_monitor{

#ifdef __cplusplus
extern "C" {
#endif

    //------------------------------------//
    // Python TuneBuffer class definition //
    // ---------------------------------- //

    // Constructor for python class. It never will be called directly.
    static PyObject* TuneBuffer_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        return (PyObject *) self;
    }

    // Initializer: tunebuffer(window)
    static int TuneBuffer_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
    {
        int window;
        if(!PyArg_ParseTuple(args, "i:__init__", &window))
            ORBIT_MPI_Finalize("monitor.tunebuffer(window): window is an integer.");
        self->cpp_obj = new TuneBuffer(window);
        ((TuneBuffer*) self->cpp_obj)->setPyWrapper((PyObject*) self);
        return 0;
    }

    // Destructor
    static void TuneBuffer_del(pyORBIT_Object* self)
    {
        TuneBuffer* cpp = (TuneBuffer*) self->cpp_obj;
        if(cpp != NULL) delete cpp;
        self->ob_type->tp_free((PyObject*)self);
    }

    // The C++ bunch behind a python bunch argument
    static Bunch* bunch_of(PyObject* pyBunch, const char* where)
    {
        PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
        if(!PyObject_IsInstance(pyBunch, pyORBIT_Bunch_Type))
            ORBIT_MPI_Finalize(where);
        return (Bunch*) ((pyORBIT_Object*) pyBunch)->cpp_obj;
    }

    // reset(bunch, window)
    static PyObject* TuneBuffer_reset(PyObject *self, PyObject *args)
    {
        TuneBuffer* cpp = (TuneBuffer*) ((pyORBIT_Object*) self)->cpp_obj;
        PyObject* pyBunch; int window;
        if(!PyArg_ParseTuple(args, "Oi:reset", &pyBunch, &window))
            ORBIT_MPI_Finalize("monitor.tunebuffer: reset(bunch, window) takes a bunch and an integer.");
        cpp->reset(bunch_of(pyBunch, "monitor.tunebuffer: reset(bunch, window) - the first argument is not a bunch."), window);
        Py_INCREF(Py_None);
        return Py_None;
    }

    // trackBunch(bunch, isfirst)
    static PyObject* TuneBuffer_trackBunch(PyObject *self, PyObject *args)
    {
        TuneBuffer* cpp = (TuneBuffer*) ((pyORBIT_Object*) self)->cpp_obj;
        PyObject* pyBunch; int isfirst;
        if(!PyArg_ParseTuple(args, "Oi:trackBunch", &pyBunch, &isfirst))
            ORBIT_MPI_Finalize("monitor.tunebuffer: trackBunch(bunch, isfirst) takes a bunch and an integer.");
        cpp->trackBunch(bunch_of(pyBunch, "monitor.tunebuffer: trackBunch(bunch, isfirst) - the first argument is not a bunch."), isfirst != 0);
        Py_INCREF(Py_None);
        return Py_None;
    }

    // The scalars
    static PyObject* TuneBuffer_nslots(PyObject *self, PyObject *args)
    {
        return Py_BuildValue("i", ((TuneBuffer*) ((pyORBIT_Object*) self)->cpp_obj)->getNSlots());
    }
    static PyObject* TuneBuffer_window(PyObject *self, PyObject *args)
    {
        return Py_BuildValue("i", ((TuneBuffer*) ((pyORBIT_Object*) self)->cpp_obj)->getWindow());
    }
    static PyObject* TuneBuffer_turn(PyObject *self, PyObject *args)
    {
        return Py_BuildValue("l", ((TuneBuffer*) ((pyORBIT_Object*) self)->cpp_obj)->getTurn());
    }
    static PyObject* TuneBuffer_latest(PyObject *self, PyObject *args)
    {
        return Py_BuildValue("i", ((TuneBuffer*) ((pyORBIT_Object*) self)->cpp_obj)->getLatest());
    }
    static PyObject* TuneBuffer_largestStep(PyObject *self, PyObject *args)
    {
        return Py_BuildValue("d", ((TuneBuffer*) ((pyORBIT_Object*) self)->cpp_obj)->getLargestStep());
    }

    // A read-write buffer object over a vector's memory, which numpy's
    // frombuffer views without copying; None for an empty vector. The memory
    // belongs to the C++ object, which must outlive the view.
    static PyObject* view_of(void* ptr, size_t nbytes)
    {
        if(nbytes == 0)
        {
            Py_INCREF(Py_None);
            return Py_None;
        }
        return PyBuffer_FromReadWriteMemory(ptr, (Py_ssize_t) nbytes);
    }
#define TUNEBUFFER_VIEW(getter) \
    { \
        TuneBuffer* cpp = (TuneBuffer*) ((pyORBIT_Object*) self)->cpp_obj; \
        if(cpp->getter().size() == 0) return view_of(NULL, 0); \
        return view_of((void*) &cpp->getter()[0], \
                       cpp->getter().size()*sizeof(cpp->getter()[0])); \
    }
    static PyObject* TuneBuffer_coords(PyObject *self, PyObject *args)
    TUNEBUFFER_VIEW(getCoords)
    static PyObject* TuneBuffer_windings(PyObject *self, PyObject *args)
    TUNEBUFFER_VIEW(getWindings)
    static PyObject* TuneBuffer_centroid(PyObject *self, PyObject *args)
    TUNEBUFFER_VIEW(getCentroid)
    static PyObject* TuneBuffer_centroidWindings(PyObject *self, PyObject *args)
    TUNEBUFFER_VIEW(getCentroidWindings)
    static PyObject* TuneBuffer_turns(PyObject *self, PyObject *args)
    TUNEBUFFER_VIEW(getTurns)
    static PyObject* TuneBuffer_isbeam(PyObject *self, PyObject *args)
    TUNEBUFFER_VIEW(getIsBeam)
    static PyObject* TuneBuffer_recordedSlots(PyObject *self, PyObject *args)
    TUNEBUFFER_VIEW(getRecordedSlots)
#undef TUNEBUFFER_VIEW

    // Declaration of methods
    static PyMethodDef TuneBufferClassMethods[] = {
        {"reset", TuneBuffer_reset, METH_VARARGS, "Assign slots to the bunch's particles and size the arrays. - reset(bunch, window)"},
        {"trackBunch", TuneBuffer_trackBunch, METH_VARARGS, "Count at this position, and record if at the start of the ring. - trackBunch(bunch, isfirst)"},
        {"nslots", TuneBuffer_nslots, METH_VARARGS, "Number of slots. - nslots()"},
        {"window", TuneBuffer_window, METH_VARARGS, "Turns held. - window()"},
        {"turn", TuneBuffer_turn, METH_VARARGS, "Turns recorded so far. - turn()"},
        {"latest", TuneBuffer_latest, METH_VARARGS, "Column of the latest turn, -1 before the first. - latest()"},
        {"largestStep", TuneBuffer_largestStep, METH_VARARGS, "Largest |phase step| seen (rad). - largestStep()"},
        {"coords", TuneBuffer_coords, METH_VARARGS, "Buffer over the coordinates, float32 (6, nslots, window). - coords()"},
        {"windings", TuneBuffer_windings, METH_VARARGS, "Buffer over the windings, int8 (2, nslots, window). - windings()"},
        {"centroid", TuneBuffer_centroid, METH_VARARGS, "Buffer over the centroid, float64 (6, window). - centroid()"},
        {"centroidWindings", TuneBuffer_centroidWindings, METH_VARARGS, "Buffer over the centroid's windings, int8 (2, window). - centroidWindings()"},
        {"turns", TuneBuffer_turns, METH_VARARGS, "Buffer over the turn of each column, long (window,). - turns()"},
        {"isbeam", TuneBuffer_isbeam, METH_VARARGS, "Buffer over the beam flags, uint8 (nslots,). - isbeam()"},
        {"recordedSlots", TuneBuffer_recordedSlots, METH_VARARGS, "Buffer over the slots in bunch order at the last record, int32. - recordedSlots()"},
        {NULL}
    };

    // Declaration of members
    static PyMemberDef TuneBufferClassMembers [] = {
        {NULL}
    };

    // Definition of PyTypeObject object
    static PyTypeObject pyORBIT_TuneBuffer_Type = {
        PyObject_HEAD_INIT(NULL)
        0, /*ob_size*/
        "tunebuffer", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) TuneBuffer_del , /*tp_dealloc*/
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
        "The TuneBuffer python wrapper.", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        TuneBufferClassMethods, /* tp_methods */
        TuneBufferClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) TuneBuffer_init, /* tp_init */
        0, /* tp_alloc */
        TuneBuffer_new, /* tp_new */
    };

    // Initialization function of the TuneBuffer class
    void initTuneBuffer(PyObject* pymodule)
    {
        if (PyType_Ready(&pyORBIT_TuneBuffer_Type) < 0) return;
        Py_INCREF(&pyORBIT_TuneBuffer_Type);
        PyModule_AddObject(pymodule, "tunebuffer", (PyObject *)&pyORBIT_TuneBuffer_Type);
    }

#ifdef __cplusplus
}
#endif

} //end of namespace wrap_monitor
