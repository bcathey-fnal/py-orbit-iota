#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_linedensityprofile_constant.hh"
#include "wrap_spacecharge.hh"

#include <iostream>

#include "LineDensityProfile.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge
{

#ifdef __cplusplus
extern "C"
{
#endif

//---------------------------------------------------------
// Python ConstantLineDensityProfile class definition
//---------------------------------------------------------


// Constructor for python class wrapping ConstantLineDensityProfile instance
// It never will be called directly

static PyObject* ConstantLineDensityProfile_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  pyORBIT_Object* self;
  self = (pyORBIT_Object *) type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  // std::cerr << "The ConstantLineDensityProfile new has been called!" << std::endl;
  return (PyObject *) self;
}


// Initializator for python  ConstantLineDensityProfile class
// This is implementation of the __init__ method

static int ConstantLineDensityProfile_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
{
  double blength_rms;
  if(!PyArg_ParseTuple(args,"d:__init__", &blength_rms))
  {
	ORBIT_MPI_Finalize("PyConstantLineDensityProfile - ConstantLineDensityProfile(bunchlength_rms) - constructor needs parameters.");
  }
  self->cpp_obj = new ConstantLineDensityProfile(blength_rms);
  ((ConstantLineDensityProfile*) self->cpp_obj)->setPyWrapper((PyObject*) self);
  // std::cerr << "The ConstantLineDensityProfile __init__ has been called!" << std::endl;
  return 0;
}


// getLocalLineDensityFactor(z)
static PyObject* ConstantLineDensityProfile_getLocalLineDensityFactor(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyConstantLineDensityProfile = (pyORBIT_Object*) self;
  ConstantLineDensityProfile* cpp_ConstantLineDensityProfile = (ConstantLineDensityProfile*) pyConstantLineDensityProfile->cpp_obj;
  double z;
  if(!PyArg_ParseTuple(args,"d:getLocalLineDensityFactor",&z))
  {
    ORBIT_MPI_Finalize("PyConstantLineDensityProfile - getLocalLineDensityFactor(z) - parameters are needed.");
  }
  return Py_BuildValue("d", cpp_ConstantLineDensityProfile->getLocalLineDensityFactor(z));
}


// setBunchLength(bunchlength_rms)
static PyObject* ConstantLineDensityProfile_setBunchLength(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyConstantLineDensityProfile = (pyORBIT_Object*) self;
  ConstantLineDensityProfile* cpp_ConstantLineDensityProfile = (ConstantLineDensityProfile*) pyConstantLineDensityProfile->cpp_obj;
  double blength_full;
  if(!PyArg_ParseTuple(args,"d:setBunchLength",&blength_full))
  {
    ORBIT_MPI_Finalize("PyConstantLineDensityProfile - setBunchLength(blength_full) - parameters are needed.");
  }
  cpp_ConstantLineDensityProfile->setBunchLength(blength_full);
  Py_INCREF(Py_None);
  return Py_None;    
}


//-----------------------------------------------------
// Destructor for python ConstantLineDensityProfile class (__del__ method).
//-----------------------------------------------------

static void ConstantLineDensityProfile_del(pyORBIT_Object* self)
{
  //std::cerr << "The ConstantLineDensityProfile __del__ has been called!" << std::endl;
  ConstantLineDensityProfile* cpp_ConstantLineDensityProfile = (ConstantLineDensityProfile*) self->cpp_obj;
  delete cpp_ConstantLineDensityProfile;
  self->ob_type->tp_free((PyObject*)self);
}


// Definition of the methods of the python ConstantLineDensityProfile wrapper class
// They will be vailable from python level

static PyMethodDef ConstantLineDensityProfileClassMethods[] =
{
  {"getLocalLineDensityFactor",  ConstantLineDensityProfile_getLocalLineDensityFactor,  METH_VARARGS, "returns the LineDensityFactor at position z"},
  {"setBunchLength",             ConstantLineDensityProfile_setBunchLength,             METH_VARARGS, "sets the bunch length"},
  {NULL}
};


// Definition of the memebers of the python ConstantLineDensityProfile wrapper class
// They will be vailable from python level

static PyMemberDef ConstantLineDensityProfileClassMembers [] =
{
  {NULL}
};

// New python ConstantLineDensityProfile wrapper type definition

static PyTypeObject pyORBIT_ConstantLineDensityProfile_Type =
{
  PyObject_HEAD_INIT(NULL)
  0, /*ob_size*/
  "ConstantLineDensityProfile", /*tp_name*/
  sizeof(pyORBIT_Object), /*tp_basicsize*/
  0, /*tp_itemsize*/
  (destructor) ConstantLineDensityProfile_del , /*tp_dealloc*/
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
  "The ConstantLineDensityProfile python wrapper", /* tp_doc */
  0, /* tp_traverse */
  0, /* tp_clear */
  0, /* tp_richcompare */
  0, /* tp_weaklistoffset */
  0, /* tp_iter */
  0, /* tp_iternext */
  ConstantLineDensityProfileClassMethods, /* tp_methods */
  ConstantLineDensityProfileClassMembers, /* tp_members */
  0, /* tp_getset */
  0, /* tp_base */
  0, /* tp_dict */
  0, /* tp_descr_get */
  0, /* tp_descr_set */
  0, /* tp_dictoffset */
  (initproc) ConstantLineDensityProfile_init, /* tp_init */
  0, /* tp_alloc */
  ConstantLineDensityProfile_new, /* tp_new */
};


//--------------------------------------------------
// Initialization function of the pyConstantLineDensityProfile class
// It will be called from SpaceCharge wrapper initialization
//--------------------------------------------------

void initConstantLineDensityProfile(PyObject* module)
{
  if (PyType_Ready(&pyORBIT_ConstantLineDensityProfile_Type) < 0) return;
  Py_INCREF(&pyORBIT_ConstantLineDensityProfile_Type);
  PyModule_AddObject(module, "ConstantLineDensityProfile", (PyObject *)&pyORBIT_ConstantLineDensityProfile_Type);
  // std::cout << "debug ConstantLineDensityProfile added!" << std::endl;
}


#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
