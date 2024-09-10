#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_linedensityprofile_gaussian.hh"
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
// Python GaussianLineDensityProfile class definition
//---------------------------------------------------------


// Constructor for python class wrapping GaussianLineDensityProfile instance
// It never will be called directly

static PyObject* GaussianLineDensityProfile_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  pyORBIT_Object* self;
  self = (pyORBIT_Object *) type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  // std::cerr << "The GaussianLineDensityProfile new has been called!" << std::endl;
  return (PyObject *) self;
}


// Initializator for python  GaussianLineDensityProfile class
// This is implementation of the __init__ method

static int GaussianLineDensityProfile_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
{
  double blength_rms;
  if(!PyArg_ParseTuple(args,"d:__init__", &blength_rms))
  {
	ORBIT_MPI_Finalize("PyGaussianLineDensityProfile - GaussianLineDensityProfile(bunchlength_rms) - constructor needs parameters.");
  }
  self->cpp_obj = new GaussianLineDensityProfile(blength_rms);
  ((GaussianLineDensityProfile*) self->cpp_obj)->setPyWrapper((PyObject*) self);
  // std::cerr << "The GaussianLineDensityProfile __init__ has been called!" << std::endl;
  return 0;
}


// getLocalLineDensityFactor(z)
static PyObject* GaussianLineDensityProfile_getLocalLineDensityFactor(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyGaussianLineDensityProfile = (pyORBIT_Object*) self;
  GaussianLineDensityProfile* cpp_GaussianLineDensityProfile = (GaussianLineDensityProfile*) pyGaussianLineDensityProfile->cpp_obj;
  double z;
  if(!PyArg_ParseTuple(args,"d:getLocalLineDensityFactor",&z))
  {
    ORBIT_MPI_Finalize("PyGaussianLineDensityProfile - getLocalLineDensityFactor(z) - parameters are needed.");
  }
  return Py_BuildValue("d", cpp_GaussianLineDensityProfile->getLocalLineDensityFactor(z));
}


// setBunchLength(bunchlength_rms)
static PyObject* GaussianLineDensityProfile_setBunchLength(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyGaussianLineDensityProfile = (pyORBIT_Object*) self;
  GaussianLineDensityProfile* cpp_GaussianLineDensityProfile = (GaussianLineDensityProfile*) pyGaussianLineDensityProfile->cpp_obj;
  double bunchlength_rms;
  if(!PyArg_ParseTuple(args,"d:setBunchLength",&bunchlength_rms))
  {
    ORBIT_MPI_Finalize("PyGaussianLineDensityProfile - setBunchLength(bunchlength_rms) - parameters are needed.");
  }
  cpp_GaussianLineDensityProfile->setBunchLength(bunchlength_rms);
  Py_INCREF(Py_None);
  return Py_None;    
}


//-----------------------------------------------------
// Destructor for python GaussianLineDensityProfile class (__del__ method).
//-----------------------------------------------------

static void GaussianLineDensityProfile_del(pyORBIT_Object* self)
{
  //std::cerr << "The GaussianLineDensityProfile __del__ has been called!" << std::endl;
  GaussianLineDensityProfile* cpp_GaussianLineDensityProfile = (GaussianLineDensityProfile*) self->cpp_obj;
  delete cpp_GaussianLineDensityProfile;
  self->ob_type->tp_free((PyObject*)self);
}


// Definition of the methods of the python GaussianLineDensityProfile wrapper class
// They will be vailable from python level

static PyMethodDef GaussianLineDensityProfileClassMethods[] =
{
  {"getLocalLineDensityFactor",  GaussianLineDensityProfile_getLocalLineDensityFactor,  METH_VARARGS, "returns the LineDensityFactor at position z"},
  {"setBunchLength",             GaussianLineDensityProfile_setBunchLength,             METH_VARARGS, "sets the bunch length"},
  {NULL}
};


// Definition of the memebers of the python GaussianLineDensityProfile wrapper class
// They will be vailable from python level

static PyMemberDef GaussianLineDensityProfileClassMembers [] =
{
  {NULL}
};

// New python GaussianLineDensityProfile wrapper type definition

static PyTypeObject pyORBIT_GaussianLineDensityProfile_Type =
{
  PyObject_HEAD_INIT(NULL)
  0, /*ob_size*/
  "GaussianLineDensityProfile", /*tp_name*/
  sizeof(pyORBIT_Object), /*tp_basicsize*/
  0, /*tp_itemsize*/
  (destructor) GaussianLineDensityProfile_del , /*tp_dealloc*/
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
  "The GaussianLineDensityProfile python wrapper", /* tp_doc */
  0, /* tp_traverse */
  0, /* tp_clear */
  0, /* tp_richcompare */
  0, /* tp_weaklistoffset */
  0, /* tp_iter */
  0, /* tp_iternext */
  GaussianLineDensityProfileClassMethods, /* tp_methods */
  GaussianLineDensityProfileClassMembers, /* tp_members */
  0, /* tp_getset */
  0, /* tp_base */
  0, /* tp_dict */
  0, /* tp_descr_get */
  0, /* tp_descr_set */
  0, /* tp_dictoffset */
  (initproc) GaussianLineDensityProfile_init, /* tp_init */
  0, /* tp_alloc */
  GaussianLineDensityProfile_new, /* tp_new */
};


//--------------------------------------------------
// Initialization function of the pyGaussianLineDensityProfile class
// It will be called from SpaceCharge wrapper initialization
//--------------------------------------------------

void initGaussianLineDensityProfile(PyObject* module)
{
  if (PyType_Ready(&pyORBIT_GaussianLineDensityProfile_Type) < 0) return;
  Py_INCREF(&pyORBIT_GaussianLineDensityProfile_Type);
  PyModule_AddObject(module, "GaussianLineDensityProfile", (PyObject *)&pyORBIT_GaussianLineDensityProfile_Type);
  // std::cout << "debug GaussianLineDensityProfile added!" << std::endl;
}


#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
