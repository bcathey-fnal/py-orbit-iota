#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_linedensityprofile_interpolated.hh"
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
// Python InterpolatedLineDensityProfile class definition
//---------------------------------------------------------


// Constructor for python class wrapping InterpolatedLineDensityProfile instance
// It never will be called directly

static PyObject* InterpolatedLineDensityProfile_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  pyORBIT_Object* self;
  self = (pyORBIT_Object *) type->tp_alloc(type, 0);
  self->cpp_obj = NULL;
  // std::cerr << "The InterpolatedLineDensityProfile new has been called!" << std::endl;
  return (PyObject *) self;
}


// Initializator for python  InterpolatedLineDensityProfile class
// This is implementation of the __init__ method
static int InterpolatedLineDensityProfile_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
{
  double z_min, z_max;
  PyObject* lambda_list;
  if(!PyArg_ParseTuple(args,"ddO:setInterpolatedLineDensityProfile", &z_min, &z_max, &lambda_list)){
	ORBIT_MPI_Finalize("pyInterpolatedLineDensityProfile - InterpolatedLineDensityProfile(z_min, z_max, lambda[]) - constructor needs parameters.");
  }
  if(!PyList_Check(lambda_list)){
	ORBIT_MPI_Finalize("pyInterpolatedLineDensityProfile - InterpolatedLineDensityProfile(z_min, z_max, lambda[]) - lambda is not a list.");
  }
  int array_length = PyList_Size(lambda_list);
  double lambda[array_length];
  PyObject* list_entry;
  for (int i=0; i < array_length; i++) {
	list_entry = PyList_GetItem(lambda_list, i);
	if (!PyFloat_Check(list_entry) & !PyInt_Check(list_entry) & !PyLong_Check(list_entry)) {
		ORBIT_MPI_Finalize("pyInterpolatedLineDensityProfile - InterpolatedLineDensityProfile(z_min, z_max, lambda[]) - lambda contains invalid entries");
	};
	lambda[i] = PyFloat_AsDouble(list_entry);
	//std::cerr << "lambda[" << i<< "] = " << lambda[i] << std::endl;
  };
  self->cpp_obj = new InterpolatedLineDensityProfile(z_min, z_max, lambda, array_length);
  ((InterpolatedLineDensityProfile*) self->cpp_obj)->setPyWrapper((PyObject*) self);
  // std::cerr << "The InterpolatedLineDensityProfile __init__ has been called!" << std::endl;
  return 0;
};


// setLineDensityProfile
static PyObject* InterpolatedLineDensityProfile_setLineDensityProfile(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyInterpolatedLineDensityProfile = (pyORBIT_Object*) self;
  InterpolatedLineDensityProfile* cpp_InterpolatedLineDensityProfile = (InterpolatedLineDensityProfile*) pyInterpolatedLineDensityProfile->cpp_obj;
  double z_min, z_max;
  PyObject* lambda_list;
  if(!PyArg_ParseTuple(args,"ddO:setLineDensityProfile", &z_min, &z_max, &lambda_list)){
	ORBIT_MPI_Finalize("pyInterpolatedLineDensityProfile - setLineDensityProfile(z_min, z_max, lambda[]) -  needs parameters.");
  }
  if(!PyList_Check(lambda_list)){
	ORBIT_MPI_Finalize("pyInterpolatedLineDensityProfile - setLineDensityProfile(z_min, z_max, lambda[]) - lambda is not a list.");
  }
  int array_length = PyList_Size(lambda_list);
  double lambda[array_length];
  PyObject* list_entry;
  for (int i=0; i < array_length; i++) {
	list_entry = PyList_GetItem(lambda_list, i);
	if (!PyFloat_Check(list_entry) & !PyInt_Check(list_entry) & !PyLong_Check(list_entry)) {
		ORBIT_MPI_Finalize("pyInterpolatedLineDensityProfile - setLineDensityProfile(z_min, z_max, lambda[]) - lambda contains invalid entries");
	};
	lambda[i] = PyFloat_AsDouble(list_entry);
	//std::cerr << "lambda[" << i<< "] = " << lambda[i] << std::endl;
  };
  cpp_InterpolatedLineDensityProfile->setLineDensityProfile(z_min, z_max, lambda, array_length);
  Py_INCREF(Py_None);
  return Py_None;   
};


// getLocalLineDensityFactor(z)
static PyObject* InterpolatedLineDensityProfile_getLocalLineDensityFactor(PyObject *self, PyObject *args)
{
  pyORBIT_Object* pyInterpolatedLineDensityProfile = (pyORBIT_Object*) self;
  InterpolatedLineDensityProfile* cpp_InterpolatedLineDensityProfile = (InterpolatedLineDensityProfile*) pyInterpolatedLineDensityProfile->cpp_obj;
  double z;
  if(!PyArg_ParseTuple(args,"d:getLocalLineDensityFactor",&z))
  {
    ORBIT_MPI_Finalize("PyInterpolatedLineDensityProfile - getLocalLineDensityFactor(z) - parameters are needed.");
  }
  return Py_BuildValue("d", cpp_InterpolatedLineDensityProfile->getLocalLineDensityFactor(z));
}



//-----------------------------------------------------
// Destructor for python InterpolatedLineDensityProfile class (__del__ method).
//-----------------------------------------------------

static void InterpolatedLineDensityProfile_del(pyORBIT_Object* self)
{
  //std::cerr << "The InterpolatedLineDensityProfile __del__ has been called!" << std::endl;
  InterpolatedLineDensityProfile* cpp_InterpolatedLineDensityProfile = (InterpolatedLineDensityProfile*) self->cpp_obj;
  delete cpp_InterpolatedLineDensityProfile;
  self->ob_type->tp_free((PyObject*)self);
}


// Definition of the methods of the python InterpolatedLineDensityProfile wrapper class
// They will be vailable from python level

static PyMethodDef InterpolatedLineDensityProfileClassMethods[] =
{
  {"getLocalLineDensityFactor",  InterpolatedLineDensityProfile_getLocalLineDensityFactor, METH_VARARGS, "returns the LineDensityFactor at position z"},
  {"setLineDensityProfile",      InterpolatedLineDensityProfile_setLineDensityProfile,     METH_VARARGS, "sets the LineDensityProfile"},  
  {NULL}
};


// Definition of the memebers of the python InterpolatedLineDensityProfile wrapper class
// They will be vailable from python level

static PyMemberDef InterpolatedLineDensityProfileClassMembers [] =
{
  {NULL}
};

// New python InterpolatedLineDensityProfile wrapper type definition

static PyTypeObject pyORBIT_InterpolatedLineDensityProfile_Type =
{
  PyObject_HEAD_INIT(NULL)
  0, /*ob_size*/
  "InterpolatedLineDensityProfile", /*tp_name*/
  sizeof(pyORBIT_Object), /*tp_basicsize*/
  0, /*tp_itemsize*/
  (destructor) InterpolatedLineDensityProfile_del , /*tp_dealloc*/
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
  "The InterpolatedLineDensityProfile python wrapper", /* tp_doc */
  0, /* tp_traverse */
  0, /* tp_clear */
  0, /* tp_richcompare */
  0, /* tp_weaklistoffset */
  0, /* tp_iter */
  0, /* tp_iternext */
  InterpolatedLineDensityProfileClassMethods, /* tp_methods */
  InterpolatedLineDensityProfileClassMembers, /* tp_members */
  0, /* tp_getset */
  0, /* tp_base */
  0, /* tp_dict */
  0, /* tp_descr_get */
  0, /* tp_descr_set */
  0, /* tp_dictoffset */
  (initproc) InterpolatedLineDensityProfile_init, /* tp_init */
  0, /* tp_alloc */
  InterpolatedLineDensityProfile_new, /* tp_new */
};


//--------------------------------------------------
// Initialization function of the pyInterpolatedLineDensityProfile class
// It will be called from SpaceCharge wrapper initialization
//--------------------------------------------------

void initInterpolatedLineDensityProfile(PyObject* module)
{
  if (PyType_Ready(&pyORBIT_InterpolatedLineDensityProfile_Type) < 0) return;
  Py_INCREF(&pyORBIT_InterpolatedLineDensityProfile_Type);
  PyModule_AddObject(module, "InterpolatedLineDensityProfile", (PyObject *)&pyORBIT_InterpolatedLineDensityProfile_Type);
  // std::cout << "debug InterpolatedLineDensityProfile added!" << std::endl;
}


#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
