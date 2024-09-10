#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"
#
#include "wrap_spacechargecalc_analytic_gaussian.hh"
#include "wrap_spacecharge.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "SpaceChargeCalcAnalyticGaussian.hh"

using namespace OrbitUtils;

namespace wrap_spacecharge{

#ifdef __cplusplus
extern "C" {
#endif

	//---------------------------------------------------------
	//Python SpaceChargeCalcAnalyticGaussian class definition
	//---------------------------------------------------------

	//constructor for python class wrapping SpaceChargeForceCalc2p5D instance
	//It never will be called directly

	static PyObject* SpaceChargeCalcAnalyticGaussian_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
	{
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		//std::cerr<<"The SpaceChargeCalcAnalyticGaussian new has been called!"<<std::endl;
		return (PyObject *) self;
	}
	
  //initializator for python SpaceChargeCalcAnalyticGaussian class
  //this is implementation of the __init__ method SpaceChargeCalcAnalyticGaussian(double intensity, double epsn_x, double epsn_y, double dpp_rms)
  static int SpaceChargeCalcAnalyticGaussian_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
  	double intensity, epsn_x, epsn_y, dpp_rms;
  	pyORBIT_Object* pyLineDensityProfile;
  	
		if(!PyArg_ParseTuple(args,"ddddO:__init__", &intensity, &epsn_x, &epsn_y, &dpp_rms, &pyLineDensityProfile)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcAnalyticGaussian - SpaceChargeCalcAnalyticGaussian(intensity, epsn_x, epsn_y, dpp_rms, LineDensityProfile) - constructor needs parameters.");
		}
		LineDensityProfile* cpp_LineDensityProfile = (LineDensityProfile*) ((pyORBIT_Object*)pyLineDensityProfile)->cpp_obj; 
		self->cpp_obj = new SpaceChargeCalcAnalyticGaussian(intensity, epsn_x, epsn_y, dpp_rms, cpp_LineDensityProfile);
		//std::cerr<<"The SpaceChargeCalcAnalyticGaussian __init__ has been called!"<<std::endl;
		return 0;
	}
	

  //setBunchParameters(double intensity, double epsn_x, double epsn_y, double dpp_rms)
  static PyObject* SpaceChargeCalcAnalyticGaussian_setBunchParameters(PyObject *self, PyObject *args){
		int nVars = PyTuple_Size(args);
		pyORBIT_Object* pySpaceChargeCalcAnalyticGaussian = (pyORBIT_Object*) self;
		SpaceChargeCalcAnalyticGaussian* cpp_SpaceChargeCalcAnalyticGaussian = (SpaceChargeCalcAnalyticGaussian*) pySpaceChargeCalcAnalyticGaussian->cpp_obj;
		double intensity, epsn_x, epsn_y, dpp_rms;
		
		if(!PyArg_ParseTuple(args,"dddd:setBunchParameters", &intensity, &epsn_x, &epsn_y, &dpp_rms)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcAnalyticGaussian.setBunchParameters(intensity, epsn_x, epsn_y, dpp_rms) - parameters are needed.");
		}
		cpp_SpaceChargeCalcAnalyticGaussian->setBunchParameters(intensity, epsn_x, epsn_y, dpp_rms);
		Py_INCREF(Py_None);
		return Py_None;  
  }


  //setLineDensityProfile(LineDensityProfile* LineDensityProfile)
  static PyObject* SpaceChargeCalcAnalyticGaussian_setLineDensityProfile(PyObject *self, PyObject *args){
		pyORBIT_Object* pySpaceChargeCalcAnalyticGaussian = (pyORBIT_Object*) self;
		SpaceChargeCalcAnalyticGaussian* cpp_SpaceChargeCalcAnalyticGaussian = (SpaceChargeCalcAnalyticGaussian*) pySpaceChargeCalcAnalyticGaussian->cpp_obj;
		pyORBIT_Object* pyLineDensityProfile;
		if(!PyArg_ParseTuple(args,"O:setLineDensityProfile", &pyLineDensityProfile)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcAnalyticGaussian.setLineDensityProfile(LineDensityProfile) - parameters are needed.");
		}
		LineDensityProfile* cpp_LineDensityProfile = (LineDensityProfile*) ((pyORBIT_Object*)pyLineDensityProfile)->cpp_obj;
		cpp_SpaceChargeCalcAnalyticGaussian->setLineDensityProfile(cpp_LineDensityProfile);
		Py_INCREF(Py_None);
		return Py_None;  
  }
  
	
  //trackBunch(Bunch* bunch, double length, [BaseBoundary2D* boundary])
  static PyObject* SpaceChargeCalcAnalyticGaussian_trackBunch(PyObject *self, PyObject *args){
		int nVars = PyTuple_Size(args);
		pyORBIT_Object* pySpaceChargeCalcAnalyticGaussian = (pyORBIT_Object*) self;
		SpaceChargeCalcAnalyticGaussian* cpp_SpaceChargeCalcAnalyticGaussian = (SpaceChargeCalcAnalyticGaussian*) pySpaceChargeCalcAnalyticGaussian->cpp_obj;
		PyObject* pyBunch;
		double length;
		
		if(!PyArg_ParseTuple(args,"Od:trackBunch",&pyBunch,&length)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcAnalyticGaussian.trackBunch(pyBunch,length) - parameters are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcAnalyticGaussian.trackBunch(pyBunch,length) - pyBunch is not Bunch.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_SpaceChargeCalcAnalyticGaussian->trackBunch(cpp_bunch,length);
		Py_INCREF(Py_None);
		return Py_None;  
  }
  
  
  //setLatticeParameters(double beta_x, double beta_y, double eta_x)
  static PyObject* SpaceChargeCalcAnalyticGaussian_setLatticeParameters(PyObject *self, PyObject *args){
  		pyORBIT_Object* pySpaceChargeCalcAnalyticGaussian = (pyORBIT_Object*) self;
		SpaceChargeCalcAnalyticGaussian* cpp_SpaceChargeCalcAnalyticGaussian = (SpaceChargeCalcAnalyticGaussian*) pySpaceChargeCalcAnalyticGaussian->cpp_obj;
		double beta_x, beta_y, eta_x, eta_y, co_x, co_y;
  		
  		if(!PyArg_ParseTuple(args,"dddddd:setLatticeParameters", &beta_x, &beta_y, &eta_x, &eta_y, &co_x, &co_y)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcAnalyticGaussian - setLatticeParameters(beta_x, beta_y, eta_x, eta_y, co_x, co_y) - parameters are needed.");
  		}
  		cpp_SpaceChargeCalcAnalyticGaussian->setLatticeParameters(beta_x, beta_y, eta_x, eta_y, co_x, co_y);
  		Py_INCREF(Py_None);
  		return Py_None;
  }
		
		
  //getLineDensityFactor(double z)
  static PyObject* SpaceChargeCalcAnalyticGaussian_getLineDensityFactor(PyObject *self, PyObject *args){
  		pyORBIT_Object* pySpaceChargeCalcAnalyticGaussian = (pyORBIT_Object*) self;
		SpaceChargeCalcAnalyticGaussian* cpp_SpaceChargeCalcAnalyticGaussian = (SpaceChargeCalcAnalyticGaussian*) pySpaceChargeCalcAnalyticGaussian->cpp_obj;
		double z;
  		if(!PyArg_ParseTuple(args,"d:getLineDensityFactor", &z)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcAnalyticGaussian - getLineDensityFactor(z) - parameters are needed.");
  		}
  		double localLineDensityFactor = cpp_SpaceChargeCalcAnalyticGaussian->getLineDensityFactor(z);
  		return Py_BuildValue("d", localLineDensityFactor);
  }  


  //BassettiErskine(double x, double y, double sigma_x, double sigma_y)
  static PyObject* SpaceChargeCalcAnalyticGaussian_BassettiErskine(PyObject *self, PyObject *args){
  		pyORBIT_Object* pySpaceChargeCalcAnalyticGaussian = (pyORBIT_Object*) self;
		SpaceChargeCalcAnalyticGaussian* cpp_SpaceChargeCalcAnalyticGaussian = (SpaceChargeCalcAnalyticGaussian*) pySpaceChargeCalcAnalyticGaussian->cpp_obj;
		double x, y, sigma_x, sigma_y, E_x, E_y;
  		
  		if(!PyArg_ParseTuple(args,"dddd:BassettiErskine", &x, &y, &sigma_x, &sigma_y)){
			ORBIT_MPI_Finalize("PySpaceChargeCalcAnalyticGaussian - BassettiErskine(x, y, sigma_x, sigma_y) - parameters are needed.");
  		}
  		cpp_SpaceChargeCalcAnalyticGaussian->BassettiErskine(x, y, sigma_x, sigma_y, E_x, E_y);
  		return Py_BuildValue("dd", E_x, E_y);
  }

  
  //-----------------------------------------------------
  //destructor for python SpaceChargeCalcAnalyticGaussian class (__del__ method).
  //-----------------------------------------------------
  static void SpaceChargeCalcAnalyticGaussian_del(pyORBIT_Object* self){
		SpaceChargeCalcAnalyticGaussian* cpp_SpaceChargeCalcAnalyticGaussian = (SpaceChargeCalcAnalyticGaussian*) self->cpp_obj;
		if(cpp_SpaceChargeCalcAnalyticGaussian != NULL){
			delete cpp_SpaceChargeCalcAnalyticGaussian;
		}
		self->ob_type->tp_free((PyObject*)self);
  }	
  
  // defenition of the methods of the python SpaceChargeCalcAnalyticGaussian wrapper class
  // they will be vailable from python level
  static PyMethodDef SpaceChargeCalcAnalyticGaussianClassMethods[] = {
		{ "setBunchParameters",    SpaceChargeCalcAnalyticGaussian_setBunchParameters,    METH_VARARGS,"set the bunch parameters - setBunchParameters(intensity, epsn_x, epsn_y, dpp_rms)"},
		{ "setLineDensityProfile", SpaceChargeCalcAnalyticGaussian_setLineDensityProfile, METH_VARARGS,"set the LineDensityProfile setLineDensityProfile(LineDensityProfile)"},
		{ "trackBunch",            SpaceChargeCalcAnalyticGaussian_trackBunch,            METH_VARARGS,"track the bunch - trackBunch(pyBunch, length, boundary)"},
		{ "setLatticeParameters",  SpaceChargeCalcAnalyticGaussian_setLatticeParameters,  METH_VARARGS,"set the lattice parameters for a space charge node - setLatticeParameters(beta_x, beta_y, eta_x)"},
		{ "getLineDensityFactor",  SpaceChargeCalcAnalyticGaussian_getLineDensityFactor,  METH_VARARGS,"get the local line density factor at z - getLineDensityProfile(z)"},
		{ "BassettiErskine",       SpaceChargeCalcAnalyticGaussian_BassettiErskine,       METH_VARARGS,"evaluate the Bassetti Erskine formula - BassettiErskine(x, y, sigma_x, sigma_y)"},
		{NULL}
  };
  
  // defenition of the memebers of the python SpaceChargeCalcAnalyticGaussian wrapper class
  // they will be vailable from python level
  static PyMemberDef SpaceChargeCalcAnalyticGaussianClassMembers [] = {
		{NULL}
  };

	//new python SpaceChargeCalcAnalyticGaussian wrapper type definition
	static PyTypeObject pyORBIT_SpaceChargeCalcAnalyticGaussian_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"SpaceChargeCalcAnalyticGaussian", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) SpaceChargeCalcAnalyticGaussian_del , /*tp_dealloc*/
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
		"The SpaceChargeCalcAnalyticGaussian python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		SpaceChargeCalcAnalyticGaussianClassMethods, /* tp_methods */
		SpaceChargeCalcAnalyticGaussianClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) SpaceChargeCalcAnalyticGaussian_init, /* tp_init */
		0, /* tp_alloc */
		SpaceChargeCalcAnalyticGaussian_new, /* tp_new */
	};	

	//--------------------------------------------------
	//Initialization function of the pySpaceChargeCalcAnalyticGaussian class
	//It will be called from SpaceCharge wrapper initialization
	//--------------------------------------------------
  void initSpaceChargeCalcAnalyticGaussian(PyObject* module){
		if (PyType_Ready(&pyORBIT_SpaceChargeCalcAnalyticGaussian_Type) < 0) return;
		Py_INCREF(&pyORBIT_SpaceChargeCalcAnalyticGaussian_Type);
		PyModule_AddObject(module, "SpaceChargeCalcAnalyticGaussian", (PyObject *)&pyORBIT_SpaceChargeCalcAnalyticGaussian_Type);
	}

#ifdef __cplusplus
}
#endif

//end of namespace wrap_spacecharge
}
