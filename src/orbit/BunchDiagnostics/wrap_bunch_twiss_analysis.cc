#include "orbit_mpi.hh"
#include "pyORBIT_Object.hh"

#include "wrap_bunch_twiss_analysis.hh"
#include "wrap_bunch.hh"

#include <iostream>

#include "BunchTwissAnalysis.hh"

namespace wrap_bunch_twiss_analysis{

  void error(const char* msg){ ORBIT_MPI_Finalize(msg); }

#ifdef __cplusplus
extern "C" {
#endif

	/** 
	    Constructor for python class wrapping c++ BunchTwissAnalysis instance.
      It never will be called directly.
	*/
	static PyObject* BunchTwissAnalysis_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
		pyORBIT_Object* self;
		self = (pyORBIT_Object *) type->tp_alloc(type, 0);
		self->cpp_obj = NULL;
		return (PyObject *) self;
	}
	
  /** This is implementation of the __init__ method */
  static int BunchTwissAnalysis_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds){
		self->cpp_obj =  new BunchTwissAnalysis();
	  ((BunchTwissAnalysis*) self->cpp_obj)->setPyWrapper((PyObject*) self);
    return 0;
  }
  
 /** Performs the Twiss analysis of the bunch */
  static PyObject* BunchTwissAnalysis_analyzeBunch(PyObject *self, PyObject *args){
	    
        int i,j;
        double T[4][4], D[4]; // Transformation matrix and dispersion vector
        bool flag_T, flag_D; // Flags designating whether valid values exist
        
        BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
		// Declare the python objects: transformation matrix, single row of matrix, bunch
        PyObject *pyT = NULL, *pyTrow = NULL, *pyBunch = NULL, *pyD = NULL;
        // Validate arguments
		if(!PyArg_ParseTuple(args,"O|O!O!:analyzeBunch",&pyBunch, &PyList_Type, &pyT, &PyList_Type, &pyD)){
			ORBIT_MPI_Finalize("BunchTwissAnalysis - analyzeBunch(bunch, T, D) - Bunch is needed for analysis. T if provided must be a 4x4 python list object. D if provided must be a list with 4 doubles.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("BunchTwissAnalysis - analyzeBunch(bunch, T, D) - bunch must be a bunch object.");
		}
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;

        // Copy transformation matrix
        if(pyT && PyList_Size(pyT) == 4) // The number of rows should be 4
        {
            for(i=0;i<4;i++)
            {
                pyTrow = PyList_GetItem(pyT, i);
                for(j=0;j<4;j++)
                {
                    T[i][j] = PyFloat_AS_DOUBLE(PyList_GetItem(pyTrow, j));
                }
            }
            flag_T = true;
        }
        else flag_T = false;
        // Copy dispersion vector
        if(pyD && PyList_Size(pyD) == 4) // The number of rows should be 4
        {
            for(i=0;i<4;i++)
                D[i] = PyFloat_AS_DOUBLE(PyList_GetItem(pyD, i));
            flag_D = true;
        }
        else flag_D = false;

        //std::cout << "flag_D = " << flag_D << " flag_T = " << flag_T << std::endl;

        // Finally call the analysis function
        if(flag_T && flag_D)
            cpp_BunchTwissAnalysis->analyzeBunch(cpp_bunch, T, D);
        else if(flag_T)
            cpp_BunchTwissAnalysis->analyzeBunch(cpp_bunch, T, NULL);
        else if(flag_D)
            cpp_BunchTwissAnalysis->analyzeBunch(cpp_bunch, NULL, D);
        else
            cpp_BunchTwissAnalysis->analyzeBunch(cpp_bunch);


		Py_INCREF(Py_None);
		return Py_None;
  }
	
	/** Returns the XY moments of the beam up to a prescribed order */
	static PyObject* BunchTwissAnalysis_computeBunchMoments(PyObject *self, PyObject *args){
		BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
		PyObject* pyBunch;
		int order;
		int dispterm;
		int emitnormterm;
		if(!PyArg_ParseTuple(args,"Oiii:computeBunchMoments",&pyBunch,&order,&dispterm,&emitnormterm)){
			ORBIT_MPI_Finalize("BunchTwissAnalysis - computeBunchMoments(Bunch* bunch, int order, int dispterm, int emitnormterm) - parameters are needed.");
		}
		PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
		if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type)){
			ORBIT_MPI_Finalize("BunchTwissAnalysis - computeBunchMoments(Bunch* bunch, int order, int dispterm, int emitnormterm) - method needs a Bunch.");
		}
		
		Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
		cpp_BunchTwissAnalysis->computeBunchMoments(cpp_bunch, order, dispterm, emitnormterm);
		Py_INCREF(Py_None);
		return Py_None;
	
	}
		
	/** It will return the centered correlation <(x-<x>)*(y-<y>)> = <x*y> - <x>*<y> for coordinates with indeces (ic,jc) */
	static PyObject* BunchTwissAnalysis_getCorrelation(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic,jc;
		if(!PyArg_ParseTuple(	args,"ii:getCorrelation",&ic,&jc)){
			error("pyBunchTwissAnalysis.getCorrelation(ic,jc) - parameters are needed");
		}
    return Py_BuildValue("d",cpp_BunchTwissAnalysis->getCorrelation(ic,jc));
	}
	
	/** It will return the (i,j) XY moment of the beam */
	static PyObject* BunchTwissAnalysis_getBunchMoment(PyObject *self, PyObject *args){
		BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
		int i,j;
		if(!PyArg_ParseTuple(	args,"ii:getBunchMoment",&i,&j)){
			error("pyBunchTwissAnalysis.getBunchMoment(i,j) - parameters are needed");
		}
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getBunchMoment(i,j));
	}

	/** Returns the average value for coordinate with index ic */
  static PyObject* BunchTwissAnalysis_getAverage(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:",&ic)){
			error("pyBunchTwissAnalysis.getAverage(ic) - parameter is needed");
		}
    return Py_BuildValue("d",cpp_BunchTwissAnalysis->getAverage(ic));
  }		
	
	/** Returns the total number of analysed macroparticles */
  static PyObject* BunchTwissAnalysis_getGlobalCount(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
    return Py_BuildValue("i",cpp_BunchTwissAnalysis->getGlobalCount());
  }		
		
	/** Returns the total macrosize */
  static PyObject* BunchTwissAnalysis_getGlobalMacrosize(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
    return Py_BuildValue("d",cpp_BunchTwissAnalysis->getGlobalMacrosize());
  }			
	
	//------------------------------------------------------------
	//Twiss functions
	//------------------------------------------------------------
 	/** It returns the emittance for index 0,1,2 - x,y,z planes*/
  static PyObject* BunchTwissAnalysis_getEmittance(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getEmittance",&ic)){
			error("pyBunchTwissAnalysis.getEmittance(ic) - parameter is needed");
		}		
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getEmittance(ic));
  }		
	
 	/** It returns the normalized emittance for index 0,1 - x,y planes*/
  static PyObject* BunchTwissAnalysis_getEmittanceNormalized(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getEmittanceNormalized",&ic)){
			error("pyBunchTwissAnalysis.getEmittanceNormalized(ic) - parameter is needed");
		}		
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getEmittanceNormalized(ic));
  }		
	
 	/** It returns the Twiss alpha for index 0,1,2 - x,y,z planes*/
  static PyObject* BunchTwissAnalysis_getAlpha(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getAlpha",&ic)){
			error("pyBunchTwissAnalysis.getAlpha(ic) - parameter is needed");
		}				
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getAlpha(ic));
  }		
	
 	/** It returns the Twiss beta for index 0,1,2 - x,y,z planes*/
  static PyObject* BunchTwissAnalysis_getBeta(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getBeta",&ic)){
			error("pyBunchTwissAnalysis.getBeta(ic) - parameter is needed");
		}				
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getBeta(ic));
  }		
	
 	/** It returns the Twiss gamma for index 0,1,2 - x,y,z planes*/
  static PyObject* BunchTwissAnalysis_getGamma(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getGamma",&ic)){
			error("pyBunchTwissAnalysis.getGamma(ic) - parameter is needed");
		}					
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getGamma(ic));
  }		

 	/** It returns Twiss dispersion function for index 0,1 - x,y planes*/
  static PyObject* BunchTwissAnalysis_getDispersion(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getDispersion",&ic)){
			error("pyBunchTwissAnalysis.getDispersion(ic) - parameter is needed");
		}					
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getDispersion(ic));
  }		
	
 	/** It returns Twiss dispersion_prime function for index 0,1 - x,y planes*/
  static PyObject* BunchTwissAnalysis_getDispersionDerivative(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getDispersionDerivative",&ic)){
			error("pyBunchTwissAnalysis.getDispersionDerivative(ic) - parameter is needed");
		}					
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getDispersionDerivative(ic));
  }		

 	/** It returns the Twiss array (alpha,beta,gamma,emittance) for index 0,1,2 - x,y,z planes*/
  static PyObject* BunchTwissAnalysis_getTwiss(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getTwiss",&ic)){
			error("pyBunchTwissAnalysis.getTwiss(ic) - parameter is needed");
		}		
		double alpha = cpp_BunchTwissAnalysis->getAlpha(ic);
		double beta = cpp_BunchTwissAnalysis->getBeta(ic);
		double gamma = cpp_BunchTwissAnalysis->getGamma(ic);
		double emitt = cpp_BunchTwissAnalysis->getEmittance(ic);
		return Py_BuildValue("(dddd)",alpha,beta,gamma,emitt);
  }			
	
	
 	/** It returns the effective emittance for index 0,1 - x,y planes*/
  static PyObject* BunchTwissAnalysis_getEffectiveEmittance(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getEffectiveEmittance",&ic)){
			error("pyBunchTwissAnalysis.getEffectiveEmittance(ic) - parameter is needed");
		}		
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getEffectiveEmittance(ic));
  }		
		
 	/** It returns the effective Twiss alpha for index 0,1 - x,y planes*/
  static PyObject* BunchTwissAnalysis_getEffectiveAlpha(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getEffectiveAlpha",&ic)){
			error("pyBunchTwissAnalysis.getEffectiveAlpha(ic) - parameter is needed");
		}				
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getEffectiveAlpha(ic));
  }		
	
 	/** It returns the effective Twiss beta for index 0,1 - x,y planes*/
  static PyObject* BunchTwissAnalysis_getEffectiveBeta(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getEffectiveBeta",&ic)){
			error("pyBunchTwissAnalysis.getEffectiveBeta(ic) - parameter is needed");
		}				
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getEffectiveBeta(ic));
  }		
	
 	/** It returns the effective Twiss gamma for index 0,1 - x,y planes*/
  static PyObject* BunchTwissAnalysis_getEffectiveGamma(PyObject *self, PyObject *args){
	  BunchTwissAnalysis* cpp_BunchTwissAnalysis = (BunchTwissAnalysis*)((pyORBIT_Object*) self)->cpp_obj;
	  int ic;
		if(!PyArg_ParseTuple(	args,"i:getEffectiveGamma",&ic)){
			error("pyBunchTwissAnalysis.getEffectiveGamma(ic) - parameter is needed");
		}					
		return Py_BuildValue("d",cpp_BunchTwissAnalysis->getEffectiveGamma(ic));
  }		
	
  //-----------------------------------------------------
  //destructor for python BunchTwissAnalysis class (__del__ method).
  //-----------------------------------------------------
  static void BunchTwissAnalysis_del(pyORBIT_Object* self){
		//std::cerr<<"The BunchTwissAnalysis __del__ has been called!"<<std::endl;
		delete ((BunchTwissAnalysis*)self->cpp_obj);
		self->ob_type->tp_free((PyObject*)self);
  }
	
	// defenition of the methods of the python BunchTwissAnalysis wrapper class
	// they will be vailable from python level
  static PyMethodDef BunchTwissAnalysisClassMethods[] = {
		{ "analyzeBunch",			BunchTwissAnalysis_analyzeBunch,			METH_VARARGS,"Performs the Twiss analysis of the bunch."},
		{ "computeBunchMoments",	BunchTwissAnalysis_computeBunchMoments,		METH_VARARGS,"Returns the XY moments of the beam up to a prescribed order"},
 		{ "getCorrelation",			BunchTwissAnalysis_getCorrelation,    		METH_VARARGS,"Returns the centered correlation <(x-<x>)*(y-<y>)> = <x*y> - <x>*<y> for coordinates with indeces (ic,jc)"},
		{ "getBunchMoment",			BunchTwissAnalysis_getBunchMoment,    		METH_VARARGS,"Returns the (i,j) xy moment of the beam"},
 		{ "getAverage",			 	BunchTwissAnalysis_getAverage,       		METH_VARARGS,"Returns the average value for coordinate with index ic"},		
 		{ "getGlobalCount",			BunchTwissAnalysis_getGlobalCount,    		METH_VARARGS,"Returns the total number of analysed macroparticles"},		
 		{ "getGlobalMacrosize", 	BunchTwissAnalysis_getGlobalMacrosize,		METH_VARARGS,"Returns the total macrosize"},		
		{ "getEmittance",			BunchTwissAnalysis_getEmittance,   			METH_VARARGS,"Returns the emittance for index 0,1,2 - x,y,z planes"},		
		{ "getEmittanceNormalized",	BunchTwissAnalysis_getEmittanceNormalized,	METH_VARARGS,"Returns the normalized emittance for index 0,1 - x,y planes"},		
		{ "getAlpha",			 	BunchTwissAnalysis_getAlpha,    			METH_VARARGS,"Returns Twiss alpha for index 0,1,2 - x,y,z planes"},		
		{ "getBeta",			 	BunchTwissAnalysis_getBeta,    				METH_VARARGS,"Returns Twiss beta for index 0,1,2 - x,y,z planes"},		
		{ "getGamma",			 	BunchTwissAnalysis_getGamma,    			METH_VARARGS,"Returns Twiss gamma for index 0,1,2 - x,y,z planes"},				
		{ "getTwiss",				BunchTwissAnalysis_getTwiss,    			METH_VARARGS,"Returns Twiss tuple (alpha,beta,gamma,emitt) for index 0,1,2 - x,y,z planes"},
		{ "getDispersion",			BunchTwissAnalysis_getDispersion,			METH_VARARGS,"Returns Twiss dispersion function for index 0,1 - x,y planes"},				
		{ "getDispersionDerivative",BunchTwissAnalysis_getDispersionDerivative,	METH_VARARGS,"Returns Twiss dispersion' function for index 0,1 - x,y planes"},				
		{ "getEffectiveEmittance",	BunchTwissAnalysis_getEffectiveEmittance,	METH_VARARGS,"Returns the effective emittance for index 0,1 - x,y planes"},		
		{ "getEffectiveAlpha",		BunchTwissAnalysis_getEffectiveAlpha,    	METH_VARARGS,"Returns effective Twiss alpha for index 0,1 - x,y planes"},		
		{ "getEffectiveBeta",		BunchTwissAnalysis_getEffectiveBeta,    	METH_VARARGS,"Returns effective Twiss beta for index 0,1 - x,y planes"},		
		{ "getEffectiveGamma",		BunchTwissAnalysis_getEffectiveGamma,		METH_VARARGS,"Returns effective Twiss gamma for index 0,1 - x,y planes"},				
		{NULL}
  };
	
	// defenition of the memebers of the python BunchTwissAnalysis wrapper class
	// they will be vailable from python level
	static PyMemberDef BunchTwissAnalysisClassMembers [] = {
		{NULL}
	};
	
	//new python BunchTwissAnalysis wrapper type definition
	static PyTypeObject pyORBIT_BunchTwissAnalysis_Type = {
		PyObject_HEAD_INIT(NULL)
		0, /*ob_size*/
		"BunchTwissAnalysis", /*tp_name*/
		sizeof(pyORBIT_Object), /*tp_basicsize*/
		0, /*tp_itemsize*/
		(destructor) BunchTwissAnalysis_del , /*tp_dealloc*/
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
		"The BunchTwissAnalysis python wrapper", /* tp_doc */
		0, /* tp_traverse */
		0, /* tp_clear */
		0, /* tp_richcompare */
		0, /* tp_weaklistoffset */
		0, /* tp_iter */
		0, /* tp_iternext */
		BunchTwissAnalysisClassMethods, /* tp_methods */
		BunchTwissAnalysisClassMembers, /* tp_members */
		0, /* tp_getset */
		0, /* tp_base */
		0, /* tp_dict */
		0, /* tp_descr_get */
		0, /* tp_descr_set */
		0, /* tp_dictoffset */
		(initproc) BunchTwissAnalysis_init, /* tp_init */
		0, /* tp_alloc */
		BunchTwissAnalysis_new, /* tp_new */
	};	
	
	
	
	//--------------------------------------------------
	//Initialization BunchTwissAnalysis of the pyBunchTwissAnalysis class
	//--------------------------------------------------
  void initbunchtwissanalysis(PyObject* module){
		if (PyType_Ready(&pyORBIT_BunchTwissAnalysis_Type) < 0) return;
		Py_INCREF(&pyORBIT_BunchTwissAnalysis_Type);
		PyModule_AddObject(module, "BunchTwissAnalysis", (PyObject *)&pyORBIT_BunchTwissAnalysis_Type);
	}

#ifdef __cplusplus
}
#endif


}
