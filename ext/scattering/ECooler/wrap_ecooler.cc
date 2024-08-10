#include "wrap_ecooler.hh"

#include <iostream>

//using namespace OrbitUtils;

namespace wrap_scattering{

#ifdef __cplusplus
extern "C" {
#endif

    //---------------------------------//
    // Python ECooler class definition //
    // ------------------------------- //

    //------------------------------------------------------------------------
    // Member functions for python class wrapping an instance of ECooler
    
    // Constructor for python class. It never will be called directly.
    static PyObject* ECooler_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
    {
        pyORBIT_Object* self;
        self = (pyORBIT_Object *) type->tp_alloc(type, 0);
        self->cpp_obj = NULL;
        // std::cout << "Creating new python ecooler object." << std::endl;
        return (PyObject *) self;
    }

    // Initializer for python ecooler class - Creates the ECooler instance.
    // This is implementation of the __init__ methods ECooler(double damp_x, double damp_y, double damp_z)
    // or ECooler(ElectronBeam *ebeam, double lmax, ECOOLER_MODEL model, double B)
    static int ECooler_init(pyORBIT_Object *self, PyObject *args, PyObject *kwds)
    {
        int Z = 1, nz, nr, nfpara, nfperp; // charge and array sizes
        // Declare the constructor arguments
        PyObject *pyEBeam, *pyvz, *pyvr, *pyfpara, *pyfperp; ECOOLER_MODEL model;
        pyEBeam = pyvz = pyvr = pyfpara = pyfperp = NULL;
        double lmax=0.1, B=0.0, damp_x=0.0, damp_y=0.0, damp_z=0.0, prefactor=1.0, fscale=1.0;
        std::vector<double> vz, vr, fpara, fperp; // Array storage
        // Declare keywords
        static char *kwlist[] = {(char*)"ebeam", (char*)"ecooler_model",
            (char*)"max_length", (char*)"B", (char*)"Z", (char*)"fscale",
            (char*)"damp_x", (char*)"damp_y", (char*)"damp_z", (char*)"prefactor",
            (char*)"vz", (char*)"vr", (char*)"fpara", (char*)"fperp", (char*)"fscale", NULL};
        // Try parsing the input parameters
        if(PyArg_ParseTupleAndKeywords(args, kwds, "Oi|ddidddddOOOO:__init__", kwlist,
                    &pyEBeam, &model, &lmax, &B, &Z, &fscale, &damp_x, &damp_y, &damp_z, &prefactor,
                    &pyvz, &pyvr, &pyfpara, &pyfperp))    
        {
            // The python electronbeam type object
            PyObject *pyORBIT_ElectronBeam_Type=wrap_scattering::getScatteringType("electronbeam");
            ElectronBeam *cpp_ebeam; // Pointer to c++ ElectronBeam class
            if(model != ECOOLER_DUMMY) // Extract the electron beam object for most models
            {
                // Check whether the python object is of the electronbeam type
                if(!PyObject_IsInstance(pyEBeam, pyORBIT_ElectronBeam_Type))
                    ORBIT_MPI_Finalize("scattering.ecooler: ecooler(ebeam,model,...)"
                                " - The first argument should be an electronbeam object.");
                // Extract the electron beam object
                cpp_ebeam = (ElectronBeam*) ((pyORBIT_Object*)pyEBeam)->cpp_obj;
            }
            switch(model) // Invoke different constructors depending on the model.
            {
                case ECOOLER_DUMMY:
                    self->cpp_obj = new ECooler(damp_x, damp_y, damp_z); // Allocate the ECooler object
                    break;
                case ECOOLER_NON_MAGNETIC:
                case ECOOLER_PARKHOMCHUK:
                    self->cpp_obj = new ECooler(cpp_ebeam, lmax, model, B, Z, fscale); // Allocate the ECooler object
                    break;
                case ECOOLER_LUT:
                    // Process list data
                    nz = import_python_list(pyvz, vz);
                    nr = import_python_list(pyvr, vr);
                    nfpara = import_python_list(pyfpara, fpara);
                    nfperp = import_python_list(pyfperp, fperp);
                    // Validate velocity and force data
                    if(nz <= 0 || nr <= 0)
                        ORBIT_MPI_Finalize("scattering.ecooler: ecooler(ebeam, model, vz, vr,...)"
                                " - vz, vr should be a lists of velocities.");
                    if(nfpara != nr*nz || nfperp != nr*nz)
                        ORBIT_MPI_Finalize("scattering.ecooler: ecooler(ebeam,"
                                " model, vz, vr, fpara, fperp,...) - fpara and"
                                " fperp should have len(vz)*len(vr) elements.");
                    // Allocate cooler object
                    self->cpp_obj = new ECooler(cpp_ebeam, prefactor, vz, vr, fpara, fperp, Z);
                    break;
                default:
                    ORBIT_MPI_Finalize("scattering.ecooler: ecooler(ebeam,model,...)"
                            " - Cooler model not recognized.");
            }
        }
        else return 1;
        
            
        ((ECooler*) self->cpp_obj)->setPyWrapper((PyObject*) self); // Point to Python object
        // std::cout << "Created new ECooler instance." << std::endl;
        return 0;
    }

    // Destructor for python ecooler class - Deletes the ECooler instance.
    static void ECooler_del(pyORBIT_Object* self)
    {
        // Get the internal C++ object
        ECooler* cpp_ECooler = (ECooler*) self->cpp_obj;
        // Delete the ECooler instance if present
        if(cpp_ECooler != NULL)
            delete cpp_ECooler;
        // Finally destroy the python object
        self->ob_type->tp_free((PyObject*)self);
        //std::cout << "Destroyed ECooler instance." << std::endl << std::flush;
    }

    // Interface to void ECooler::setUsageELKick(bool flag_elkick)
    static PyObject* ECooler_setUsageELKick(PyObject *self, PyObject *args)
    {
        int flag_elkick=1;
        pyORBIT_Object* pyECooler = (pyORBIT_Object*) self; // Get the python ECooler object
        ECooler* cpp_ECooler = (ECooler*) pyECooler->cpp_obj; // Get the internal ECooler instance
        if(!PyArg_ParseTuple(args,"i:setUsageELKick",&flag_elkick)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("scattering.ecooler: setUsageELKick(flag_elkick) takes 1 argument.");
        cpp_ECooler->setUsageELKick((bool)flag_elkick);
        
        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Interface to void ECooler::setUsageECool(bool flag_ecool)
    static PyObject* ECooler_setUsageECool(PyObject *self, PyObject *args)
    {
        int flag_ecool=1;
        pyORBIT_Object* pyECooler = (pyORBIT_Object*) self; // Get the python ECooler object
        ECooler* cpp_ECooler = (ECooler*) pyECooler->cpp_obj; // Get the internal ECooler instance
        if(!PyArg_ParseTuple(args,"i:setUsageECool",&flag_ecool)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("scattering.ecooler: setUsageECool(flag_ecool) takes 1 argument.");
        cpp_ECooler->setUsageECool((bool)flag_ecool);

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Interface to void ECooler::setForceScaleFactor(double fscale)
    static PyObject* ECooler_setForceScaleFactor(PyObject *self, PyObject *args)
    {
        double fscale=1.0;
        pyORBIT_Object* pyECooler = (pyORBIT_Object*) self; // Get the python ECooler object
        ECooler* cpp_ECooler = (ECooler*) pyECooler->cpp_obj; // Get the internal ECooler instance
        if(!PyArg_ParseTuple(args,"d:setForceScaleFactor",&fscale)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("scattering.ecooler: setForceScaleFactor(fscale) takes 1 argument.");
        cpp_ECooler->setForceScaleFactor(fscale);

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Interface to void ECooler::setErrors(double dE, double dtheta, double dphi)
    static PyObject* ECooler_setErrors(PyObject *self, PyObject *args, PyObject *kwds)
    {
        double dE, dtheta, dphi; // Function arguments
		dE = dtheta = dphi = 0.0; // By default all errors are 0
		// Declare keywords
        static char *kwlist[] = {(char*)"dE", (char*)"dtheta", (char*)"dphi", NULL};
        pyORBIT_Object* pyECooler = (pyORBIT_Object*) self; // Get the python ECooler object
        ECooler* cpp_ECooler = (ECooler*) pyECooler->cpp_obj; // Get the internal ECooler instance
        // Try parsing the input parameters
        if(!PyArg_ParseTupleAndKeywords(args, kwds, "|ddd:setErrors", kwlist,
                    &dE, &dtheta, &dphi))
            ORBIT_MPI_Finalize("scattering.ecooler: setErrors(dE, dtheta, dphi) Could not parse.");

        // Set the errors 
        cpp_ECooler->setErrors(dE, dtheta, dphi);

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Interface to void ECooler::getErrors(double &dE, double &dtheta, double &dphi)
    static PyObject* ECooler_getErrors(PyObject *self, PyObject *args)
    {
        double dE, dtheta, dphi;
        pyORBIT_Object* pyECooler = (pyORBIT_Object*) self; // Get the python ECooler object
        ECooler* cpp_ECooler = (ECooler*) pyECooler->cpp_obj; // Get the internal ECooler instance
        cpp_ECooler->getErrors(dE, dtheta, dphi); // Get the errors
        return Py_BuildValue("ddd", dE, dtheta, dphi); // Return the errors
    }

    // Interface to int ECooler::getcoolingforce(double x, double y, double vx, double vy, double vz,
    //            double &fx, double &fy, double &fz)
    static PyObject* ECooler_getcoolingforce(PyObject *self, PyObject *args)
    {
        double r, vz, vr, fx, fy, fz;
        pyORBIT_Object* pyECooler = (pyORBIT_Object*) self; // Get the python ECooler object
        ECooler* cpp_ECooler = (ECooler*) pyECooler->cpp_obj; // Get the internal ECooler instance
        if(!PyArg_ParseTuple(args,"ddd:getforce",&r,&vz,&vr)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("scattering.ecooler: getcoolingforce(r,vz,vr) takes 3 arguments.");
        cpp_ECooler->getcoolingforce(r, 0.0, vr, 0.0, vz, fx, fy, fz);
        return Py_BuildValue("dd", fz, fx);
    }

    // Interface to ECooler::trackBunch(Bunch* bunch, double length)
    static PyObject* ECooler_trackBunch(PyObject *self, PyObject *args)
    {
        pyORBIT_Object* pyECooler = (pyORBIT_Object*) self; // Get the python ECooler object
        ECooler* cpp_ECooler = (ECooler*) pyECooler->cpp_obj; // Get the internal ECooler instance
        PyObject* pyBunch; // Declare a pointer to the python bunch object
        double length, xc, yc; // We will also need the length through which the particles are tracked.

        // Validate arguments
        if(!PyArg_ParseTuple(args,"Oddd:trackBunch",&pyBunch,&length,&xc,&yc)) // Try to obtain the arguments
            ORBIT_MPI_Finalize("scattering.ecooler: trackBunch(bunch,length,xc,yc) takes 4 arguments.");
        // Obtain the python bunch type
        PyObject* pyORBIT_Bunch_Type = wrap_orbit_bunch::getBunchType("Bunch");
        // Check whether the python object is of the bunch type
        if(!PyObject_IsInstance(pyBunch,pyORBIT_Bunch_Type))
            ORBIT_MPI_Finalize("scattering.ecooler: trackBunch(bunch,length) - The first argument is not a bunch.");
        Bunch* cpp_bunch = (Bunch*) ((pyORBIT_Object*)pyBunch)->cpp_obj;
        cpp_ECooler->trackBunch(cpp_bunch,length,xc,yc); // Finally call the tracker!

        // Return the python object none
        Py_INCREF(Py_None);
        return Py_None;
    }

    // Interface to int ECooler::copyLUT(double &prefactor, std::vector<double> &vz,
    // std::vector<double> &vr, std::vector<double> &fpara, std::vector<double> &fperp)    
    static PyObject* ECooler_copyLUT(PyObject *self, PyObject *args)
    {
        pyORBIT_Object* pyECooler = (pyORBIT_Object*) self; // Get the python ECooler object
        ECooler* cpp_ECooler = (ECooler*) pyECooler->cpp_obj; // Get the internal ECooler instance
        // LUT information
        double prefactor; std::vector<double> vz, vr, fpara, fperp;
        if(cpp_ECooler->copyLUT(prefactor, vz, vr,fpara, fperp))
            return Py_BuildValue("OOOOd", export_python_list(vz), export_python_list(vr),
                    export_python_list(fpara), export_python_list(fperp), prefactor);
        // Return the python object none if we weren't successful
        Py_INCREF(Py_None);
        return Py_None;
    }
    
    
    //------------------------------------------------------------------------
    // scattering.ecooler class manifest - lists methods, members and other
    // properties of the python object
    // Declaration of methods
    static PyMethodDef ECoolerClassMethods[] = {
        {"setUsageELKick", ECooler_setUsageELKick, METH_VARARGS, "Set whether the electron lens kick is applied. - setUsageELKick(flag)"},
        {"setUsageECool", ECooler_setUsageECool, METH_VARARGS, "Set whether the electron cooling force is applied. - setUsageECool(flag)"},
        {"setForceScaleFactor", ECooler_setForceScaleFactor, METH_VARARGS, "Set the scale of the cooler force. 1.0 gives the actual theoretical force. - setForceScaleFactor(fscale)"},
        {"setErrors", (PyCFunction)ECooler_setErrors, METH_VARARGS | METH_KEYWORDS, "Set errors for the electron beam. setErrors(dE, dtheta, dphi)"},
        {"getErrors", ECooler_getErrors, METH_VARARGS, "Get errors for the electron beam. dE, dtheta, dphi = getErrors()"},
        {"getcoolingforce", ECooler_getcoolingforce, METH_VARARGS, "Get the cooling force - fpara, fperp = getforce(vz,vr)"},
        {"trackBunch", ECooler_trackBunch, METH_VARARGS, "Track the bunch - trackBunch(bunch,length)"},
        {"copyLUT", ECooler_copyLUT, METH_VARARGS, "Copy the stored LUT - vz, vr, fpara, fperp = copyLUT()"},
        {NULL}
        //{NULL, NULL}
    };

    // Declaration of members
    static PyMemberDef ECoolerClassMembers [] = {
        {NULL}
    };

    // Definition of PyTypeObject object
    static PyTypeObject pyORBIT_ECooler_Type = {
        PyObject_HEAD_INIT(NULL)
        0, /*ob_size*/
        "ecooler", /*tp_name*/
        sizeof(pyORBIT_Object), /*tp_basicsize*/
        0, /*tp_itemsize*/
        (destructor) ECooler_del , /*tp_dealloc*/
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
        "The ECooler python wrapper", /* tp_doc */
        0, /* tp_traverse */
        0, /* tp_clear */
        0, /* tp_richcompare */
        0, /* tp_weaklistoffset */
        0, /* tp_iter */
        0, /* tp_iternext */
        ECoolerClassMethods, /* tp_methods */
        ECoolerClassMembers, /* tp_members */
        0, /* tp_getset */
        0, /* tp_base */
        0, /* tp_dict */
        0, /* tp_descr_get */
        0, /* tp_descr_set */
        0, /* tp_dictoffset */
        (initproc) ECooler_init, /* tp_init */
        0, /* tp_alloc */
        ECooler_new, /* tp_new */
    };

    // Initialization function of the ECooler class
    // This makes it known to python
    void initECooler(PyObject* module)
    {
        if (PyType_Ready(&pyORBIT_ECooler_Type) < 0) return;
        Py_INCREF(&pyORBIT_ECooler_Type);
        PyModule_AddObject(module, "ecooler", (PyObject *)&pyORBIT_ECooler_Type);
    }

#ifdef __cplusplus
}
#endif

} //end of namespace wrap_scattering

// Some helper functions
// Import python list, flatten multi-dimensional arrays
int import_python_list(PyObject *pylist, std::vector<double> &vec)
{
    if(!pylist)return -1; // Null pointer
    if(!PyList_Check(pylist))return -1; // If the object is not a list
    int i, net_size=0, pylist_size = PyList_Size(pylist); // Find the size of the list
    //std::cout << "(" << pylist_size << ")["; // Debug for now
    for(i=0;i<pylist_size;i++) // Iterate of all elements
    {
        PyObject *pyval = PyList_GetItem(pylist, i); // Get the item
        if(PyFloat_Check(pyval)) // Check whether it's a float
        {
            double val = PyFloat_AsDouble(pyval); // Obtain the C++ double
            vec.push_back(val); // Push this into the array
            net_size++;
            //std::cout << val << ", ";
        }
        else if(PyInt_Check(pyval)) // Check whether it's an int
        {
            int val = PyInt_AsLong(pyval); // Obtain the C++ int
            vec.push_back(val); // Push this into the array
            net_size++;
            //std::cout << val << ", ";
        }
        else if(PyList_Check(pyval))
        {
            // If the item is list, then import values from the list
            int size_of_listitem = import_python_list(pyval, vec);
            if(size_of_listitem <= 0)return -1; // Something went wrong inside
            net_size += size_of_listitem; // Increment the count
            //std::cout << ", ";
        }
        else return -1;
    }
    //std::cout << "]";
    return net_size;
}

// Export double array to a flat python list
PyObject* export_python_list(std::vector<double> &vec)
{
    int i, npts = vec.size(); // Number of points
    PyObject *pylist = PyList_New(npts); // Allocate a new list
    for(i=0; i<npts; i++) // Populate the list
        PyList_SetItem(pylist, i, PyFloat_FromDouble(vec[i]));
    return pylist;
}
