/*****************************************************************************/
/* File: ndim_integrator.hh
 * Description: This file declares a class to calculate integrals of defined
 * functions over n dimensions. The class is most useful when the integral
 * needs to be evaulated repeatedly in a program. It uses gsl_integration
 * routines to compute the 1D integrals.
 */
/*****************************************************************************/

#ifndef NDIM_INTEGRATOR_H
#define NDIM_INTEGRATOR_H

// Includes
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

// Some constant definitions
#define GSL_WORKSPACE_SIZE 1000 // Maximum size of gsl workspaces
#define INTEGRAND_DUMP_FILE "ndim_integrand.dump"
#define INTEGRAND_DUMP_NPTS 100

/*****************************************************************************/
// Define some useful types
typedef double (*user_func_ptr)(double*,void*);
typedef int (*user_singularity_ptr)(double*,void*,double*);
typedef std::vector<double> dblvec;

// Define some constants to decide integrator type
typedef enum {NDIM_INT_QAG, NDIM_INT_QAGS, NDIM_INT_QAGP} NDIM_INT_TYPE;

// Integrand function wrapper
double integrand_function_wrapper(double x, void *ndim_params);

// Declare n dimensional integrator class
class ndim_integrator;

/*****************************************************************************/
// Parameter structure for the integrand function wrapper
typedef struct
{
    // Pointer to the integrator object
    ndim_integrator *integrator;
    // Which variable is being integrated?
    int integral_number;
} ndim_integrator_params;

/*****************************************************************************/
// Definition of integrator class
class ndim_integrator
{
    int ndim, verbosity; // Number of dimensions of the integral

    // Relative tolerance of integral
    double relative_tolerance;

    // Integrator object name
    std::string integrator_name;

    // Vector of variables
    dblvec all_variables;
    dblvec all_lower_limits;
    dblvec all_upper_limits;

    // Integrator types
    std::vector<int> gsl_integrator_types;
    // Pointer to user specified params
    void *user_specified_params;
    // Vector of user specified integrands
    std::vector<user_func_ptr> user_specified_integrands;
    // Vector of singularity calculation functions
    std::vector<user_singularity_ptr> user_specified_singularities;

    // Parameter structs to pass to the wrapper function 
    std::vector<ndim_integrator_params> gsl_parameters;
    // GSL workspaces
    std::vector<gsl_integration_workspace*> gsl_workspaces;
    // GSL functions
    std::vector<gsl_function> gsl_functions;

    public:

    // Constructor sets number of variables and integration limits
    ndim_integrator(int nvars, dblvec ll, dblvec ul, double epsrel = 1e-7);
    // Sets each partial integrand
    int set_integrand(int integral_number, user_func_ptr func);
    // Sets each integration type
    int set_integrator_type(int integral_number, NDIM_INT_TYPE integral_type);
    // Sets the functions used to calculate singularities
    int set_singularity(int integral_number, user_singularity_ptr singularity_func);
    // Set the relative tolerance
    void set_reltol(double epsrel);
    // Set verbosity
    void set_verbosity(int v);
    // Set integrator name
    void set_name(std::string name);
    // Dump integrand
    void dump_integrand(int integral_number);
    // Evaluates an integrand at a given point
    double evaluate_integrand(double x, int integral_number);
    // Evaluates the integral for given input parameters
    double evaluate(void *input_params);
    // Destructor. Clears all allocations
    ~ndim_integrator();

};

#endif // End of header
