/*****************************************************************************/
/* File: ndim_integrator.cc
 * Description: This file defines a class to calculate integrals of defined
 * functions over n dimensions. The class is most useful when the integral
 * needs to be evaulated repeatedly in a program. It uses gsl_integration
 * routines to compute the 1D integrals.
 */
/*****************************************************************************/

// Includes
#include "ndim_integrator.hh"

// Integrand function wrapper
double integrand_function_wrapper(double x, void *ndim_params)
{
    // Extract the struct
    ndim_integrator_params *params_struct = (ndim_integrator_params*)ndim_params;
    
    // Call the evaluation method of the integrator object
    return params_struct->integrator->evaluate_integrand(x,
            params_struct->integral_number);
}

/*****************************************************************************/
// ndim_integrator class methods
// Constructor sets number of variables and integration limits
ndim_integrator::ndim_integrator(int nvars, dblvec ll, dblvec ul, double epsrel)
{
    int i; // Use for looping
    ndim = nvars; // Set the number of dimensions of the integral
    // Assign limits
    all_lower_limits = ll;
    all_upper_limits = ul;
    // Check whether we have valid arguments.
    if(ll.size() != ndim || ul.size() != ndim)
        throw std::length_error("ndim_integrator::ndim_integrator(int nvars,"
            " dblvec ll, dblvec ul) The length of the vectors ll and ul"
            " should equal nvars.");

    // Allocate space for all vectors
    all_variables.resize(ndim);
    gsl_integrator_types.resize(ndim);
    user_specified_integrands.resize(ndim);
    user_specified_singularities.resize(ndim);
    gsl_parameters.resize(ndim);
    gsl_workspaces.resize(ndim);
    gsl_functions.resize(ndim);

    // Set up parameters, workspaces and functions
    for(i=0; i<ndim; i++)
    {
        all_variables[i] = 0.0; // Initialize all variables to zero
        gsl_parameters[i].integrator = this; // Pointer to this object
        gsl_parameters[i].integral_number = i; // Set the integral number
        gsl_functions[i].function = integrand_function_wrapper; // Pointer to wrapper
        gsl_functions[i].params = (void*)(&gsl_parameters[i]); // Set the params
        // Allocate workspace
        gsl_workspaces[i] = gsl_integration_workspace_alloc(GSL_WORKSPACE_SIZE);
        // Set interand functions and singularities to NULL
        user_specified_integrands[i] = NULL;
        user_specified_singularities[i] = NULL;
        // By default use the QAGS integrator
        gsl_integrator_types[i] = NDIM_INT_QAGS;
    }

    // Set some defaults
    relative_tolerance = epsrel; // Tolerance
    verbosity = 0; // GSL will call an abort if there's a problem.
    integrator_name = "NO_NAME";

}

// Destructor clears all allocated memory
ndim_integrator::~ndim_integrator()
{
    // All the std::vector objects are freed automatically.
    int i;
    // Free the workspace memory
    for(i=0;i<ndim;i++)
        if(gsl_workspaces[i])
            gsl_integration_workspace_free(gsl_workspaces[i]);
}

// Function to set an integrand
// integral_number = 0 (outermost) to ndim - 1 (innermost)
int ndim_integrator::set_integrand(int integral_number, user_func_ptr func)
{
    if(integral_number >= 0 && integral_number < ndim) // Check validity
    {
        // Assign the function pointer to the integrand
        user_specified_integrands[integral_number] = func;
        return GSL_SUCCESS;
    }
    else return GSL_FAILURE;
}

// Function to set integrator type
// integral_number = 0 (outermost) to ndim - 1 (innermost)
int ndim_integrator::set_integrator_type(int integral_number, NDIM_INT_TYPE integral_type)
{
    if(integral_number >= 0 && integral_number < ndim) // Check validity
    {
        // Assign the integrator type
        gsl_integrator_types[integral_number] = integral_type;
        return GSL_SUCCESS;
    }
    else return GSL_FAILURE;
}

// Function to set singularities
// integral_number = 0 (outermost) to ndim - 1 (innermost)
int ndim_integrator::set_singularity(int integral_number, user_singularity_ptr singularity_func)
{
    if(integral_number >= 0 && integral_number < ndim) // Check validity
    {   
        // Assign the function pointer to the integrand
        user_specified_singularities[integral_number] = singularity_func;
        return GSL_SUCCESS;
    }
    else return GSL_FAILURE;
}

// Set the relative tolerance
void ndim_integrator::set_reltol(double epsrel)
{
    relative_tolerance = epsrel; // Set the relative tolerance
}

void ndim_integrator::set_verbosity(int v)
{
    verbosity = v;
}

void ndim_integrator::set_name(std::string name)
{
    integrator_name = name;
}

// Run the integrator
double ndim_integrator::evaluate_integrand(double x, int integral_number)
{
    double integrand = 1.0; // Initial value of integrand
    int success = GSL_SUCCESS, next_integral = integral_number + 1; // Next integral
    if(integral_number > -1)
    {
        all_variables[integral_number] = x; // Store current value
        // Get pointer to user specified integrand
        user_func_ptr intfunc = user_specified_integrands[integral_number];
        // Evaluate the integrand if a function is specified.
        if(intfunc) integrand *= intfunc(all_variables.data(), user_specified_params);
        if(abs(integrand) < 1e-30)return 0.0;
    }

    if(next_integral < ndim)
    {
        // Useful declarations
        double result=0.0, error, singularity_pts[16];
        int singularity_npts = 0;
        user_singularity_ptr singfunc;
        // Integrate using the selected algorithm
        switch(gsl_integrator_types[next_integral])
        {
            case NDIM_INT_QAG: // Relatively fast integration
                success=gsl_integration_qag(&gsl_functions[next_integral],
                        all_lower_limits[next_integral],
                        all_upper_limits[next_integral],
                        0.0, relative_tolerance, GSL_WORKSPACE_SIZE,
                        GSL_INTEG_GAUSS31, gsl_workspaces[next_integral],
                        &result, &error);
                break;
            case NDIM_INT_QAGS: // Singularity tolerant
                success=gsl_integration_qags(&gsl_functions[next_integral],
                        all_lower_limits[next_integral],
                        all_upper_limits[next_integral],
                        0.0, relative_tolerance, GSL_WORKSPACE_SIZE,
                        gsl_workspaces[next_integral], &result, &error);
                break;
            case NDIM_INT_QAGP: // Take care of known integrable singularities
                // Extract the sigularity function
                singfunc = user_specified_singularities[next_integral];
                // Use the function to get the singularities
                singularity_npts = singfunc(all_variables.data(),
                        user_specified_params, singularity_pts+1);
                // Set the lower and upper limits
                singularity_pts[0] = all_lower_limits[next_integral];
                singularity_pts[singularity_npts+1] = all_upper_limits[next_integral];
                // Finally call the integrator
                success=gsl_integration_qagp(&gsl_functions[next_integral], singularity_pts,
                        singularity_npts+2, 0.0, relative_tolerance, GSL_WORKSPACE_SIZE,
                        gsl_workspaces[next_integral], &result, &error);
                break;
            default:
                throw std::invalid_argument(integrator_name+"::evaluate The"
                        "integrator type is invalid.");
        }
        integrand *= result;

        if(success != GSL_SUCCESS)
            switch(verbosity)
            {
                case 1: // Do nothing
                    break;
                case 2: // Print variables but don't stop
                    std::cout << integrator_name << ": Exception " << gsl_strerror (success) <<
                        " occured in integral " << next_integral << " at coords: {";
                    for(auto x: all_variables) std::cout << x << ", ";
                    std::cout << "}" << std::endl << std::flush;
                    break;
                case 3: // Dump some points of the integrand and kill the program
                    std::cout << integrator_name << ": Exception " << gsl_strerror (success) <<
                        " occured in integral " << next_integral << " at coords: {";
                    for(auto x: all_variables) std::cout << x << ", ";
                    std::cout <<  "}. Dumping integrand in " << INTEGRAND_DUMP_FILE << std::endl;
                    dump_integrand(next_integral);
                    throw std::runtime_error(integrator_name+"::evaluate Integral did not converge.");
                    break;
                default:
                    break;

            }
    }
    return integrand;
}

// Function to evaluate the integral
double ndim_integrator::evaluate(void *input_params)
{
    // If we are verbose turn off the error handler
    if(verbosity)gsl_set_error_handler_off();
    user_specified_params = input_params; // Set the parameters
    return evaluate_integrand(0.0, -1); // This will do all the integrations
}

// Function to dump integral
void ndim_integrator::dump_integrand(int integral_number)
{
    // Step size for points
    double x, dx = (all_upper_limits[integral_number]-
            all_lower_limits[integral_number])/(INTEGRAND_DUMP_NPTS-1);
    std::fstream fout(INTEGRAND_DUMP_FILE, std::fstream::out); // Open file
    fout << "# x f(x)" << std::endl;
    for(x=all_lower_limits[integral_number]; x<=all_upper_limits[integral_number]; x+=dx)
        fout << x << " " << evaluate_integrand(x, integral_number) << std::endl;
    fout.close(); // Close the file
}
