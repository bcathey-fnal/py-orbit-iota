#ifndef WRAP_MCMILLAN_H
#define WRAP_MCMILLAN_H

#include "Python.h"

namespace wrap_mcmillan{

#ifdef __cplusplus
extern "C" {
#endif
    // Register module so that python knows about it.
    void initmcmillan(void);
    // Register the McMillan class
    void initclassMcMillan(PyObject* module);

#ifdef __cplusplus
}
#endif  // __cplusplus

} // end of namespace wrap_mcmillan

#endif // WRAP_MCMILLAN_H
