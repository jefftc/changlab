/* cgeofnsmodule.c
 * 071028  created
 */

#include "Python.h"
#include <math.h>


/* Functions in this module. */


/* Module definition stuff */

static PyMethodDef CGeofnsMethods[] = {
  {NULL, NULL, 0, NULL}
};

static char cgeofns__doc__[] =
"This module provides optimized replacement functions.\n\
";

PyMODINIT_FUNC initcgeofns(void)
{
  (void) Py_InitModule3("cgeofns", CGeofnsMethods, cgeofns__doc__);
}
