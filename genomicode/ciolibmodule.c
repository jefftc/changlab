/* ciolibmodule.c
 * 071028  created
 */

#include "Python.h"
#include <math.h>

#if PYTHON_API_VERSION <= 1012  /* Python 2.3 is 1012.  */
#define Py_ssize_t int
#endif


/* Functions in this module. */

static char ciolib_split_tdf__doc__[] = 
"XXX\n";

static PyObject *ciolib_split_tdf(
  PyObject *self, PyObject *args)
{
    char *string;
    int start, end;
    PyObject *py_matrix, *py_row, *py_item;

    py_matrix = py_row = py_item = NULL;
    if(!PyArg_ParseTuple(args, "s", &string))
	return NULL;
    if(!(py_matrix = PyList_New(0)))
	return NULL;

    start = end = 0;
    while(string[end]) {
	start = end;
	/* Advance end to the next word. */
	while(string[end] && string[end] != '\n' && string[end] != '\r' && 
	      string[end] != '\t')
	    end++;

	/* Create the new word. */
	if((!(py_item=PyString_FromStringAndSize(&string[start], end-start))))
	    goto split_tdf_cleanup;

	/* Insert the new word into a row, creating it if necessary. */
	if(!py_row && !(py_row = PyList_New(0)))
	    goto split_tdf_cleanup;
	if(PyList_Append(py_row, py_item) == -1)
	    goto split_tdf_cleanup;
	Py_DECREF(py_item); py_item = NULL;
	
	/* Check if I'm at the end of the row or file.  If so, add the
	   row to the matrix */
	if(!string[end] || string[end] == '\n' || string[end] == '\r') {
	    if(PyList_Append(py_matrix, py_row) == -1)
		goto split_tdf_cleanup;
	    Py_DECREF(py_row); py_row = NULL;
	}

	/* Advance to the next word. */
	if(string[end] == '\t')
	    end++;
	else if(string[end] == '\n' || string[end] == '\r') {
	    /* BUG: Will skip blank lines. */
	    while(string[end] == '\n' || string[end] == '\r')
		end++;
	} else if(string[end]) {
	    PyErr_SetString(PyExc_AssertionError, "Invalid string");
	    goto split_tdf_cleanup;
	}
    }

    /* Add the last row to the matrix. */
    if(py_row) {
	if(PyList_Append(py_matrix, py_row) == -1)
	    goto split_tdf_cleanup;
	Py_DECREF(py_row); py_row = NULL;
    }

 split_tdf_cleanup:
    if(py_item) {
	Py_DECREF(py_item);
    }
    if(py_row) {
	Py_DECREF(py_row);
    }
    if(PyErr_Occurred()) {
	if(py_matrix) {
	    Py_DECREF(py_matrix);
	}
	return NULL;
    }
    return py_matrix;
}

static char ciolib_strip_each__doc__[] = 
"XXX\n";

static PyObject *ciolib_strip_each(PyObject *self, PyObject *args)
{
    PyObject *py_L, *py_retval;
    Py_ssize_t len_L;

    Py_ssize_t i;
    PyObject *py_item, *py_sitem;
    char *item, *sitem;
    Py_ssize_t len_item, len_sitem;

    py_item = py_sitem = NULL;
    py_retval = NULL;

    if(!PyArg_ParseTuple(args, "O", &py_L))
	return NULL;
    if(!PySequence_Check(py_L)) {
	PyErr_SetString(PyExc_AssertionError, "arg must be a sequence");
	goto strip_each_cleanup;
    }
    if((len_L = PySequence_Size(py_L)) == -1) {
	PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
	goto strip_each_cleanup;
    }
    if(!(py_retval = PyList_New(len_L)))
	goto strip_each_cleanup;

    for(i=0; i<len_L; i++) {
	if(!(py_item = PySequence_GetItem(py_L, i)))  /* New reference. */
	    goto strip_each_cleanup;
	if((PyString_AsStringAndSize(py_item, &item, &len_item)) == -1)
	    /* May raise TypeError if py_item is not a string. */
	    goto strip_each_cleanup;

	sitem = item;
	len_sitem = len_item;

	/* Skip the spaces at the beginning. */
	while(len_sitem && isspace(sitem[0])) {
	    sitem++;
	    len_sitem--;
	}
	/* Skip the spaces at the end. */
	while(len_sitem && isspace(sitem[len_sitem-1])) {
	    len_sitem--;
	}

	if(sitem == item && len_sitem == len_item) {
	    /* Nothing to strip, reuse the same string. */
	    py_sitem = py_item;
	    Py_INCREF(py_sitem);
	} else {
	    /* Make a new string. */
	    if(!(py_sitem = PyString_FromStringAndSize(sitem, len_sitem)))
		goto strip_each_cleanup;
	}
       
	Py_DECREF(py_item);
	py_item = NULL;
	PyList_SET_ITEM(py_retval, i, py_sitem);  /* Steals the reference. */
	py_sitem = NULL;
    }

    
 strip_each_cleanup:
    if(py_item) {
	Py_DECREF(py_item);
    }
    if(py_sitem) {
	Py_DECREF(py_sitem);
    }
    if(PyErr_Occurred()) {
	if(py_retval) {
	    Py_DECREF(py_retval);
	}
	return NULL;
    }
    return py_retval;
}


/* Module definition stuff */

static PyMethodDef CIOlibMethods[] = {
  {"split_tdf", ciolib_split_tdf, METH_VARARGS, ciolib_split_tdf__doc__},
  {"strip_each", ciolib_strip_each, METH_VARARGS, ciolib_strip_each__doc__},
  {NULL, NULL, 0, NULL}
};

static char ciolib__doc__[] =
"This module provides optimized replacement functions.\n\
";

PyMODINIT_FUNC initciolib(void)
{
  (void) Py_InitModule3("ciolib", CIOlibMethods, ciolib__doc__);
}
