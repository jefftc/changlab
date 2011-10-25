/* ciolibmodule.c
 * 071028  created
 */

#include "Python.h"
#include <math.h>

#if PYTHON_API_VERSION <= 1012  /* Python 2.3 is 1012.  */
#define Py_ssize_t int
#endif


/* Functions in this module. */


// Returns a new char * that must be freed.
char *charstar_FromPyObject(PyObject *py_obj)
{
    PyObject *py_str;
    char *c_str;
    Py_ssize_t length;
    char *str;

    // TODO: If py_obj is already a string, then don't do this.
    py_str = NULL;
    str = NULL;
    if(!(py_str = PyObject_Str(py_obj)))
	goto charstar_FromPyObject_cleanup;
    if(PyString_AsStringAndSize(py_str, &c_str, &length) == -1)
	goto charstar_FromPyObject_cleanup;
    if(!(str = (char *)malloc((length+1)*sizeof(*str)))) {
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
	goto charstar_FromPyObject_cleanup;
    }
    strcpy(str, c_str);

    return str;

 charstar_FromPyObject_cleanup:
    if(py_str)
	Py_DECREF(py_str);
    if(str)
	free(str);
    return NULL;
}


/* Convert Python object to a matrix of strings. */
// NOT TESTED YET.
int py2c_smatrix(PyObject *py_X, char **X, int *nrow, int *ncol)
{
    int i, j, index;
    int x;
    PyObject *py_X_row, *py_Fast_X_row, *py_item;

    *X = NULL;
    *nrow = *ncol = 0;
    py_X_row = py_Fast_X_row = py_item = NULL;

    if((*nrow = PySequence_Size(py_X)) == -1)
	goto py2c_matrix_err;
    if(*nrow < 1) {
	PyErr_SetString(PyExc_AssertionError, "empty matrix");
	goto py2c_matrix_err;
    }

    /* Get the number of columns in the matrix. */
    if(!(py_X_row = PySequence_GetItem(py_X, 0)))
	goto py2c_matrix_err;
    if((x = PySequence_Size(py_X_row)) == -1)
	goto py2c_matrix_err;
    Py_DECREF(py_X_row); py_X_row = NULL;
    *ncol = x;

    if(!(*X = (char *)malloc(*nrow**ncol*sizeof(**X)))) {
	PyErr_SetString(PyExc_MemoryError, "out of memory");
	goto py2c_matrix_err;
    }
    memset((void *)X, 0, *nrow**ncol*sizeof(**X));

    index = 0;
    for(i=0; i<*nrow; i++) {
	if(!(py_X_row = PySequence_GetItem(py_X, i)))
	    goto py2c_matrix_err;
	if(!(py_Fast_X_row = PySequence_Fast(py_X_row, "not sequence")))
	    goto py2c_matrix_err;
	Py_DECREF(py_X_row); py_X_row = NULL;
	x = PySequence_Fast_GET_SIZE(py_Fast_X_row);
	if(*ncol != x) {
	    PyErr_SetString(PyExc_AssertionError, "matrix length mismatch");
	    goto py2c_matrix_err;
	}
	for(j=0; j<*ncol; j++) {
	    /* borrowed ref */
	    py_item = PySequence_Fast_GET_ITEM(py_Fast_X_row, j);
	    X[index++] = (char *)charstar_FromPyObject(py_item);
	}
	py_item = NULL;
	Py_DECREF(py_Fast_X_row); py_Fast_X_row = NULL;
	if(PyErr_Occurred()) {
	    /* Optimization: This might be set in charstar_FromPyObject,
	       but check here instead of inside loop. */
	    goto py2c_matrix_err;
	}
    }

    return 1;

  py2c_matrix_err:
    if(py_item) { Py_DECREF(py_item); }
    if(py_X_row) { Py_DECREF(py_X_row); }
    if(py_Fast_X_row) { Py_DECREF(py_Fast_X_row); }
    if(*X) { 
	for(index=0; index<*nrow**ncol; index++) {
	    if(X[index] != NULL) {
		free(X[index]);
	    }
        }
	free(*X); *X = NULL; 
    }
    *nrow = *ncol = 0;
    return 0;
}



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


static char ciolib_cleanwrite__doc__[] = 
"XXX\n";

#define MAXLEN 1024*1024
static PyObject *ciolib_cleanwrite(
    PyObject *self, PyObject *args, PyObject *keywds)
{
    PyObject *py_data, *py_outhandle;
    char *delim;

    PyObject *py_data_row, *py_item, *py_str;
    char *c_str;
    Py_ssize_t length;
    int nrow, ncol;
    int r, c;
    int i;
    char buffer[MAXLEN];
    char *clean;

    static char *kwlist[] = {"data", "outhandle", "delim", NULL};

    py_data = py_outhandle = NULL;
    py_data_row = py_item = py_str = NULL;
    delim = "\t";
    if(!PyArg_ParseTupleAndKeywords(args, keywds, "OO|s", kwlist, 
				    &py_data, &py_outhandle, &delim))
	return NULL;

    if(!PySequence_Check(py_data)) {
	PyErr_SetString(PyExc_AssertionError, "not a matrix");
	goto cleanwrite_cleanup;
    }
    if(!PyFile_Check(py_outhandle)) {
	PyErr_SetString(PyExc_AssertionError, "no file");
	goto cleanwrite_cleanup;
    }

    if((nrow = PySequence_Size(py_data)) == -1)
	goto cleanwrite_cleanup;
    for(r=0; r<nrow; r++) {
	/* New reference. */
	if(!(py_data_row = PySequence_GetItem(py_data, r)))
	    goto cleanwrite_cleanup;
	if(!PySequence_Check(py_data_row)) {
	    PyErr_SetString(PyExc_AssertionError, "not a matrix");
	    goto cleanwrite_cleanup;
	}
	if((ncol = PySequence_Size(py_data_row)) == -1)
	    goto cleanwrite_cleanup;
	for(c=0; c<ncol; c++) {
	    // Make a copy of the string.
	    /* New reference. */
	    if(!(py_item = PySequence_GetItem(py_data_row, c)))
		goto cleanwrite_cleanup;
	    if(PyString_Check(py_item)) {
		if(PyString_AsStringAndSize(py_item, &c_str, &length) == -1)
		    goto cleanwrite_cleanup;
	    } else {
		if(!(py_str = PyObject_Str(py_item)))
		    goto cleanwrite_cleanup;
		if(PyString_AsStringAndSize(py_str, &c_str, &length) == -1)
		    goto cleanwrite_cleanup;
	    }
	    if(length >= MAXLEN) {
		PyErr_SetString(PyExc_AssertionError, "string too long");
		goto cleanwrite_cleanup;
	    }
	    clean = buffer;
	    strcpy(clean, c_str);

	    Py_DECREF(py_item);
	    py_item = NULL;
	    if(py_str) {
		Py_DECREF(py_str);
		py_str = NULL;
	    }

	    // Clean up the string.
	    // Strip the spaces at the beginning.
	    while(isspace(*clean)) {
		clean++;
	    }

	    // Strip the spaces from the end.
	    i = strlen(clean)-1;
	    while(i >= 0 && isspace(clean[i])) {
		clean[i--] = 0;
	    }
	    
	    // Convert spaces to " ".
	    i = 0;
	    while(clean[i]) {
		if(isspace(clean[i]))
		    clean[i] = ' ';
		i += 1;
	    }

	    // Write it to the file.
	    if(PyFile_WriteString(clean, py_outhandle) == -1)
		goto cleanwrite_cleanup;
	    if(c < ncol-1) {
		if(PyFile_WriteString(delim, py_outhandle) == -1)
		    goto cleanwrite_cleanup;
	    }
	}
	Py_DECREF(py_data_row);
	py_data_row = NULL;

	if(PyFile_WriteString("\n", py_outhandle) == -1)
	    goto cleanwrite_cleanup;
    }
    

 cleanwrite_cleanup:
    if(py_data_row)
	Py_DECREF(py_data_row);
    if(py_item)
	Py_DECREF(py_item);
    if(py_str)
	Py_DECREF(py_str);
    if(PyErr_Occurred())
	return NULL;
    Py_INCREF(Py_None);
    return Py_None;
}


/* Module definition stuff */

static PyMethodDef CIOlibMethods[] = {
  {"split_tdf", ciolib_split_tdf, METH_VARARGS, ciolib_split_tdf__doc__},
  {"strip_each", ciolib_strip_each, METH_VARARGS, ciolib_strip_each__doc__},
  {"cleanwrite", (PyCFunction)ciolib_cleanwrite, METH_VARARGS|METH_KEYWORDS, 
   ciolib_cleanwrite__doc__},
  {NULL, NULL, 0, NULL}
};

static char ciolib__doc__[] =
"This module provides optimized replacement functions.\n\
";

PyMODINIT_FUNC initciolib(void)
{
  (void) Py_InitModule3("ciolib", CIOlibMethods, ciolib__doc__);
}
