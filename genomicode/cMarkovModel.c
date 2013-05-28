/* cMarkovModel.c
 * 130520  created
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include "arrayobject.h"


#if PYTHON_API_VERSION <= 1012  /* Python 2.3 is 1012.  */
#define Py_ssize_t int
#endif



static char cMarkovModel__max_sum_and_index__doc__[] = 
"XXX\n";

static PyObject *cMarkovModel__max_sum_and_index(PyObject *self, 
  PyObject *args)
{
    PyArrayObject *py_object1, *py_object2;
    PyArrayObject *py_array1=NULL, *py_array2=NULL;
    double *array1, *array2;
    int len_array1, len_array2;
    int i, max_i;
    double sum, max_sum;
    PyObject *py_max_sum=NULL, *py_max_i=NULL;
    PyObject *py_retval=NULL;

    if(!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &py_object1, 
                         &PyArray_Type, &py_object2)) 
        return NULL;

    /* Make sure original object is double, to prevent conversions. */
    if(PyArray_TYPE(py_object1) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_AssertionError, "array1 not double");
        goto _max_sum_and_index_cleanup;
    }
    if(PyArray_TYPE(py_object2) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_AssertionError, "array2 not double");
        goto _max_sum_and_index_cleanup;
    }

    /* Make sure the data is contiguous. */
    if(!(py_array1 = (PyArrayObject *)
         PyArray_ContiguousFromAny((PyObject *)py_object1, NPY_DOUBLE, 1, 1)))
        goto _max_sum_and_index_cleanup;
    if(!(py_array2 = (PyArrayObject *)
         PyArray_ContiguousFromAny((PyObject *)py_object2, NPY_DOUBLE, 1, 1)))
         goto _max_sum_and_index_cleanup;

    /* Make sure both are one dimensional. */
    if(PyArray_NDIM(py_array1) != 1) {
        PyErr_SetString(PyExc_AssertionError, "array1 not 1 dimensional");
        goto _max_sum_and_index_cleanup;
    }
    if(PyArray_NDIM(py_array2) != 1) {
        PyErr_SetString(PyExc_AssertionError, "array2 not 1 dimensional");
        goto _max_sum_and_index_cleanup;
    }

    /* Make sure the lengths are the same. */
    len_array1 = (int)(PyArray_DIMS(py_array1)[0]);
    len_array2 = (int)(PyArray_DIMS(py_array2)[0]);
    if(len_array1 != len_array2) {
        PyErr_SetString(PyExc_AssertionError, "arrays not same length");
        goto _max_sum_and_index_cleanup;
    }
    if(len_array1 <= 0) {
        PyErr_SetString(PyExc_AssertionError, "empty arrays");
        goto _max_sum_and_index_cleanup;
    }


    /* Get the maximum sum and maximum index. */
    array1 = PyArray_GETPTR1(py_array1, 0);
    array2 = PyArray_GETPTR1(py_array2, 0);
    max_i = 0;
    max_sum = array1[0] + array2[0];
    for(i=1; i<len_array1; i++) {
        sum = array1[i] + array2[i];
        if(sum > max_sum) {
            max_sum = sum;
            max_i = i;
        }
    }

    /* Make Python objects for max_sum and max_i. */
    if(!(py_max_sum = PyFloat_FromDouble(max_sum)))
        goto _max_sum_and_index_cleanup;
    if(!(py_max_i = PyInt_FromLong(max_i)))
        goto _max_sum_and_index_cleanup;

    if(!(py_retval = PyTuple_New(2)))
        goto _max_sum_and_index_cleanup;
    PyTuple_SET_ITEM(py_retval, 0, py_max_sum); /* steals ref */
    PyTuple_SET_ITEM(py_retval, 1, py_max_i);
    py_max_sum = NULL;
    py_max_i = NULL;

    
  _max_sum_and_index_cleanup:
    if(py_array1) { Py_DECREF(py_array1); }
    if(py_array2) { Py_DECREF(py_array2); }
    if(py_max_sum) { Py_DECREF(py_max_sum); }
    if(py_max_i) { Py_DECREF(py_max_i); }
    if(PyErr_Occurred()) {
	if(py_retval) { Py_DECREF(py_retval); }
	return NULL;
    }
    return py_retval;
}


static char cMarkovModel__viterbi_h__doc__[] = 
"XXX\n";

static PyObject *cMarkovModel__viterbi_h(PyObject *self, 
  PyObject *args)
{
    PyArrayObject *py_scores, *py_lp_emission;
    PyObject *py_backtrace, *py_output;
    PyObject *py_lp_trans_j;
    int N, T;

    int i, j;
    int t;
    int k;
    PyObject *py_k;
    double *scores_prev=NULL;
    PyArrayObject *py_array=NULL, **py_arrays=NULL;  /* for py_lp_trans */
    PyObject *py_lp_trans_j_j;
    double *array;
    int max_i;
    PyObject *py_max_i=NULL;
    double sum, max_sum;
    double score;
    PyObject *py_temp=NULL;
    PyObject *py_backtrace_j;

    if(!PyArg_ParseTuple(args, "O!O!iiO!O!O!", 
                         &PyArray_Type, &py_scores, 
                         &PyList_Type, &py_backtrace, 
                         &N, &T, 
                         &PyList_Type, &py_output, 
                         &PyList_Type, &py_lp_trans_j,
                         &PyArray_Type, &py_lp_emission))
        return NULL;

    if(PyArray_NDIM(py_scores) != 2) {
        PyErr_SetString(PyExc_AssertionError, "scores not 2 dimensional");
        goto _viterbi_h_cleanup;
    }
    if(PyArray_TYPE(py_scores) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_AssertionError, "scores not double");
        goto _viterbi_h_cleanup;
    }
    if(PyList_GET_SIZE(py_backtrace) != N) {
        PyErr_SetString(PyExc_AssertionError, "size of backtrace");
        goto _viterbi_h_cleanup;
    }
    if(PyList_GET_SIZE(py_output) != T) {
        PyErr_SetString(PyExc_AssertionError, "size of output");
        goto _viterbi_h_cleanup;
    }
    if(PyList_GET_SIZE(py_lp_trans_j) != N) {
        PyErr_SetString(PyExc_AssertionError, "size of lp_trans_j");
        goto _viterbi_h_cleanup;
    }
    if(PyArray_NDIM(py_lp_emission) != 2) {
        PyErr_SetString(PyExc_AssertionError, "lp_emission not 2 dimensional");
        goto _viterbi_h_cleanup;
    }
    if(PyArray_TYPE(py_lp_emission) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_AssertionError, "lp_emission not double");
        goto _viterbi_h_cleanup;
    }

    for(i=0; i<N; i++) {
        py_backtrace_j = PyList_GetItem(py_backtrace, i); /* borrowed */
        if(!PyList_Check(py_backtrace_j)) {
            PyErr_SetString(PyExc_AssertionError, "not list");
            goto _viterbi_h_cleanup;
        }
        if(PyList_Size(py_backtrace_j) != T) {
            PyErr_SetString(PyExc_AssertionError, "list size");
            goto _viterbi_h_cleanup;
        }
    }

    if(!(scores_prev = malloc(N*sizeof(*scores_prev)))) {
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        goto _viterbi_h_cleanup;
    }

    if(!(py_arrays = malloc(N*sizeof(PyArrayObject *)))) {
        PyErr_SetString(PyExc_MemoryError, "out of memory");
        goto _viterbi_h_cleanup;
    }
    for(i=0; i<N; i++)
        py_arrays[i] = NULL;
    for(i=0; i<N; i++) {
        /* borrowed ref */
        if(!(py_lp_trans_j_j = PyList_GetItem(py_lp_trans_j, i)))
            goto _viterbi_h_cleanup;
        if(!PyArray_Check(py_lp_trans_j_j)) {
            PyErr_SetString(PyExc_AssertionError, "lp_trans_j[j] not array");
            goto _viterbi_h_cleanup;
        }
        if(PyArray_TYPE((PyArrayObject *)py_lp_trans_j_j) != NPY_DOUBLE) {
            PyErr_SetString(PyExc_AssertionError, "lp_trans_j[j] not double");
            goto _viterbi_h_cleanup;
        }
        if(!(py_array = (PyArrayObject *)
             PyArray_ContiguousFromAny(py_lp_trans_j_j, NPY_DOUBLE, 1, 1)))
            goto _viterbi_h_cleanup;
        if(PyArray_NDIM(py_array) != 1) {
            PyErr_SetString(PyExc_AssertionError, "not 1 dimensional");
            goto _viterbi_h_cleanup;
        }
        if((int)PyArray_DIMS(py_array)[0] != N) {
            PyErr_SetString(PyExc_AssertionError, "array length");
            goto _viterbi_h_cleanup;
        }
        py_arrays[i] = py_array;
        py_array = NULL;
    }


    for(t=1; t<T; t++) {
        py_k = PyList_GetItem(py_output, t); /* borrowed ref */
        k = PyInt_AsLong(py_k);
        if(PyErr_Occurred())
            goto _viterbi_h_cleanup;

        for(i=0; i<N; i++)
            scores_prev[i] = *(double *)PyArray_GETPTR2(py_scores, i, t-1);

        for(j=0; j<N; j++) {
            /* Get the maximum sum and maximum index. */
            array = PyArray_GETPTR1(py_arrays[j], 0);

            max_i = 0;
            max_sum = scores_prev[0] + array[0];
            for(i=1; i<N; i++) {
                sum = scores_prev[i] + array[i];
                if(sum > max_sum) {
                    max_sum = sum;
                    max_i = i;
                }
            }

            score = max_sum + *(double *)PyArray_GETPTR2(py_lp_emission, j, k);
            *(double *)PyArray_GETPTR2(py_scores, j, t) = score;

            if(!(py_max_i = PyInt_FromLong(max_i)))
                goto _viterbi_h_cleanup;
            if(!(py_temp = PyList_New(1)))
                goto _viterbi_h_cleanup;
            PyList_SET_ITEM(py_temp, 0, py_max_i); /* steals ref */
            py_max_i = NULL;

            py_backtrace_j = PyList_GET_ITEM(py_backtrace, j); /* borrowed */
            PyList_SetItem(py_backtrace_j, t, py_temp); /* steals ref */
            py_temp = NULL;
        }
    }

  _viterbi_h_cleanup:
    if(scores_prev) { free(scores_prev); }
    if(py_arrays) { 
        for(i=0; i<N; i++)
            if(py_arrays[i] != NULL)
                Py_DECREF(py_arrays[i]);
        free(py_arrays);
    }
    if(py_array) { Py_DECREF(py_array); }
    if(py_max_i) { Py_DECREF(py_max_i); }
    if(py_temp) { Py_DECREF(py_temp); }
    if(PyErr_Occurred()) { return NULL; }

    Py_INCREF(Py_None);
    return Py_None;
}


/* Module definition stuff */

static PyMethodDef cMarkovModelMethods[] = {
  {"_max_sum_and_index", cMarkovModel__max_sum_and_index, METH_VARARGS, 
   cMarkovModel__max_sum_and_index__doc__},
  {"_viterbi_h", cMarkovModel__viterbi_h, METH_VARARGS, 
   cMarkovModel__viterbi_h__doc__},
  {NULL, NULL, 0, NULL}
};

static char cMarkovModel__doc__[] =
"This module provides optimized replacement functions.\n\
";

PyMODINIT_FUNC initcMarkovModel(void)
{
  (void) Py_InitModule3("cMarkovModel", cMarkovModelMethods, 
                        cMarkovModel__doc__);
  import_array();
}
