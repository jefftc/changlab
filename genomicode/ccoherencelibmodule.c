/* ccoherencelibmodule.c
 * 071024  created
 */

#include "Python.h"
#include "mathlib.h"

#if PYTHON_API_VERSION <= 1012  /* Python 2.3 is 1012.  */
#define Py_ssize_t int
#endif

double _convert_z_to_R(double Z_CUTOFF, int N, double DELTA)
{
    int i;
    double R;
    double min_R, max_R;
    double z;

    i = 0;
    min_R = 0.0; max_R = 1.0;
    while(min_R < max_R) {
        R = (min_R + max_R) / 2;
        z = fisher_z(R, N);
        if(fabs(z-Z_CUTOFF) < DELTA)
            break;
	else if(z > Z_CUTOFF)
            max_R = R - 1E-5;
        else if(z < Z_CUTOFF)
            min_R = R + 1E-5;
        i++;
    }
    return R;
}

//int _find_medoid(double *CORS, int nrow, int ncol, double CUTOFF)
int _find_medoid(float *CORS, int nrow, int ncol, float CUTOFF)
{
    int i, j;
    int *pcors, *ncors;
    int num_cors;
    int max_i, max_cors;
    float *x;

    pcors = ncors = NULL;
    max_cors = max_i = -1;

    if(!(pcors = (int *)malloc(nrow*sizeof(*pcors)))) {
	PyErr_SetString(PyExc_MemoryError, "out of memory");
	goto _find_medoid_cleanup;
    }
    if(!(ncors = (int *)malloc(nrow*sizeof(*ncors)))) {
	PyErr_SetString(PyExc_MemoryError, "out of memory");
	goto _find_medoid_cleanup;
    }
    memset(pcors, 0, nrow*sizeof(*pcors));
    memset(ncors, 0, nrow*sizeof(*ncors));
    
    for(i=1; i<nrow; i++) {
	x = CORS + i*ncol;
	for(j=i+1; j<nrow; j++) {
	    if(x[j] >= CUTOFF) {
		pcors[i] += 1;
		pcors[j] += 1;
	    } else if(x[j] <= -CUTOFF) {
		ncors[i] += 1;
		ncors[j] += 1;
	    }
	}
    }

    max_cors = max_i = -1;
    for(i=0; i<nrow; i++) {
	num_cors = (pcors[i] > ncors[i]) ? pcors[i] : ncors[i];
	if(max_cors == -1 || num_cors > max_cors) {
	    max_cors = num_cors;
	    max_i = i;
	}
    }

 _find_medoid_cleanup:
    if(pcors) { free(pcors); }
    if(ncors) { free(ncors); }
    return max_i;
}

/* int _find_medoid(float *CORS, int nrow, int ncol, float CUTOFF)
{
    int i, j;
    int pcors, ncors, num_cors;
    int max_i, max_cors;
    float *x;

    max_cors = max_i = -1;
    for(i=0; i<nrow; i++) {
	x = CORS + i*ncol;
	ncors = pcors = 0;
	for(j=0; j<ncol; j++) {
	    if(x[j] >= CUTOFF)
		pcors += 1;
	    else if(x[j] <= -CUTOFF)
		ncors += 1;
	    //if(x[j] >= CUTOFF || x[j] <= -CUTOFF)
	    //num_cors++;
	}
	num_cors = (pcors > ncors) ? pcors : ncors;
	//num_cors = pcors + ncors;
	if(max_cors == -1 || num_cors > max_cors) {
	    max_cors = num_cors;
	    max_i = i;
	}
    }
    return max_i;
    } */

/* Functions in this module. */

static char ccoherencelib__find_medoid__doc__[] = 
"XXX\n";

static PyObject *ccoherencelib__find_medoid(PyObject *self, PyObject *args)
{
    PyObject *py_CORS;
    float CUTOFF;
    float *CORS;
    int nrow, ncol;
    int medoid_i;

    if(!PyArg_ParseTuple(args, "Of", &py_CORS, &CUTOFF))
	return NULL;
    if(!py2c_fmatrix(py_CORS, &CORS, &nrow, &ncol))
	return NULL;
    medoid_i = _find_medoid(CORS, nrow, ncol, CUTOFF);
    free(CORS);
    return PyInt_FromLong(medoid_i);
}

static char ccoherencelib__calc_coherence_score_h__doc__[] = 
"XXX\n";

static PyObject *ccoherencelib__calc_coherence_score_h(
  PyObject *self, PyObject *args)
{
    PyObject *py_X;
    float Z_CUTOFF;

    float *X;
    int nrow, ncol;

    float *CORS;
    float DELTA = 0.001;
    float R_CUTOFF;

    int MEDOID_I;
    float *cors_med;

    int i;
    int *PCORS_I, *NCORS_I;
    int num_pcors, num_ncors;
    int *TMP_I;
    int num_tmp;
    float R;

    PyObject *py_retval;
    PyObject *py_NROW, *py_NCOL;
    PyObject *py_CORS, *py_ZCORS;
    PyObject *py_MEDOID_I;
    PyObject *py_PCORS_I, *py_NCORS_I;

    X = NULL;
    CORS = NULL;
    PCORS_I = NCORS_I = NULL;

    py_retval = NULL;
    py_NROW = py_NCOL = NULL;
    py_CORS = py_ZCORS = py_MEDOID_I = NULL;
    py_PCORS_I = py_NCORS_I = NULL;

    if(!PyArg_ParseTuple(args, "Of", &py_X, &Z_CUTOFF))
	return NULL;
    if(!py2c_fmatrix(py_X, &X, &nrow, &ncol))
	return NULL;
    if(ncol < 2) {
	PyErr_SetString(PyExc_AssertionError, "at least 2 samples required");
	return NULL;
    }

    if(!(CORS = cor_byrow_f(X, nrow, ncol, 1))) {
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
	goto _calc_coherence_score_cleanup;
    }

    /* Find the R-value that corresponds to Z_CUTOFF. */
    R_CUTOFF = _convert_z_to_R(Z_CUTOFF, ncol, DELTA);

    /* Find the medoid gene (correlated with the most neighbors). */
    MEDOID_I = _find_medoid(CORS, nrow, nrow, R_CUTOFF);
    if(MEDOID_I < 0)
	goto _calc_coherence_score_cleanup;

    /* Find genes correlated with the centroid. */
    if(!(PCORS_I = (int *)malloc(nrow*sizeof(*PCORS_I)))) {
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
	goto _calc_coherence_score_cleanup;
    }
    if(!(NCORS_I = (int *)malloc(nrow*sizeof(*NCORS_I)))) {
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
	goto _calc_coherence_score_cleanup; 
    }
    num_pcors = num_ncors = 0;
    cors_med = CORS + MEDOID_I*nrow;
    for(i=0; i<nrow; i++) {
	R = *cors_med++;
	if(R >= R_CUTOFF)
	    PCORS_I[num_pcors++] = i;
	else if(R <= -R_CUTOFF)
	    NCORS_I[num_ncors++] = i;
    }

    /* Arbitrarily set positive correlation as the group with the most
       genes. */
    if(num_pcors < num_ncors) {
	TMP_I = PCORS_I; PCORS_I = NCORS_I; NCORS_I = TMP_I;
	num_tmp = num_pcors; num_pcors = num_ncors; num_ncors = num_tmp;
	} 

    /* Format the return values. */
    if(!(py_NROW = PyInt_FromLong((long)nrow)))
	goto _calc_coherence_score_cleanup;
    if(!(py_NCOL = PyInt_FromLong((long)ncol)))
	goto _calc_coherence_score_cleanup;
    if(!(py_CORS = c2py_fmatrix(CORS, nrow, nrow)))
	goto _calc_coherence_score_cleanup;
    Py_INCREF(Py_None);
    py_ZCORS = Py_None;
    if(!(py_MEDOID_I = PyInt_FromLong(MEDOID_I)))
	goto _calc_coherence_score_cleanup;
    if(!(py_PCORS_I = c2py_ivector(PCORS_I, num_pcors)))
	goto _calc_coherence_score_cleanup;
    if(!(py_NCORS_I = c2py_ivector(NCORS_I, num_ncors)))
	goto _calc_coherence_score_cleanup;
    if(!(py_retval = PyList_New(7)))
	goto _calc_coherence_score_cleanup;
    PyList_SET_ITEM(py_retval, 0, py_NROW);  /* steals references */
    PyList_SET_ITEM(py_retval, 1, py_NCOL);  /* steals references */
    PyList_SET_ITEM(py_retval, 2, py_CORS);  /* steals references */
    PyList_SET_ITEM(py_retval, 3, py_ZCORS);
    PyList_SET_ITEM(py_retval, 4, py_MEDOID_I);
    PyList_SET_ITEM(py_retval, 5, py_PCORS_I);
    PyList_SET_ITEM(py_retval, 6, py_NCORS_I);
    py_NROW = py_NCOL = NULL;
    py_CORS = py_ZCORS = py_MEDOID_I = NULL;
    py_PCORS_I = py_NCORS_I = NULL;

 _calc_coherence_score_cleanup:
    if(X) { free(X); }
    if(CORS) { free(CORS); }
    if(PCORS_I) { free(PCORS_I); }
    if(NCORS_I) { free(NCORS_I); }
    if(py_NROW) { Py_DECREF(py_NROW); }
    if(py_NCOL) { Py_DECREF(py_NCOL); }
    if(py_CORS) { Py_DECREF(py_CORS); }
    if(py_ZCORS) { Py_DECREF(py_ZCORS); }
    if(py_MEDOID_I) { Py_DECREF(py_MEDOID_I); }
    if(py_PCORS_I) { Py_DECREF(py_PCORS_I); }
    if(py_NCORS_I) { Py_DECREF(py_NCORS_I); }
    return py_retval;
}

static char ccoherencelib__calc_consensus_matrix__doc__[] = 
"XXX\n";

static PyObject *ccoherencelib__calc_consensus_matrix(
  PyObject *self, PyObject *args)
{
    PyObject *py_X;

    float *X, *X_i, *X_j;
    int nrow, ncol;

    float *C;
    int i, j, k;
    int count;

    PyObject *py_retval;

    X = NULL;
    C = NULL;
    py_retval = NULL;

    if(!PyArg_ParseTuple(args, "O", &py_X))
	return NULL;
    if(!py2c_fmatrix(py_X, &X, &nrow, &ncol))
	return NULL;
    if(nrow < 1 || ncol < 1) {
	PyErr_SetString(PyExc_AssertionError, "empty matrix");
	return NULL;
    }

    if(!(C = (float *)malloc(nrow*nrow*sizeof(*C))))
	goto _calc_consensus_matrix_cleanup;
    memset(C, 0, nrow*nrow*sizeof(*C));

    for(i=0; i<nrow; i++) {
	for(j=i; j<nrow; j++) {
	    X_i = X + i*ncol;
	    X_j = X + j*ncol;
	    count = 0;
	    for(k=0; k<ncol; k++) {
		count += *X_i++ == *X_j++;
                /* count += X[i*ncol+k] == X[j*ncol+k]; */
	    }
	    C[i*nrow+j] = (float)count / ncol;
	}
    }
    if(!(py_retval = c2py_fmatrix(C, nrow, nrow)))
	goto _calc_consensus_matrix_cleanup;

 _calc_consensus_matrix_cleanup:
    if(X) { free(X); }
    if(C) { free(C); }
    return py_retval;
}

static char ccoherencelib__get_unique_indexes__doc__[] = 
"XXX\n";

static PyObject *ccoherencelib__get_unique_indexes(
  PyObject *self, PyObject *args)
{
    PyObject *py_names, *py_good_names;
    Py_ssize_t len_names;
    PyObject *py_names_fast;

    PyObject *py_seen;
    int *indexes;
    int num_indexes;

    Py_ssize_t i;
    PyObject *py_name;
    //int contains;

    PyObject *py_retval;

    py_names_fast = NULL;
    py_seen = NULL;
    indexes = NULL;
    py_retval = NULL;

    if(!PyArg_ParseTuple(args, "OO", &py_names, &py_good_names))
	return NULL;

    /* Set up the names variable. */
    if(!PySequence_Check(py_names)) {
	PyErr_SetString(PyExc_AssertionError, "arg must be a sequence");
	goto _get_unique_indexes_cleanup;
    }
    if((len_names = PySequence_Size(py_names)) == -1) {
	PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
	goto _get_unique_indexes_cleanup;
    }
    if(!(py_names_fast = PySequence_Fast(py_names, "bad fast sequence")))
	goto _get_unique_indexes_cleanup;

    /* Set up the good_names variable. */
    if(!PyDict_Check(py_good_names)) {
	PyErr_SetString(PyExc_AssertionError, "arg must be a dictionary");
	goto _get_unique_indexes_cleanup;
    }

    if(!(py_seen = PyDict_New()))
	goto _get_unique_indexes_cleanup;
    num_indexes = 0;
    if(!(indexes = (int *)malloc(len_names*sizeof(*indexes)))) {
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
	goto _get_unique_indexes_cleanup;
    }

    for(i=0; i<len_names; i++) {
	/* borrowed ref */
	if(!(py_name = PySequence_Fast_GET_ITEM(py_names_fast, i)))
	    goto _get_unique_indexes_cleanup;
	if(!(PyString_Check(py_name))) {
	    PyErr_SetString(PyExc_AssertionError, "item must be a string");
	    goto _get_unique_indexes_cleanup;
	}
	
	/* Make sure name in good_names. */
	if(!PyDict_GetItem(py_good_names, py_name))
	    continue;
	//contains = PyDict_Contains(py_good_names, py_name);
	//if(contains == -1)
	//    goto _get_unique_indexes_cleanup;
	//if(!contains)
	//    continue;

	/* Make sure not in seen. */
	//contains = PyDict_Contains(py_seen, py_name);
	//if(contains == -1)
	//    goto _get_unique_indexes_cleanup;
	//if(contains)
	//    continue;
	if(PyDict_GetItem(py_seen, py_name))
	    continue;
	if(PyDict_SetItem(py_seen, py_name, Py_None) == -1)
	    goto _get_unique_indexes_cleanup;

	if(num_indexes >= len_names) {
	    PyErr_SetString(PyExc_AssertionError, "no more indexes");
	    goto _get_unique_indexes_cleanup;
	}
	indexes[num_indexes++] = (int)i;
    }

    if(!(py_retval = c2py_ivector(indexes, num_indexes)))
	goto _get_unique_indexes_cleanup;

 _get_unique_indexes_cleanup:
    if(py_names_fast) { Py_DECREF(py_names_fast); }
    if(py_seen) { Py_DECREF(py_seen); }
    if(indexes) { free(indexes); }
    return py_retval;
}

/* Module definition stuff */

static PyMethodDef cCoherenceLibMethods[] = {
  {"_find_medoid", ccoherencelib__find_medoid, METH_VARARGS, 
   ccoherencelib__find_medoid__doc__},
  {"_calc_coherence_score_h", ccoherencelib__calc_coherence_score_h, 
   METH_VARARGS, ccoherencelib__calc_coherence_score_h__doc__},
  /* Out of date. */
  /*  {"_calc_consensus_matrix", ccoherencelib__calc_consensus_matrix, 
      METH_VARARGS, ccoherencelib__calc_consensus_matrix__doc__}, */
  {"_get_unique_indexes", ccoherencelib__get_unique_indexes, 
   METH_VARARGS, ccoherencelib__get_unique_indexes__doc__},
  {NULL, NULL, 0, NULL}
};

static char ccoherencelib__doc__[] =
"This module provides optimized replacement functions.\n\
";

PyMODINIT_FUNC initccoherencelib(void)
{
  (void) Py_InitModule3("ccoherencelib", cCoherenceLibMethods, 
      ccoherencelib__doc__);
}
