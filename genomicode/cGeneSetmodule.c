/* cGeneSetmodule.c
 * 080714  created
 */

#include "Python.h"

#if PYTHON_API_VERSION <= 1012  /* Python 2.3 is 1012.  */
#define Py_ssize_t int
#endif

int is_valid_gene_id(const char *id)
{
    int i;

    i = 0;

    /* Skip to the first non-space character. */
    while(id[i] && isspace(id[i]))
	i++;

    if(!id[i])                         /* "" */
	return 0;
    if(id[i] == '0' && !id[i+1])       /* "0"; BUG: will accept "0 " too. */
	return 0;
    if(id[i] == '-')                   /* Starts with "-" */
	return 0;
    if(strstr(id+i, "///"))            /* Contains "///" */
	return 0;
    return 1;
}

static char cGeneSet_is_valid_gene_id__doc__[] = 
"XXX\n";

static PyObject *cGeneSet_is_valid_gene_id(PyObject *self, PyObject *args)
{
    PyObject *py_id;
    char *id;
    Py_ssize_t len_id;

    if(!PyArg_ParseTuple(args, "O", &py_id))
	return NULL;
    /* None is not a valid gene id. */
    if(py_id == Py_None) 
	return PyInt_FromLong(0);

    if(!PyString_Check(py_id)) {
	PyErr_SetString(PyExc_TypeError, "id should be a string");
 	return NULL;
    }
    if(PyString_AsStringAndSize(py_id, &id, &len_id) == -1)
	return NULL;
    return PyInt_FromLong(is_valid_gene_id(id));
}



static char cGeneSet_select_valid_gene_id__doc__[] = 
"XXX\n";

static PyObject *cGeneSet_select_valid_gene_id(
    PyObject *self, PyObject *args)
{
    PyObject *py_L;
    PyObject *py_L_fast;
    Py_ssize_t len_L;

    int *I;
    int len_I;
    PyObject *py_L_new;

    Py_ssize_t i;
    PyObject *py_item;
    char *item;
    Py_ssize_t len_item;

    py_L_fast = NULL;
    py_item = NULL;
    I = NULL;
    py_L_new = NULL;

    if(!PyArg_ParseTuple(args, "O", &py_L))
	return NULL;
    if(!PySequence_Check(py_L)) {
	PyErr_SetString(PyExc_AssertionError, "arg must be a sequence");
	goto select_valid_gene_id_cleanup;
    }
    if(!(py_L_fast = PySequence_Fast(py_L, "unable to get fast sequence")))
	goto select_valid_gene_id_cleanup;
    if((len_L = PySequence_Fast_GET_SIZE(py_L_fast)) == -1) {
	PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
	goto select_valid_gene_id_cleanup;
    }

    if(!(I = (int *)malloc(len_L*sizeof(*I))))
	goto select_valid_gene_id_cleanup;
    len_I = 0;

    for(i=0; i<len_L; i++) {
	if(!(py_item = PySequence_Fast_GET_ITEM(py_L_fast, i)))  /* borrowed */
	    goto select_valid_gene_id_cleanup;
	if(py_item == Py_None) {
	    /* None is not valid gene id, so skip. */
	    py_item = NULL;
	    continue;
	}

	if((PyString_AsStringAndSize(py_item, &item, &len_item)) == -1) {
	    /* May raise TypeError if py_item is not a string. */
	    py_item = NULL;
	    goto select_valid_gene_id_cleanup;
	}
	if(is_valid_gene_id(item))
	    I[len_I++] = i;
	py_item = NULL;
    }

    if(!(py_L_new = PyList_New(len_I)))
	goto select_valid_gene_id_cleanup;
    for(i=0; i<len_I; i++) {
	/* borrowed */
	if(!(py_item = PySequence_Fast_GET_ITEM(py_L_fast, I[i])))
	    goto select_valid_gene_id_cleanup;
	Py_INCREF(py_item);
	PyList_SET_ITEM(py_L_new, i, py_item);  /* Steals reference. */
	py_item = NULL;
    }

 select_valid_gene_id_cleanup:
    if(py_L_fast) {
	Py_DECREF(py_L_fast);
    }
    if(py_item) {
	Py_DECREF(py_item);
    }
    if(I) {
	free(I);
    }
    if(PyErr_Occurred()) {
	if(py_L_new) {
	    Py_DECREF(py_L_new);
	}
	return NULL;
    }
    return py_L_new;
}

static char cGeneSet_select_valid_gene_id_I__doc__[] = 
"XXX\n";

static PyObject *cGeneSet_select_valid_gene_id_I(
    PyObject *self, PyObject *args)
{
    PyObject *py_L;
    Py_ssize_t len_L;

    int *I;
    int len_I;
    PyObject *py_I;

    Py_ssize_t i;
    PyObject *py_item;
    char *item;
    Py_ssize_t len_item;

    py_item = NULL;
    I = NULL;
    py_I = NULL;

    if(!PyArg_ParseTuple(args, "O", &py_L))
	return NULL;
    if(!PySequence_Check(py_L)) {
	PyErr_SetString(PyExc_AssertionError, "arg must be a sequence");
	goto select_valid_gene_id_I_cleanup;
    }
    if((len_L = PySequence_Size(py_L)) == -1) {
	PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
	goto select_valid_gene_id_I_cleanup;
    }

    if(!(I = (int *)malloc(len_L*sizeof(*I))))
	goto select_valid_gene_id_I_cleanup;
    len_I = 0;

    for(i=0; i<len_L; i++) {
	if(!(py_item = PySequence_GetItem(py_L, i)))  /* New reference. */
	    goto select_valid_gene_id_I_cleanup;
	if(py_item == Py_None) {
	    /* None is not valid gene id, so skip. */
	    Py_DECREF(py_item);
	    py_item = NULL;
	    continue;
	}

	if((PyString_AsStringAndSize(py_item, &item, &len_item)) == -1)
	    /* May raise TypeError if py_item is not a string. */
	    goto select_valid_gene_id_I_cleanup;
	if(is_valid_gene_id(item))
	    I[len_I++] = i;
	
	Py_DECREF(py_item);
	py_item = NULL;
    }

    if(!(py_I = PyList_New(len_I)))
	goto select_valid_gene_id_I_cleanup;
    for(i=0; i<len_I; i++) {
	if(!(py_item = PyInt_FromLong(I[i])))
	    goto select_valid_gene_id_I_cleanup;
	PyList_SET_ITEM(py_I, i, py_item);  /* Steals reference. */
	py_item = NULL;
    }

 select_valid_gene_id_I_cleanup:
    if(py_item) {
	Py_DECREF(py_item);
    }
    if(I) {
	free(I);
    }
    if(PyErr_Occurred()) {
	if(py_I) {
	    Py_DECREF(py_I);
	}
	return NULL;
    }
    return py_I;
}


int seeded = 0;

int set_random_seed() 
{
    int i;
    FILE *handle;
    char random_bytes[4];
    int seed;

    if(!(handle = fopen("/dev/urandom", "r")))
	return 0;
    fread((void *)random_bytes, 4, 1, handle);
    seed = 0;
    for(i=0; i<4; i++) {
	seed = seed << 8;
	seed += (int)random_bytes[i];
    }
    //printf("Seed: %u\n", seed);
    srand(seed);
    return 1;
}

static char cGeneSet_sample__doc__[] = 
"XXX\n";

static PyObject *cGeneSet_sample(PyObject *self, PyObject *args)
{
    PyObject *py_L;
    Py_ssize_t len_L;
    int n;

    int i, j;
    unsigned char *chosen;
    int num_chosen;

    PyObject *py_L_fast;
    PyObject *py_item, *py_retval;

    chosen = NULL;
    py_L_fast = NULL;
    py_item = NULL;
    py_retval = NULL;

    if(!PyArg_ParseTuple(args, "Oi", &py_L, &n))
	return NULL;
    if((len_L = PySequence_Size(py_L)) == -1)
	return NULL;
    if(len_L < n) {
	PyErr_SetString(PyExc_AssertionError, "too few objects");
	return NULL;
    }

    if(!(chosen = (unsigned char *)malloc(len_L*sizeof(*chosen))))
	return NULL;
    memset(chosen, 0, len_L*sizeof(*chosen));

    if(!seeded) {
	// Bug: should detect error here.
	set_random_seed();
	seeded = 1;
    }

    num_chosen = 0;
    while(num_chosen < n) {
	// Make sure this doesn't overflow.
	j = (int)(floor(rand()/(RAND_MAX+1.0)*len_L));
	if(chosen[j]) 
	    continue;
	chosen[j] = 1;
	num_chosen += 1;
    }

    if(!(py_L_fast = PySequence_Fast(py_L, "broken")))
	goto sample_cleanup;
    if(!(py_retval = PyList_New(n)))
	goto sample_cleanup;
    j = 0;
    for(i=0; i<len_L; i++) {
	if(!chosen[i])
	    continue;
	py_item = PySequence_Fast_GET_ITEM(py_L_fast, i);  /* borrowed REF */
	PyList_SET_ITEM(py_retval, j, py_item);  /* steals reference */
	j++;
	Py_INCREF(py_item);  /* compensate for stolen reference */
	py_item = NULL;
    }
    
 sample_cleanup:
    if(chosen) {
	free(chosen);
    }
    if(py_L_fast) {
	Py_DECREF(py_L_fast);
    }
    if(py_item) {
	Py_DECREF(py_item);
    }
    if(PyErr_Occurred()) {
	if(py_retval) {
	    Py_DECREF(py_retval);
	    py_retval = NULL;
	}
    }
    return py_retval;
}


/* Module definition stuff */

static PyMethodDef cGeneSetMethods[] = {
  {"is_valid_gene_id", cGeneSet_is_valid_gene_id, METH_VARARGS, 
   cGeneSet_is_valid_gene_id__doc__},
  {"select_valid_gene_id", cGeneSet_select_valid_gene_id, METH_VARARGS, 
   cGeneSet_select_valid_gene_id__doc__},
  {"select_valid_gene_id_I", cGeneSet_select_valid_gene_id_I, METH_VARARGS, 
   cGeneSet_select_valid_gene_id_I__doc__},
  {"sample", cGeneSet_sample, METH_VARARGS, 
   cGeneSet_sample__doc__},
  {NULL, NULL, 0, NULL}
};

static char cGeneSet__doc__[] =
"This module provides optimized replacement functions.\n\
";

PyMODINIT_FUNC initcGeneSet(void)
{
  (void) Py_InitModule3("cGeneSet", cGeneSetMethods, cGeneSet__doc__);
}
