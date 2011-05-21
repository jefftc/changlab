/* cMatrixmodule.c
 * 080715  created
 */

#include "Python.h"
//#include "uthash.h"

#if PYTHON_API_VERSION <= 1012  /* Python 2.3 is 1012.  */
#define Py_ssize_t int
#endif

#define ERR_SIZE 1024
char *ERR_STRING[ERR_SIZE];

/* Returns -1 on error. */
int normalize_index(int i, int length) 
{
    if(i < 0)
	i += length;
    if(i < 0 || i >= length)
	return -1;
    return i;
}

static char cMatrix_normalize_index__doc__[] = 
"XXX\n";

static PyObject *cMatrix_normalize_index(PyObject *self, PyObject *args)
{
    int i, length;
    int i_old;

    if(!PyArg_ParseTuple(args, "ii", &i, &length))
	return NULL;

    i_old = i;
    if((i = normalize_index(i, length)) == -1) {
	snprintf((char *)ERR_STRING, ERR_SIZE, 
		 "matrix index out of range [%d:%d]", i_old, length);
	PyErr_SetString(PyExc_IndexError, (char *)ERR_STRING);
	return NULL;
    }
    return PyInt_FromLong(i);
}


static char cMatrix_normalize_indexes__doc__[] = 
"XXX\n";

static PyObject *cMatrix_normalize_indexes(
    PyObject *self, PyObject *args)
{
    PyObject *py_L;
    Py_ssize_t len_L;
    int i, length;

    PyObject *py_item;
    int item_old, item;
    PyObject *py_retval;

    py_item = NULL;
    py_retval = NULL;

    if(!PyArg_ParseTuple(args, "Oi", &py_L, &length))
	return NULL;
    if(!PySequence_Check(py_L)) {
	PyErr_SetString(PyExc_AssertionError, "arg must be a sequence");
	goto normalize_indexes_cleanup;
    }
    if((len_L = PySequence_Size(py_L)) == -1) {
	PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
	goto normalize_indexes_cleanup;
    }

    if(!(py_retval = PyList_New(len_L)))
	goto normalize_indexes_cleanup;

    for(i=0; i<len_L; i++) {
	/* Pull the item out of the list. */
	if(!(py_item = PySequence_GetItem(py_L, i)))  /* New reference. */
	    goto normalize_indexes_cleanup;
	item = (int)PyInt_AsLong(py_item);
	if(item == -1 && PyErr_Occurred())
	    goto normalize_indexes_cleanup;
	Py_DECREF(py_item);
	py_item = NULL;

	/* Normalize the item. */
	item_old = item;
	if((item = normalize_index(item, length)) == -1) {
	    snprintf((char *)ERR_STRING, ERR_SIZE, 
		     "matrix index out of range [%d:%d]", item_old, length);
	    PyErr_SetString(PyExc_IndexError, (char *)ERR_STRING);
	    //PyErr_SetString(PyExc_AssertionError, "index out of range");
	    goto normalize_indexes_cleanup;
	}

	/* Set the item into py_retval. */
	if(!(py_item = PyInt_FromLong(item)))
	    goto normalize_indexes_cleanup;
	PyList_SET_ITEM(py_retval, i, py_item);  /* Steals reference. */
	py_item = NULL;
    }

 normalize_indexes_cleanup:
    if(py_item) { Py_DECREF(py_item); }
    if(PyErr_Occurred()) {
	if(py_retval) { Py_DECREF(py_retval); }
	return NULL;
    }
    return py_retval;
}

//#define MAX_INDEXES 64
//struct word2indexes_struct {
//    char *word;
//    Py_ssize_t indexes[MAX_INDEXES];
//    int num_indexes;
//    UT_hash_handle hh;
//};

static char cMatrix__index_strings__doc__[] = 
"XXX\n";

// This is slower than the pure Python implementation.
//static PyObject *cMatrix__index_strings_slow(
//    PyObject *self, PyObject *args)
//{
//    PyObject *py_words;
//    Py_ssize_t len_words;
//    PyObject *py_words_fast;
//
//    struct word2indexes_struct *word2indexes;
//
//    Py_ssize_t i;
//    PyObject *py_word;
//    char *word;
//    Py_ssize_t word_len;
//    struct word2indexes_struct *w2i, *tmp;
//
//    PyObject *py_retval;
//    PyObject *py_value, *py_item;
//
//    py_words_fast = NULL;
//    word2indexes = NULL;
//    py_retval = NULL;
//    py_value = NULL;
//    py_item = NULL;
//
//    if(!PyArg_ParseTuple(args, "O", &py_words))
//	return NULL;
//
//    /* Set up the words variable. */
//    if(!PySequence_Check(py_words)) {
//	PyErr_SetString(PyExc_AssertionError, "arg must be a sequence");
//	goto _index_strings_cleanup;
//    }
//    if((len_words = PySequence_Size(py_words)) == -1) {
//	PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
//	goto _index_strings_cleanup;
//    }
//    if(!(py_words_fast = PySequence_Fast(py_words, "bad fast sequence")))
//	goto _index_strings_cleanup;
//
//    for(i=0; i<len_words; i++) {
//	//printf("%d.  %d\n", (int)i, HASH_COUNT(word2indexes));
//	/* Borrowed ref. */
//	if(!(py_word = PySequence_Fast_GET_ITEM(py_words, i)))
//	    goto _index_strings_cleanup;
//	if(!PyString_Check(py_word)) {
//	    PyErr_SetString(PyExc_AssertionError, "list must contain strings");
//	    goto _index_strings_cleanup;
//	}
//	if(PyString_AsStringAndSize(py_word, &word, &word_len) == -1)
//	    goto _index_strings_cleanup;
//
//	HASH_FIND_STR(word2indexes, word, w2i);
//	if(w2i == NULL) {
//	    //printf("'%s' not found.\n", word);
//	    /* Was not found in the hash. */
//	    if(!(w2i = (struct word2indexes_struct *)malloc(sizeof(*w2i)))) {
//		PyErr_SetString(PyExc_MemoryError, "out of memory");
//		goto _index_strings_cleanup;
//	    }
//	    w2i->word = word;
//	    w2i->num_indexes = 0;
//	    HASH_ADD_KEYPTR(hh, word2indexes, w2i->word, word_len, w2i);
//	}
//	if(w2i->num_indexes >= MAX_INDEXES) {
//	    PyErr_SetString(PyExc_AssertionError, "too many indexes");
//	    goto _index_strings_cleanup;
//	}
//	w2i->indexes[w2i->num_indexes++] = i;
//    }
//
//    if(!(py_retval = PyDict_New()))
//	goto _index_strings_cleanup;
//    HASH_ITER(hh, word2indexes, w2i, tmp) {
//	//printf("ITER %s %d\n", w2i->word, HASH_COUNT(word2indexes));
//	if(!(py_value = PyList_New((Py_ssize_t)w2i->num_indexes)))
//	    goto _index_strings_cleanup;
//	for(i=0; i<w2i->num_indexes; i++) {
//	    if(!(py_item = PyInt_FromLong(w2i->indexes[i])))
//		goto _index_strings_cleanup;
//	    /* Steals reference. */
//	    PyList_SET_ITEM(py_value, i, py_item);
//	    py_item = NULL;
//	}
//	if(PyDict_SetItemString(py_retval, w2i->word, py_value) == -1)
//	    goto _index_strings_cleanup;
//	Py_DECREF(py_value);
//	py_value = NULL;
//    }
//    
//  _index_strings_cleanup:
//    if(word2indexes) { 
//	HASH_ITER(hh, word2indexes, w2i, tmp) {
//	    HASH_DEL(word2indexes, w2i);
//	    free(w2i);
//	}
//	HASH_CLEAR(hh, word2indexes);
//    }
//    if(py_words_fast) { Py_DECREF(py_words_fast); }
//    if(py_value) { Py_DECREF(py_value); }
//    if(py_item) { Py_DECREF(py_item); }
//    if(PyErr_Occurred()) {
//	if(py_retval) { Py_DECREF(py_retval); }
//	return NULL;
//    }
//    return py_retval;
//}

// This is a little faster than the pure Python implementation.
static PyObject *cMatrix__index_strings(
    PyObject *self, PyObject *args)
{
    PyObject *py_words;
    Py_ssize_t len_words;
    PyObject *py_words_fast;

    PyObject *py_word2indexes;

    Py_ssize_t i;
    PyObject *py_word, *py_indexes;
    PyObject *py_item;
    //int contains;

    py_words_fast = NULL;
    py_word2indexes = NULL;
    py_indexes = NULL;
    py_item = NULL;

    if(!PyArg_ParseTuple(args, "O", &py_words))
	return NULL;

    /* Set up the words variable. */
    if(!PySequence_Check(py_words)) {
	PyErr_SetString(PyExc_AssertionError, "arg must be a sequence");
	goto _index_strings_cleanup;
    }
    if((len_words = PySequence_Size(py_words)) == -1) {
	PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
	goto _index_strings_cleanup;
    }
    if(!(py_words_fast = PySequence_Fast(py_words, "bad fast sequence")))
	goto _index_strings_cleanup;

    if(!(py_word2indexes = PyDict_New()))
	goto _index_strings_cleanup;
    for(i=0; i<len_words; i++) {
	/* Borrowed ref. */
	if(!(py_word = PySequence_Fast_GET_ITEM(py_words, i)))
	    goto _index_strings_cleanup;
	if(!PyString_Check(py_word)) {
	    PyErr_SetString(PyExc_AssertionError, "list must contain strings");
	    goto _index_strings_cleanup;
	}

	if(!(py_item = PyInt_FromLong((long)i)))
	    goto _index_strings_cleanup;
	
	//contains = PyDict_Contains(py_word2indexes, py_word);
	//if(contains == -1)
	//    goto _index_strings_cleanup;
	//if(!contains) {
	if(!PyDict_GetItem(py_word2indexes, py_word)) {
	    if(!(py_indexes = PyList_New(1)))
		goto _index_strings_cleanup;
	    /* Steals reference to py_item. */
	    PyList_SET_ITEM(py_indexes, 0, py_item);
	    py_item = NULL;

	    if(PyDict_SetItem(py_word2indexes, py_word, py_indexes) == -1)
		goto _index_strings_cleanup;
	    Py_DECREF(py_indexes);
	    py_indexes = NULL;
	} else {
	    /* borrowed ref */
	    if(!(py_indexes = PyDict_GetItem(py_word2indexes, py_word))) {
		PyErr_SetString(PyExc_AssertionError, "missing word");
		goto _index_strings_cleanup;
	    }
	    if(PyList_Append(py_indexes, py_item) == -1) {
		py_indexes = NULL;
		goto _index_strings_cleanup;
	    }
	    py_indexes = NULL;
	    Py_DECREF(py_item);
	    py_item = NULL;
	}
    }

  _index_strings_cleanup:
    if(py_words_fast) { Py_DECREF(py_words_fast); }
    if(py_indexes) { Py_DECREF(py_indexes); }
    if(py_item) { Py_DECREF(py_item); }
    if(PyErr_Occurred()) {
	if(py_word2indexes) { Py_DECREF(py_word2indexes); }
	return NULL;
    }
    return py_word2indexes;
}


/* Module definition stuff */

static PyMethodDef cMatrixMethods[] = {
  {"normalize_index", cMatrix_normalize_index, METH_VARARGS, 
   cMatrix_normalize_index__doc__},
  {"normalize_indexes", cMatrix_normalize_indexes, METH_VARARGS, 
   cMatrix_normalize_indexes__doc__},
  {"_index_strings", cMatrix__index_strings, METH_VARARGS, 
   cMatrix__index_strings__doc__},
  {NULL, NULL, 0, NULL}
};

static char cMatrix__doc__[] =
"This module provides optimized replacement functions.\n\
";

PyMODINIT_FUNC initcMatrix(void)
{
  (void) Py_InitModule3("cMatrix", cMatrixMethods, cMatrix__doc__);
}
