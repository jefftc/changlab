/* cbie3module.c
 * 150826  created
 */

#include "Python.h"
#include <math.h>


/* Functions in this module. */

struct int_list_struct {
    int len;
    int *values;
};

struct int_list_struct *new_int_list(int length) {
    int *values;
    struct int_list_struct *list;

    if(!(list = (struct int_list_struct *)malloc(
        sizeof(struct int_list_struct))))
        return NULL;
    if(!(values = (int *)malloc(length*sizeof(int)))) {
        free(list);
        return NULL;
    }
    list->len = length;
    list->values = values;
    return list;
}


void free_int_list(struct int_list_struct *list) {
    if(list == NULL)
        return;
    if(list->values)
        free(list->values);
    free(list);
}

void free_list_of_int_list(struct int_list_struct **list_of_lists,
                           int length) {
    int i;
    
    if(!list_of_lists)
        return;
    for(i=0; i<length; i++) {
        if(list_of_lists[i])
            free_int_list(list_of_lists[i]);
    }
    free(list_of_lists);
}


struct int_list_struct **PyTupleList_to_C(PyObject *py_tup_list) {
    // Return a pointer to an array of int_list_structs.  It has the
    // same number of int_list_structs as py_tup_list.
    int i, j;
    int len_tup_list;
    struct int_list_struct **c_tup_list;  // Pointer to pointers to c_tup_list;

    PyObject *py_tuple;
    PyObject *py_number, *py_int;
    int len_tup;
    struct int_list_struct *c_tup;

    len_tup_list = 0;
    c_tup_list = NULL;
    py_tuple = NULL;
    py_number = NULL;
    py_int = NULL;
    c_tup = NULL;

    if(!PySequence_Check(py_tup_list)) {
        PyErr_SetString(PyExc_AssertionError, "tup_list1 must be a sequence");
        goto py_tuplelist_to_c_cleanup;
    }
    if((len_tup_list = PySequence_Size(py_tup_list)) == -1) {
        PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
        goto py_tuplelist_to_c_cleanup;
    }
    i = len_tup_list * sizeof(struct int_list_struct *);
    if(!(c_tup_list = (struct int_list_struct **)malloc(i))) {
        PyErr_SetString(PyExc_AssertionError, "malloc");
        goto py_tuplelist_to_c_cleanup;
    }
    for(i=0; i<len_tup_list; i++)
        c_tup_list[i] = NULL;

    // Copy py_tuple_list to c_tup_list.
    for(i=0; i<len_tup_list; i++) {
        /* New reference. */
        if(!(py_tuple = PySequence_GetItem(py_tup_list, i)))
            goto py_tuplelist_to_c_cleanup;
        if(!PyTuple_Check(py_tuple)) {
            PyErr_SetString(PyExc_AssertionError, "not tuple");
            goto py_tuplelist_to_c_cleanup;
        }
        len_tup = PyTuple_Size(py_tuple);

        if(!(c_tup = new_int_list(len_tup))) {
            PyErr_SetString(PyExc_AssertionError, "memory");
            goto py_tuplelist_to_c_cleanup;
        }
        for(j=0; j<len_tup; j++) {
            // Borrowed reference.
            if(!(py_number = PyTuple_GetItem(py_tuple, j))) {
                goto py_tuplelist_to_c_cleanup;
            }
            if(!PyNumber_Check(py_number)) {
                //PyObject_Print(py_number, stdout, Py_PRINT_RAW);
                PyErr_SetString(PyExc_AssertionError, "not number");
                py_number = NULL;
                goto py_tuplelist_to_c_cleanup;
            }
            if(!(py_int = PyNumber_Int(py_number))) {
                py_number = NULL;
                goto py_tuplelist_to_c_cleanup;
            }
            c_tup->values[j] = (int)PyInt_AS_LONG(py_int);
            py_number = NULL;
            Py_DECREF(py_int);
            py_int = NULL;
        }

        c_tup_list[i] = c_tup;
        c_tup = NULL;
        
        Py_DECREF(py_tuple);
        py_tuple = NULL;
    }


 py_tuplelist_to_c_cleanup:
    if(py_tuple) { Py_DECREF(py_tuple); }
    if(py_number) { Py_DECREF(py_number); }
    if(py_int) { Py_DECREF(py_int); }
    if(c_tup) { free_int_list(c_tup); }
    if(PyErr_Occurred()) {
        if(c_tup_list) {
            free_list_of_int_list(c_tup_list, len_tup_list);
        }
        return NULL;
    }
    return c_tup_list;
}


int compare_int(const void *a, const void *b) {
    int *a_int = (int *)a;
    int *b_int = (int *)b;
    return *a_int - *b_int;
}


static void uniq_int_list(int *values, int *num_values) {
    // Change values in place.
    int i, j;

    // 90% of the calls from bie3._merge_paths have only 2 values.
    // Special code these cases.  The rest only have 3.
    if(*num_values <= 1) {
    } else if(*num_values == 2) {
        if(values[0] == values[1])
            *num_values = 1;
    } else if(*num_values == 3) {
        // Values X1, X2, X3.
        // Compare X2 and X3.
        if(values[1] == values[2])
            *num_values -= 1;
        // Now list is either (X1, X2, X3) or (X1, X2).
        // Compare X1 and X2.
        if(values[0] == values[1]) {
            values[1] = values[2];  // May not be necessary, but doesn't hurt.
            *num_values -= 1;
        }
        // Now (X1, X2, X3), (X1, X3), (X1, X2), (X1).
        // Compare X1 and X3.
        if((*num_values == 3) && (values[0] == values[2]))
            *num_values -= 1;
        else if((*num_values == 2) && (values[0] == values[1]))
            *num_values -= 1;
    } else {
        qsort(values, *num_values, sizeof(int), compare_int);
        i = 0;
        while(i < *num_values-1) {
            if(values[i] == values[i+1]) {
                for(j=i+1; j<*num_values-1; j++)
                    values[j] = values[j+1];
                (*num_values)--;
            } else {
                i++;
            }
        }
    }
}

/* Module definition stuff */

static char cbie3__product_and_chain_h__doc__[] =
    "XXX\n";

#define MAX_VALUES 1024 * 8
static PyObject *cbie3__product_and_chain_h(PyObject *self, PyObject *args)
{
    int i, j, k, l, m;
    PyObject *py_tup_list1, *py_tup_list2;
    int len_tup_list1, len_tup_list2;
    int max_length;

    struct int_list_struct **c_tup_list1, **c_tup_list2;

    //struct int_list_struct *DEBUG;

    PyObject *py_results;
    PyObject *py_tuple;
    PyObject *py_item;

    int len;
    int len1, len2;
    int *values;
    int v;
    //int v1, v2;

    c_tup_list1 = NULL;
    c_tup_list2 = NULL;
    py_results = NULL;
    py_tuple = NULL;
    py_item = NULL;
    values = NULL;
    
    if(!PyArg_ParseTuple(args, "OOi", &py_tup_list1, &py_tup_list2,
                         &max_length))
        return NULL;
    if(!PySequence_Check(py_tup_list1)) {
        PyErr_SetString(PyExc_AssertionError, "tup_list1 must be a sequence");
        goto _product_and_chain_h_cleanup;
    }
    if(!PySequence_Check(py_tup_list2)) {
        PyErr_SetString(PyExc_AssertionError, "tup_list2 must be a sequence");
        goto _product_and_chain_h_cleanup;
    }
    if((len_tup_list1 = PySequence_Size(py_tup_list1)) == -1) {
        PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
        goto _product_and_chain_h_cleanup;
    }
    if((len_tup_list2 = PySequence_Size(py_tup_list2)) == -1) {
        PyErr_SetString(PyExc_AssertionError, "unable to get sequence size");
        goto _product_and_chain_h_cleanup;
    }

    if(!(c_tup_list1 = PyTupleList_to_C(py_tup_list1)))
        goto _product_and_chain_h_cleanup;
    if(!(c_tup_list2 = PyTupleList_to_C(py_tup_list2)))
        goto _product_and_chain_h_cleanup;
    
    if(!(py_results = PyDict_New()))
        goto _product_and_chain_h_cleanup;

    if(!(values = (int *)malloc(MAX_VALUES*sizeof(int)))) {
        PyErr_SetString(PyExc_AssertionError, "out of memory");
        goto _product_and_chain_h_cleanup;
    }
    
    for(i=0; i<len_tup_list1; i++) {
        for(j=0; j<len_tup_list2; j++) {
            // Merge the two lists.
            len1 = c_tup_list1[i]->len;
            len2 = c_tup_list2[j]->len;
            if(len1+len2 > MAX_VALUES) {
                PyErr_SetString(PyExc_AssertionError, "too many values");
                goto _product_and_chain_h_cleanup;
            }

            //// DEBUG
            //for(k=0; k<len1-1; k++) {
            //    v1 = c_tup_list1[i]->values[k];
            //    v2 = c_tup_list1[i]->values[k+1];
            //    if(v1 >= v2) {
            //        PyErr_SetString(PyExc_AssertionError, "out of order 1");
            //        goto _product_and_chain_h_cleanup;
            //    }
            // }
            //// DEBUG
            //for(k=0; k<len2-1; k++) {
            //    v1 = c_tup_list2[j]->values[k];
            //    v2 = c_tup_list2[j]->values[k+1];
            //    if(v1 >= v2) {
            //        PyErr_SetString(PyExc_AssertionError, "out of order 2");
            //        goto _product_and_chain_h_cleanup;
            //    }
            // }

            // The vast majority of the time, len1 and len2 are at the
            // limit.
            
            // Add the members from list1.
            for(k=0; k<len1; k++)
                values[k] = c_tup_list1[i]->values[k];
            len = len1;

            // Add each member from list2 if it's not already in
            // list1.  Assume that list1 and list2 are already sorted.
            // Maintain a sorted order in values.
            l = 0;
            for(k=0; k<len2; k++) {
                // Find the right place to add the new value.
                v = c_tup_list2[j]->values[k];
                for(; l<len; l++) {
                    if(v <= values[l])
                        break;
                }
                // No duplicates.
                if((l < len) && (v == values[l]))
                    continue;
                // Move the other values out of the way.
                for(m=len; m>l; m--)
                    values[m] = values[m-1];
                // Add this one.
                values[l] = v;
                len++;
                if(len > max_length)
                    break;
                //values[len1+k] = c_tup_list2[j]->values[k];
            }

            //// DEBUG
            //for(k=0; k<len-1; k++) {
            //    if(values[k] >= values[k+1]) {
            //        PyErr_SetString(PyExc_AssertionError, "out of order 3");
            //        goto _product_and_chain_h_cleanup;
            //    }
            // }

            //// Sort the values.
            //// Optimize: don't sort if there is <= 2 values?
            //qsort(values, len, sizeof(int), compare_int);
            //  
            //// No duplicate values.  Optimize this.
            //k = 0;
            //while(k < len-1) {
            //    if(values[k] == values[k+1]) {
            //        for(l=k+1; l<len-1; l++)
            //            values[l] = values[l+1];
            //        len -= 1;
            //    } else {
            //        k += 1;
            //    }
            // }

            if(len > max_length)
                continue;
            
            /* Make a tuple of the values. */
            if(!(py_tuple = PyTuple_New(len)))    /* New Reference */
                goto _product_and_chain_h_cleanup;
            for(k=0; k<len; k++) {
                if(!(py_item = PyInt_FromLong(values[k])))
                    goto _product_and_chain_h_cleanup;
                // Returns 0 on success.
                if(PyTuple_SetItem(py_tuple, k, py_item) != 0)
                    goto _product_and_chain_h_cleanup;
                py_item = NULL;  /* Reference stolen by SetItem */
            }

            /* Add these values to py_results. */
            /* Does not steal reference to Py_None.  Will increment itself. */
            if(PyDict_SetItem(py_results, py_tuple, Py_None) != 0)
                goto _product_and_chain_h_cleanup;
            Py_DECREF(py_tuple);
            py_tuple = NULL;
        }
    }
    
 _product_and_chain_h_cleanup:
    if(c_tup_list1) { free_list_of_int_list(c_tup_list1, len_tup_list1); }
    if(c_tup_list2) { free_list_of_int_list(c_tup_list2, len_tup_list2); }
    if(values) { free(values); }
    if(py_tuple) { Py_DECREF(py_tuple); }
    if(py_item) { Py_DECREF(py_item); }
    if(PyErr_Occurred()) {
        if(py_results) {
            Py_DECREF(py_results);
        }
        return NULL;
    }
    return py_results;
}


static char cbie3__uniq_intlist__doc__[] =
    "XXX\n";

static PyObject *cbie3__uniq_intlist(PyObject *self, PyObject *args)
{
    int i;
    PyObject *py_seq, *py_number, *py_int;
    int len_py_seq;
    int *values;
    int num_values;

    PyObject *py_results;
    PyObject *py_item;

    values = NULL;
    py_number = NULL;
    py_int = NULL;
    py_results = NULL;
    py_item = NULL;

    if(!PyArg_ParseTuple(args, "O", &py_seq))
        return NULL;
    if(!PyList_Check(py_seq)) {
        PyErr_SetString(PyExc_AssertionError, "seq must be a list");
        goto _uniq_intlist_cleanup;
    }
    if((len_py_seq = PyList_Size(py_seq)) == -1) {
        PyErr_SetString(PyExc_AssertionError, "unable to get list size");
        goto _uniq_intlist_cleanup;
    }

    if(!(values = (int *)malloc(len_py_seq*sizeof(int)))) {
        PyErr_SetString(PyExc_AssertionError, "out of memory");
        goto _uniq_intlist_cleanup;
    }

    // Extract the list into values.
    num_values = len_py_seq;
    for(i=0; i<len_py_seq; i++) {
        // Borrowed reference.
        if(!(py_number = PyList_GetItem(py_seq, i))) {
            goto _uniq_intlist_cleanup;
        }
        if(!PyNumber_Check(py_number)) {
            //PyObject_Print(py_number, stdout, Py_PRINT_RAW);
            PyErr_SetString(PyExc_AssertionError, "not number");
            py_number = NULL;
            goto _uniq_intlist_cleanup;
        }
        // New reference.
        if(!(py_int = PyNumber_Int(py_number))) {
            py_number = NULL;
            goto _uniq_intlist_cleanup;
        }
        values[i] = (int)PyInt_AS_LONG(py_int);
        py_number = NULL;
        Py_DECREF(py_int);
        py_int = NULL;
    }

    // Get rid of duplicate values.
    if(num_values >= 2) {
        uniq_int_list(values, &num_values);
    }

    // Save the results.
    if(!(py_results = PyList_New(num_values)))
        goto _uniq_intlist_cleanup;
    for(i=0; i<num_values; i++) {
        // New reference.
        if(!(py_item = PyInt_FromLong(values[i])))
            goto _uniq_intlist_cleanup;
        // Returns 0 on success.  Steals reference.
        // SET_ITEM leads to segfault.  Maybe trying to decref NULL?
        if(PyList_SetItem(py_results, i, py_item) != 0)
            goto _uniq_intlist_cleanup;
        py_item = NULL;
    }

 _uniq_intlist_cleanup:
    if(values) { free(values); }
    if(py_number) { Py_DECREF(py_number); }
    if(py_int) { Py_DECREF(py_int); }
    if(py_item) { Py_DECREF(py_item); }
    if(PyErr_Occurred()) {
        if(py_results) {
            Py_DECREF(py_results);
        }
        return NULL;
    }
    return py_results;
}


static char cbie3__uniq_flatten1_intlist__doc__[] =
    "XXX\n";

static PyObject *cbie3__uniq_flatten1_intlist(PyObject *self, PyObject *args)
{
    int i, j;
    PyObject *py_seq, *py_list, *py_number, *py_int;
    int len_py_seq, len_py_list;
    int *values;
    int num_values;

    PyObject *py_results;
    PyObject *py_item;

    values = NULL;
    py_number = NULL;
    py_list = NULL;
    py_int = NULL;
    py_results = NULL;
    py_item = NULL;

    if(!PyArg_ParseTuple(args, "O", &py_seq))
        return NULL;
    if(!PyList_Check(py_seq)) {
        PyErr_SetString(PyExc_AssertionError, "seq must be a list");
        goto _uniq_flatten1_intlist_cleanup;
    }
    if((len_py_seq = PyList_Size(py_seq)) == -1) {
        PyErr_SetString(PyExc_AssertionError, "unable to get list size");
        goto _uniq_flatten1_intlist_cleanup;
    }

    if(!(values = (int *)malloc(MAX_VALUES*sizeof(int)))) {
        PyErr_SetString(PyExc_AssertionError, "out of memory");
        goto _uniq_flatten1_intlist_cleanup;
    }
    num_values = 0;

    // Extract the list into values.
    for(i=0; i<len_py_seq; i++) {
        // Borrowed reference.
        if(!(py_list = PyList_GetItem(py_seq, i))) {
            goto _uniq_flatten1_intlist_cleanup;
        }
        if(!PyList_Check(py_list)) {
            PyErr_SetString(PyExc_AssertionError, "each item must be a list");
            py_list = NULL;
            goto _uniq_flatten1_intlist_cleanup;
        }
        if((len_py_list = PyList_Size(py_list)) == -1) {
            PyErr_SetString(PyExc_AssertionError, "unable to get list size");
            py_list = NULL;
            goto _uniq_flatten1_intlist_cleanup;
        }
        for(j=0; j<len_py_list; j++) {
            // Borrowed reference.
            if(!(py_number = PyList_GetItem(py_list, j))) {
                py_list = NULL;
                goto _uniq_flatten1_intlist_cleanup;
            }
            if(!PyNumber_Check(py_number)) {
                PyErr_SetString(PyExc_AssertionError, "not number");
                py_list = NULL;
                py_number = NULL;
                goto _uniq_flatten1_intlist_cleanup;
            }
            // New reference.
            if(!(py_int = PyNumber_Int(py_number))) {
                py_list = NULL;
                py_number = NULL;
                goto _uniq_flatten1_intlist_cleanup;
            }
            values[num_values++] = (int)PyInt_AS_LONG(py_int);
            py_number = NULL;
            Py_DECREF(py_int);
            py_int = NULL;
        }
        py_list = NULL;
    }

    // Get rid of duplicate values.
    if(num_values >= 2) {
        uniq_int_list(values, &num_values);
    }

    // Save the results.
    if(!(py_results = PyList_New(num_values)))
        goto _uniq_flatten1_intlist_cleanup;
    for(i=0; i<num_values; i++) {
        // New reference.
        if(!(py_item = PyInt_FromLong(values[i])))
            goto _uniq_flatten1_intlist_cleanup;
        // Returns 0 on success.  Steals reference.
        // SET_ITEM leads to segfault.  Maybe trying to decref NULL?
        if(PyList_SetItem(py_results, i, py_item) != 0)
            goto _uniq_flatten1_intlist_cleanup;
        py_item = NULL;
    }

 _uniq_flatten1_intlist_cleanup:
    if(values) { free(values); }
    if(py_list) { Py_DECREF(py_list); }
    if(py_number) { Py_DECREF(py_number); }
    if(py_int) { Py_DECREF(py_int); }
    if(py_item) { Py_DECREF(py_item); }
    if(PyErr_Occurred()) {
        if(py_results) {
            Py_DECREF(py_results);
        }
        return NULL;
    }
    return py_results;
}

//static char cbie3__get_attribute_type__doc__[] =
//    "XXX\n";
//
//#define TYPE_ATOM 100
//#define TYPE_ENUM 101
//static PyObject *cbie3__get_attribute_type(PyObject *self, PyObject *args)
//{
//    PyObject *py_name;
//
//    if(!PyArg_ParseTuple(args, "O", &py_name))
//        return NULL;
//    if(PyString_Check(py_name)) {
//        return PyInt_FromLong(TYPE_ATOM);
//    }
//    if(PyList_Check(py_name) || PyTuple_Check(py_name)) {
//        return PyInt_FromLong(TYPE_ENUM);
//    }
//    PyErr_SetString(PyExc_AssertionError, "Unknown attribute type");
//    return NULL;
//}



static PyMethodDef cbie3Methods[] = {
  {"_product_and_chain_h", cbie3__product_and_chain_h, METH_VARARGS,
   cbie3__product_and_chain_h__doc__},
  {"_uniq_intlist", cbie3__uniq_intlist, METH_VARARGS,
   cbie3__uniq_intlist__doc__},
  {"_uniq_flatten1_intlist", cbie3__uniq_flatten1_intlist, METH_VARARGS,
   cbie3__uniq_flatten1_intlist__doc__},
  // Does not result in much speed-up.
  //{"_get_attribute_type", cbie3__get_attribute_type, METH_VARARGS,
  // cbie3__get_attribute_type__doc__},
  {NULL, NULL, 0, NULL}
};


static char cbie3__doc__[] =
"This module provides optimized replacement functions.\n\
";

PyMODINIT_FUNC initcbie3(void)
{
  (void) Py_InitModule3("cbie3", cbie3Methods, cbie3__doc__);
}
