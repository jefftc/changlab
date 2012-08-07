/* cjmathmodule.c
 * 080830  created, from ciofnsmodule, ccoherencefnsmodule.
 */

#include "Python.h"
#include "mathlib.h"


#if PYTHON_API_VERSION <= 1012  /* Python 2.3 is 1012.  */
#define Py_ssize_t int
#endif

#define PI 3.1415926535897931



/* Functions in this module. */

static char cjmath_safe_int__doc__[] = 
"XXX\n";

static PyObject *cjmath_safe_int(PyObject *self, PyObject *args)
{
    PyObject *py_x;
    char *buffer;
    Py_ssize_t length;

    if(!PyArg_ParseTuple(args, "O", &py_x))
	return NULL;
    if(py_x == Py_None) {
	Py_INCREF(Py_None);
	return Py_None;
    }
    if(PyString_Check(py_x)) {
	if(PyString_AsStringAndSize(py_x, &buffer, &length) == -1)
	    return NULL;
	if(length == 0 || 
           strcasecmp(buffer, "na") == 0 ||
           strcasecmp(buffer, "null") == 0 ||
	   strcasecmp(buffer, "nan") == 0) {
	    Py_INCREF(Py_None);
	    return Py_None;
	}
    }

    return PyNumber_Int(py_x);
}

static char cjmath_safe_float__doc__[] = 
"XXX\n";

static PyObject *cjmath_safe_float(PyObject *self, PyObject *args)
{
    PyObject *py_x;
    PyObject *py_nan_str, *py_nan_float;
    char *buffer;
    Py_ssize_t length;

    if(!PyArg_ParseTuple(args, "O", &py_x))
	return NULL;
    if(py_x == Py_None) {
	Py_INCREF(Py_None);
	return Py_None;
    }
    if(PyString_Check(py_x)) {
	if(PyString_AsStringAndSize(py_x, &buffer, &length) == -1)
	    return NULL;
	if(length == 0 || 
           strcasecmp(buffer, "na") == 0 ||
           strcasecmp(buffer, "null") == 0) {
	    Py_INCREF(Py_None);
	    return Py_None;
	}
	if(strcasecmp(buffer, "nan") == 0) {
	    if(!(py_nan_str = PyString_FromString("nan")))
		return NULL;
	    py_nan_float = PyFloat_FromString(py_nan_str, NULL);
	    Py_DECREF(py_nan_str);
	    return py_nan_float;
	}

    }
    return PyNumber_Float(py_x);
}

static char cjmath_mean_list__doc__[] = 
"XXX\n";

static PyObject *cjmath_mean_list(PyObject *self, PyObject *args)
{
    PyObject *py_X;
    float *X;
    int length;
    float m;

    if(!PyArg_ParseTuple(args, "O", &py_X))
	return NULL;
    if(!(py2c_fvector(py_X, &X, &length)))
	return NULL; 
    
    m = mean(X, length);
    free(X);
    return PyFloat_FromDouble((double)m);
}

static char cjmath_var_list__doc__[] = 
"XXX\n";

static PyObject *cjmath_var_list(PyObject *self, PyObject *args)
{
    PyObject *py_X;
    float *X;
    int length;
    float m;
    
    if(!PyArg_ParseTuple(args, "O", &py_X))
	return NULL;
    if(!(py2c_fvector(py_X, &X, &length)))
	return NULL;
    
    m = var(X, length);
    free(X);
    return PyFloat_FromDouble((double)m);
}

static char cjmath_cov_matrix__doc__[] = 
"XXX\n";

static PyObject *cjmath_cov_matrix(PyObject *self, PyObject *args, 
    PyObject *keywds)
{
    PyObject *py_X, *py_cov;
    float *X, *cov;
    int nrow, ncol;
    int byrow = 1;

    static char *kwlist[] = {"X", "byrow", NULL};

    X = cov = NULL;
    py_cov = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, keywds, "O|i", kwlist, 
        &py_X, &byrow))
	return NULL;
    
    /* Convert X to a C matrix. */
    if(!py2c_fmatrix(py_X, &X, &nrow, &ncol))
	return NULL;
    if(!byrow) {
	PyErr_SetString(PyExc_NotImplementedError, "col cov not implemented");
    } else {
	if((cov = cov_byrow(X, nrow, ncol)))
	    py_cov = c2py_fmatrix(cov, nrow, nrow);
    }

    if(cov) { free(cov); }
    if(X) { free(X); }
    return py_cov;
}

static char cjmath_cor_matrix__doc__[] = 
"XXX\n";

static PyObject *cjmath_cor_matrix(PyObject *self, PyObject *args, 
    PyObject *keywds)
{
    PyObject *py_X, *py_cor;
    float *X, *cor;
    //double *X, *cor;
    int i, j;
    int nrow, ncol;
    int byrow = 1;
    int safe = 0;
    int abs_ = 0;

    static char *kwlist[] = {"X", "byrow", "safe", "abs", NULL};

    X = cor = NULL;
    py_cor = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, keywds, "O|iii", kwlist, 
				    &py_X, &byrow, &safe, &abs_))
	return NULL;

    /* Convert X to a C matrix. */
    if(!py2c_fmatrix(py_X, &X, &nrow, &ncol))
	goto _cor_matrix_cleanup;
    if(!byrow) {
	PyErr_SetString(PyExc_NotImplementedError, "col cor not implemented");
	goto _cor_matrix_cleanup;
    } 

    if(!(cor = cor_byrow_f(X, nrow, ncol, safe))) {
	PyErr_SetString(PyExc_AssertionError, "cor error");
	goto _cor_matrix_cleanup;
    }
    if(abs_) {
	for(i=0; i<nrow; i++) {
	    for(j=i+1; j<nrow; j++) {
		cor[i*nrow+j] = fabs(cor[i*nrow+j]);
	    }
	}
    }

    py_cor = c2py_fmatrix(cor, nrow, nrow);
    //py_cor = c2py_dmatrix(cor, nrow, nrow);

 _cor_matrix_cleanup:
    if(cor) { free(cor); }
    if(X) { free(X); }
    return py_cor;
}


static char cjmath_log_add__doc__[] = 
"_logadd(logx, logy)\n";

static PyObject *cjmath_log_add(PyObject *self, PyObject *args)
{
    double logx, logy;
    double retval;

    if(!PyArg_ParseTuple(args, "dd", &logx, &logy))
	return NULL;
    if(logy - logx > 100) {
	retval = logy;
    } else if(logx - logy > 100) {
	retval = logx;
    } else {
	double minxy = (logx < logy) ? logx : logy;
	retval = minxy + log(exp(logx-minxy) + exp(logy-minxy));
    }
    return PyFloat_FromDouble(retval);
}

static char cjmath_fisher_z_item__doc__[] = 
"XXX\n";

static PyObject *cjmath_fisher_z_item(PyObject *self, PyObject *args)
{
    double R;
    int N;
    double z;

    if(!PyArg_ParseTuple(args, "di", &R, &N))
	return NULL;
    if(R < -1 || R > 1) {
	PyErr_SetString(PyExc_AssertionError, "R out of range.");
	return NULL;
    }
    if(N < 2) {
	PyErr_SetString(PyExc_AssertionError, "N must be >= 2.");
	return NULL;
    }

    z = fisher_z(R, N);
    return PyFloat_FromDouble(z);
}

static char cjmath_fisher_z_matrix__doc__[] = 
"XXX\n";

static PyObject *cjmath_fisher_z_matrix(PyObject *self, PyObject *args)
{
    PyObject *py_cor, *py_zcor;
    int N;
    float *cor, *zcor;
    int nrow, ncol;

    float z;

    int i, j;
    float r;
    float MAX, MIN;

    py_zcor = NULL;
    cor = zcor = NULL;
    if(!PyArg_ParseTuple(args, "Oi", &py_cor, &N))
	return NULL;
    if(!py2c_fmatrix(py_cor, &cor, &nrow, &ncol))
	return NULL;
    if(nrow != ncol) {
	PyErr_SetString(PyExc_AssertionError, "Invalid matrix");
	goto fisher_z_transform_cleanup;
    }
    if(N < 2) {
	PyErr_SetString(PyExc_AssertionError, "N must be at least 2.");
	goto fisher_z_transform_cleanup;
    }

    /* Do a fisher z-transform on the correlations. */
    if(!(zcor = (float *)malloc(nrow*ncol*sizeof(*zcor)))) {
	PyErr_SetString(PyExc_MemoryError, "Out of memory");
	goto fisher_z_transform_cleanup;
    }
    MAX = 1-1E-10;
    MIN = -1+1E-10;
    for(i=0; i<nrow; i++) {
	for(j=i; j<ncol; j++) {
	    r = cor[i*ncol+j];
	    if(r > MAX) r = MAX;
	    if(r < MIN) r = MIN;
	    z = (float)fisher_z((double)r, N);
	    zcor[i*ncol+j] = zcor[j*ncol+i] = z;
	}
    }

    if(!(py_zcor = c2py_fmatrix(zcor, nrow, ncol)))
	goto fisher_z_transform_cleanup;

 fisher_z_transform_cleanup:
    if(cor) { free(cor); }
    if(zcor) { free(zcor); }
    if(PyErr_Occurred()) {
	if(py_zcor) { Py_DECREF(py_zcor); }
	return NULL;
    }
    return py_zcor;
}

static char cjmath_dnorm_item__doc__[] = 
"XXX\n";

static PyObject *cjmath_dnorm_item(PyObject *self, PyObject *args, 
				   PyObject *keywds)
{
    PyObject *py_x;
    double x;
    double mean = 0.0, sd = 1.0;

    static char *kwlist[] = {"x", "mean", "sd", NULL};

    if(!PyArg_ParseTupleAndKeywords(args, keywds, "O|dd", kwlist,
				    &py_x, &mean, &sd))
	return NULL;
    if(sd <= 0) {
	PyErr_SetString(PyExc_AssertionError, "sd must be positive");
	return NULL;
    }

    x = double_FromPyObject(py_x);
    if(PyErr_Occurred())
	return NULL;
    return PyFloat_FromDouble(dnorm(x, mean, sd));
}

static char cjmath_dparzen__doc__[] = 
"XXX\n";

static PyObject *cjmath_dparzen(PyObject *self, PyObject *args, 
    PyObject *keywds)
{
    float x, h;
    PyObject *py_X_obs, *py_X_count;
    float *X_obs;
    int *X_count;
    int X_obs_len, X_count_len;

    float density;
    PyObject *py_d;

    static char *kwlist[] = {"x", "X_obs", "X_count", "h", NULL};

    h = 1.0;
    py_X_obs = NULL; py_X_count = NULL;
    X_obs = NULL; X_count = NULL;
    py_d = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, keywds, "fOO|f", kwlist, 
        &x, &py_X_obs, &py_X_count, &h))
	return NULL;
    
    /* Convert to C vectors. */
    if(!py2c_fvector(py_X_obs, &X_obs, &X_obs_len))
	goto _dparzen_cleanup;
    if(!py2c_ivector(py_X_count, &X_count, &X_count_len))
	goto _dparzen_cleanup;
    if(X_obs_len != X_count_len) {
	PyErr_SetString(PyExc_AssertionError, "unaligned lists");
	goto _dparzen_cleanup;
    }
    density = dparzen(x, X_obs, X_count, X_count_len, h);
    py_d = PyFloat_FromDouble((double)density);

 _dparzen_cleanup:
    if(X_obs) { free(X_obs); }
    if(X_count) { free(X_count); }
    return py_d;
}


static char cjmath_int_simpson__doc__[] = 
"XXX\n";

static PyObject *cjmath_int_simpson(PyObject *self, PyObject *args)
{
    PyObject *py_fn;
    PyObject *py_arglist, *py_result;
    float h3;
    float a, b, c;
    float f_a, f_b, f_c;

    int i;
    float x[3];
    float f_x[3];
    float s;

    s = 0;
    py_arglist = py_result = NULL;
    if(!PyArg_ParseTuple(args, "Off", &py_fn, &a, &b))
	return NULL;
    if(!PyCallable_Check(py_fn)) {
	PyErr_SetString(PyExc_TypeError, "fn must be callable");
	return NULL;
    }

    c = (a+b)/2.0;
    h3 = fabs(b-a)/6.0;

    x[0] = a; x[1] = b; x[2] = c;
    for(i=0; i<3; i++) {
	if(!(py_arglist = Py_BuildValue("(f)", x[i])))
	    goto int_simpson_cleanup;
	if(!(py_result = PyEval_CallObject(py_fn, py_arglist)))
	    goto int_simpson_cleanup;
	f_x[i] = (float)double_FromPyObject(py_result);
	if(PyErr_Occurred())
	    goto int_simpson_cleanup;
	Py_DECREF(py_arglist); py_arglist=NULL;
	Py_DECREF(py_result); py_result=NULL;
    }
    f_a = f_x[0]; f_b = f_x[1]; f_c = f_x[2];

    s = h3 * (f_a + 4*f_c + f_b);
    /*printf("CSTUFF %g %g %g %g %g %g %g %g\n", 
      a, b, c, f_a, f_b, f_c, s, h3);*/

 int_simpson_cleanup:
    if(py_arglist) { Py_DECREF(py_arglist); }
    if(py_result) { Py_DECREF(py_result); }
    if(PyErr_Occurred()) { return NULL; }
    return PyFloat_FromDouble((double)s);
}

float call_simpson_fn(PyObject *py_fn, float x)
{
    PyObject *py_arglist, *py_result;
    float f_x;

    py_arglist = py_result = NULL;
    f_x = 0;

    if(!(py_arglist = Py_BuildValue("(f)", x)))
	goto call_simpson_fn_cleanup;
    if(!(py_result = PyEval_CallObject(py_fn, py_arglist)))
	goto call_simpson_fn_cleanup;
    f_x = (float)double_FromPyObject(py_result);
    if(PyErr_Occurred())
	goto call_simpson_fn_cleanup;
    
 call_simpson_fn_cleanup:
    if(py_arglist) { Py_DECREF(py_arglist); }
    if(py_result) { Py_DECREF(py_result); }
    return f_x;
}

static char cjmath_int_simpson_adaptive__doc__[] = 
"XXX\n";

#define MAX_STACK 100
static PyObject *cjmath_int_simpson_adaptive(
  PyObject *self, PyObject *args)
{
    PyObject *py_fn;
    float eps;
    float a, b, c, d, e;
    float f_a, f_b, f_c, f_d, f_e;
    float s, left, right, whole;
    float sum;
    float stack[MAX_STACK][8];
    int i, stack_size;

    sum = 0.0;
    if(!PyArg_ParseTuple(args, "Offf", &py_fn, &a, &b, &eps))
	return NULL;
    if(!PyCallable_Check(py_fn)) {
	PyErr_SetString(PyExc_TypeError, "fn must be callable");
	return NULL;
    }

    c = (a+b)/2.0;
    f_a = call_simpson_fn(py_fn, a); 
    f_b = call_simpson_fn(py_fn, b); 
    f_c = call_simpson_fn(py_fn, c);
    s = (b-a)/6.0*(f_a + 4*f_c + f_b);
    if(PyErr_Occurred())
	goto int_simpson_adaptive_cleanup;

    stack[0][0] = a; stack[0][1] = c; stack[0][2] = b;
    stack[0][3] = f_a; stack[0][4] = f_c; stack[0][5] = f_b;
    stack[0][6] = eps; stack[0][7] = s;
    stack_size = 1;

    sum = 0.0;
    while(stack_size) {
	/* a  b  c  d  e */
	/* Pop the stack. */
	i = stack_size-1;
	a = stack[i][0]; c = stack[i][1]; e = stack[i][2];
	f_a = stack[i][3]; f_c = stack[i][4]; f_e = stack[i][5];
	eps = stack[i][6]; whole = stack[i][7];
	stack_size--;

	b = (a+c)/2.0; d = (c+e)/2.0;
	f_b = call_simpson_fn(py_fn, b);
	f_d = call_simpson_fn(py_fn, d);

	left = (c-a)/6.0 * (f_a + 4*f_b + f_c);
	right = (e-c)/6.0 * (f_c + 4*f_d + f_e);
	if(fabs(left+right-whole) <= 15*eps) {
	    sum += left+right+(left+right-whole)/15.0;
	    continue;
	}

	/* Push left and right onto the stack. */
        /* Check for stack overflow. */
        if(stack_size >= MAX_STACK-2) {
            PyErr_SetString(PyExc_AssertionError, "stack overflow");
            goto int_simpson_adaptive_cleanup;
        }
	i = stack_size;
        stack[i][0] = a; stack[i][1] = b; stack[i][2] = c;
        stack[i][3] = f_a; stack[i][4] = f_b; stack[i][5] = f_c;
        stack[i][6] = eps/2.0; stack[i][7] = left;
	stack_size++;
	
	i = stack_size;
        stack[i][0] = c; stack[i][1] = d; stack[i][2] = e;
        stack[i][3] = f_c; stack[i][4] = f_d; stack[i][5] = f_e;
        stack[i][6] = eps/2.0; stack[i][7] = right;
	stack_size++;
    }

 int_simpson_adaptive_cleanup:
    if(PyErr_Occurred()) { return NULL; }
    return PyFloat_FromDouble((double)sum);
}

static char cjmath_pparzen__doc__[] = 
"XXX\n";

static PyObject *cjmath_pparzen(PyObject *self, 
    PyObject *args, PyObject *keywds)
{
    float p;
    PyObject *py_X_obs, *py_X_count;
    float *X_obs;
    int *X_count;
    int X_obs_len, X_count_len;
    float h, eps;

    float a, b, c, d, e;
    float f_a, f_b, f_c, f_d, f_e;
    float s, left, right, whole;
    float sum;
    float stack[MAX_STACK][8];
    int i, stack_size;

    int N, max_X_count, min_X_obs;
    float A;

    static char *kwlist[] = {"p", "X_obs", "X_count", "h", "eps", NULL};

    py_X_obs = NULL; py_X_count = NULL;
    X_obs = NULL; X_count = NULL;
    h = 1.0;
    eps = 1E-10;    /* Maximum error for dparzen estimate. */
    sum = 0.0;

    if(!PyArg_ParseTupleAndKeywords(args, keywds, "fOO|ff", kwlist, 
        &p, &py_X_obs, &py_X_count, &h, &eps))
	return NULL;
    /* Convert to C vectors. */
    if(!py2c_fvector(py_X_obs, &X_obs, &X_obs_len))
	goto _pparzen_cleanup;
    if(!py2c_ivector(py_X_count, &X_count, &X_count_len))
	goto _pparzen_cleanup;
    if(X_obs_len != X_count_len) {
	PyErr_SetString(PyExc_AssertionError, "unaligned lists");
	goto _pparzen_cleanup;
    }

    /* See the Python source for the calculation of the lower
       bound. */
    N = 0;
    max_X_count = X_count[0];
    min_X_obs = X_obs[0];
    for(i=0; i<X_count_len; i++) {
	N += X_count[i];
	if(X_count[i] > max_X_count)
	    max_X_count = X_count[i];
	if(X_obs[i] < min_X_obs)
	    min_X_obs = X_obs[i];
    }
    A = eps*N*h*sqrt(2*PI) / max_X_count;
    a = min_X_obs - h*sqrt(-2*log(A));
    b = p;

    c = (a+b)/2.0;
    f_a = dparzen_fast(a, X_obs, X_count, X_count_len, h, N, max_X_count, eps);
    f_b = dparzen_fast(b, X_obs, X_count, X_count_len, h, N, max_X_count, eps);
    f_c = dparzen_fast(c, X_obs, X_count, X_count_len, h, N, max_X_count, eps);
    s = (b-a)/6.0*(f_a + 4*f_c + f_b);
    if(PyErr_Occurred())
	goto _pparzen_cleanup;

    stack[0][0] = a; stack[0][1] = c; stack[0][2] = b;
    stack[0][3] = f_a; stack[0][4] = f_c; stack[0][5] = f_b;
    stack[0][6] = eps; stack[0][7] = s;
    stack_size = 1;

    sum = 0.0;
    while(stack_size) {
	/* a  b  c  d  e */
	/* Pop the stack. */
	i = stack_size-1;
	a = stack[i][0]; c = stack[i][1]; e = stack[i][2];
	f_a = stack[i][3]; f_c = stack[i][4]; f_e = stack[i][5];
	eps = stack[i][6]; whole = stack[i][7];
	stack_size--;

	b = (a+c)/2.0; d = (c+e)/2.0;
	f_b = dparzen_fast(b, X_obs, X_count, X_count_len, h, 
            N, max_X_count, eps);
	f_d = dparzen_fast(d, X_obs, X_count, X_count_len, h,
            N, max_X_count, eps);
	left = (c-a)/6.0 * (f_a + 4*f_b + f_c);
	right = (e-c)/6.0 * (f_c + 4*f_d + f_e);
	if(fabs(left+right-whole) <= 15*eps) {
	    sum += left+right+(left+right-whole)/15.0;
	    continue;
	}

	/* Push left and right onto the stack. */
        /* Check for stack overflow. */
	/*printf("Stack: %d %g %g %g %g %g\n", stack_size, left, right, whole, 
	  fabs(left+right-whole), 15*eps); */
        if(stack_size >= MAX_STACK-2) {
            PyErr_SetString(PyExc_AssertionError, "stack overflow");
            goto _pparzen_cleanup;
        }
	i = stack_size;
        stack[i][0] = a; stack[i][1] = b; stack[i][2] = c;
        stack[i][3] = f_a; stack[i][4] = f_b; stack[i][5] = f_c;
        stack[i][6] = eps/2.0; stack[i][7] = left;
	stack_size++;
	
	i = stack_size;
        stack[i][0] = c; stack[i][1] = d; stack[i][2] = e;
        stack[i][3] = f_c; stack[i][4] = f_d; stack[i][5] = f_e;
        stack[i][6] = eps/2.0; stack[i][7] = right;
	stack_size++;
    }

 _pparzen_cleanup:
    if(X_obs) { free(X_obs); }
    if(X_count) { free(X_count); }
    if(PyErr_Occurred()) { return NULL; }
    return PyFloat_FromDouble((double)sum);
}

static char cjmath__pparzen_ul__doc__[] = 
"XXX\n";

static PyObject *cjmath__pparzen_ul(PyObject *self, 
    PyObject *args, PyObject *keywds)
{
    float p_l, p_u;
    PyObject *py_X_obs, *py_X_count;
    float *X_obs;
    int *X_count;
    int X_obs_len, X_count_len;
    float h, eps;

    float a, b, c, d, e;
    float f_a, f_b, f_c, f_d, f_e;
    float s, left, right, whole;
    float sum;
    float stack[MAX_STACK][8];
    int i, stack_size;

    int N, max_X_count, min_X_obs;

    static char *kwlist[] = {"p_l", "p_u", "X_obs", "X_count", "h", "eps", 
			     NULL};

    py_X_obs = NULL; py_X_count = NULL;
    X_obs = NULL; X_count = NULL;
    h = 1.0;
    eps = 1E-10;    /* Maximum error for dparzen estimate. */
    sum = 0.0;

    if(!PyArg_ParseTupleAndKeywords(args, keywds, "ffOO|ff", kwlist, 
	&p_l, &p_u, &py_X_obs, &py_X_count, &h, &eps))
	return NULL;
    /* Convert to C vectors. */
    if(!py2c_fvector(py_X_obs, &X_obs, &X_obs_len))
	goto _pparzen_ul_cleanup;
    if(!py2c_ivector(py_X_count, &X_count, &X_count_len))
	goto _pparzen_ul_cleanup;
    if(X_obs_len != X_count_len) {
	PyErr_SetString(PyExc_AssertionError, "unaligned lists");
	goto _pparzen_ul_cleanup;
    }

    N = 0;
    max_X_count = X_count[0];
    min_X_obs = X_obs[0];
    for(i=0; i<X_count_len; i++) {
	N += X_count[i];
	if(X_count[i] > max_X_count)
	    max_X_count = X_count[i];
	if(X_obs[i] < min_X_obs)
	    min_X_obs = X_obs[i];
    }
    a = p_l;
    b = p_u;

    c = (a+b)/2.0;
    f_a = dparzen_fast(a, X_obs, X_count, X_count_len, h, N, max_X_count, eps);
    f_b = dparzen_fast(b, X_obs, X_count, X_count_len, h, N, max_X_count, eps);
    f_c = dparzen_fast(c, X_obs, X_count, X_count_len, h, N, max_X_count, eps);
    s = (b-a)/6.0*(f_a + 4*f_c + f_b);
    if(PyErr_Occurred())
	goto _pparzen_ul_cleanup;

    stack[0][0] = a; stack[0][1] = c; stack[0][2] = b;
    stack[0][3] = f_a; stack[0][4] = f_c; stack[0][5] = f_b;
    stack[0][6] = eps; stack[0][7] = s;
    stack_size = 1;

    sum = 0.0;
    while(stack_size) {
	/* a  b  c  d  e */
	/* Pop the stack. */
	i = stack_size-1;
	a = stack[i][0]; c = stack[i][1]; e = stack[i][2];
	f_a = stack[i][3]; f_c = stack[i][4]; f_e = stack[i][5];
	eps = stack[i][6]; whole = stack[i][7];
	stack_size--;

	b = (a+c)/2.0; d = (c+e)/2.0;
	f_b = dparzen_fast(b, X_obs, X_count, X_count_len, h, 
            N, max_X_count, eps);
	f_d = dparzen_fast(d, X_obs, X_count, X_count_len, h,
            N, max_X_count, eps);
	left = (c-a)/6.0 * (f_a + 4*f_b + f_c);
	right = (e-c)/6.0 * (f_c + 4*f_d + f_e);
	if(fabs(left+right-whole) <= 15*eps) {
	    sum += left+right+(left+right-whole)/15.0;
	    continue;
	}

	/* Push left and right onto the stack. */
        /* Check for stack overflow. */
	/*printf("Stack: %d %g %g %g %g %g\n", stack_size, left, right, whole, 
	  fabs(left+right-whole), 15*eps); */
        if(stack_size >= MAX_STACK-2) {
            PyErr_SetString(PyExc_AssertionError, "stack overflow");
            goto _pparzen_ul_cleanup;
        }
	i = stack_size;
        stack[i][0] = a; stack[i][1] = b; stack[i][2] = c;
        stack[i][3] = f_a; stack[i][4] = f_b; stack[i][5] = f_c;
        stack[i][6] = eps/2.0; stack[i][7] = left;
	stack_size++;
	
	i = stack_size;
        stack[i][0] = c; stack[i][1] = d; stack[i][2] = e;
        stack[i][3] = f_c; stack[i][4] = f_d; stack[i][5] = f_e;
        stack[i][6] = eps/2.0; stack[i][7] = right;
	stack_size++;
    }

 _pparzen_ul_cleanup:
    if(X_obs) { free(X_obs); }
    if(X_count) { free(X_count); }
    if(PyErr_Occurred()) { return NULL; }
    return PyFloat_FromDouble((double)sum);
}

/* Module definition stuff */

static PyMethodDef CJMathMethods[] = {
  {"safe_int", cjmath_safe_int, METH_VARARGS, cjmath_safe_int__doc__},
  {"safe_float", cjmath_safe_float, METH_VARARGS, cjmath_safe_float__doc__},
  {"mean_list", cjmath_mean_list, METH_VARARGS, cjmath_mean_list__doc__},
  {"var_list", cjmath_var_list, METH_VARARGS, cjmath_var_list__doc__},
  {"cov_matrix", (PyCFunction)cjmath_cov_matrix, METH_VARARGS | METH_KEYWORDS, 
   cjmath_cov_matrix__doc__},
  {"cor_matrix", (PyCFunction)cjmath_cor_matrix, METH_VARARGS | METH_KEYWORDS, 
   cjmath_cor_matrix__doc__},
  {"log_add", (PyCFunction)cjmath_log_add, METH_VARARGS, 
   cjmath_log_add__doc__},
  {"fisher_z_item", cjmath_fisher_z_item, METH_VARARGS, 
   cjmath_fisher_z_item__doc__},
  {"fisher_z_matrix", cjmath_fisher_z_matrix, METH_VARARGS, 
   cjmath_fisher_z_matrix__doc__},
  {"dnorm_item", (PyCFunction)cjmath_dnorm_item, METH_VARARGS|METH_KEYWORDS, 
   cjmath_dnorm_item__doc__},
  {"dparzen", (PyCFunction)cjmath_dparzen, METH_VARARGS | METH_KEYWORDS, 
   cjmath_dparzen__doc__},
  {"pparzen", (PyCFunction)cjmath_pparzen, METH_VARARGS | METH_KEYWORDS, 
   cjmath_pparzen__doc__},
  {"_pparzen_ul", (PyCFunction)cjmath__pparzen_ul, 
   METH_VARARGS | METH_KEYWORDS, cjmath__pparzen_ul__doc__},
  {"int_simpson", cjmath_int_simpson, METH_VARARGS, cjmath_int_simpson__doc__},
  {"int_simpson_adaptive", cjmath_int_simpson_adaptive, METH_VARARGS, 
   cjmath_int_simpson_adaptive__doc__},
  {NULL, NULL, 0, NULL},
};

static char cjmath__doc__[] =
"This module provides optimized replacement functions.\n\
";

PyMODINIT_FUNC initcjmath(void)
{
  (void) Py_InitModule3("cjmath", CJMathMethods, cjmath__doc__);
}
