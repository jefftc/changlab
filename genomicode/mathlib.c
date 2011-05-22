/* mathlib.c
 * 080831  created
 */

#include <Python.h>
//#include <pmmintrin.h>  // For __mm128 stuff.
#include "mathlib.h"

#if PYTHON_API_VERSION <= 1012  /* Python 2.3 is 1012.  */
#define Py_ssize_t int
#endif

#define PI 3.1415926535897931

/* typedef double v4sf __attribute__ ((vector_size (32))); */
/* typedef float v4sf __attribute__ ((mode(V4SF))); */

/* 
typedef float v4sf __attribute__ ((vector_size (16)));
union f4vector 
{
    v4sf v;
    float f[4];
    };
*/


float mean(float *X, int length)
{
    int i;
    float sum;
    float sum0, sum1, sum2, sum3;

    sum = 0;
    sum0 = sum1 = sum2 = sum3 = 0;
    for(i=0; i<=length-4; i+=4) {
	sum0 += X[i];
	sum1 += X[i+1];
	sum2 += X[i+2];
	sum3 += X[i+3];
    }
    for(; i<length; i++) {
	sum += X[i];
    }
    sum += (sum0+sum1) + (sum2+sum3);
    return sum / length;
}

float var(float *X, int length)
{
    int i;
    float sum;
    float m, x;

    m = mean(X, length);
    sum = 0;
    for(i=0; i<length; i++) {
	/* Can optimize this loop. */
	x = X[i] - m;
	sum += x*x;
    }
    return sum / (length-1);
}

float *cov_byrow(float *X, int nrow, int ncol) 
{
    int i, j;
    int i1, i2;
    float sum;
    float *cov;
    float *means;

    if((X == NULL) || (nrow < 1 || ncol < 1))
	return NULL;

    cov = means = NULL;
    if(!(cov = (float *)malloc(nrow*nrow*sizeof(*cov))))
	return NULL;
    if(!(means = (float *)malloc(nrow*sizeof(*means)))) {
	free(cov); return NULL;
    }

    for(i=0; i<nrow; i++) {
	sum = 0;
	for(j=0; j<ncol; j++)
	    sum += X[i*ncol+j];
	means[i] = sum/ncol;
	/*printf("%d %f\n", i, means[i]);*/
    }

    for(i1=0; i1<nrow; i1++) {
	for(i2=i1; i2<nrow; i2++) {
	    /* Covariance = E((X-mean(X))(Y-mean(Y))) */
	    sum = 0;
	    /* Can optimize this loop. */
	    for(j=0; j<ncol; j++)
		sum += (X[i1*ncol+j]-means[i1]) * (X[i2*ncol+j]-means[i2]);
	    /* Results for (ncol-1) matches numpy and R */
	    cov[i1*nrow+i2] = cov[i2*nrow+i1] = sum/(ncol-1);
	}
    }

    free(means);
    return cov;
}

// Use the SSE Intrinsics.
// Generates slightly different results as the C version.
//float dot_intrin(float *X, float *Y, int n) {
//    int i;
//    float sum;
//    __m128 n1, n2, n3, n4;
//
//    sum = 0.0;
//    n4 = _mm_setzero_ps();  // set to 0
//    for(i=0; i<=n-4; i+=4) {
//	n1 = _mm_loadu_ps(X+i);
//	n2 = _mm_loadu_ps(Y+i);
//	n3 = _mm_mul_ps(n1, n2);
//	n3 = _mm_hadd_ps(n3, n3);
//	n4 = _mm_add_ps(n4, n3);
//    }
//    for(; i<n; i++)
//	sum += X[i] * Y[i];
//    n4 = _mm_hadd_ps(n4, n4);
//    _mm_store_ss(&sum, n4);
//
//    return sum;
//}

double dot_C_d(double *X, double *Y, int n) {
    int i;
    double sum;
    double sum0, sum1, sum2, sum3;
    /*union f4vector *a, *b, c;*/
    /*v4sf a, b, c;
      float temp[4];*/
    
    sum = 0;
    sum0 = sum1 = sum2 = sum3 = 0;
    /*c.f[0] = c.f[1] = c.f[2] = c.f[3] = 0;*/
    /*c = __builtin_ia32_xorps(c, c); // zero the accumulator */
    for(i=0; i<=n-4; i+=4) {
	/* Implementations in order of fastest to slowest. */
	sum0 += X[i]   * Y[i];
	sum1 += X[i+1] * Y[i+1];
	sum2 += X[i+2] * Y[i+2];
	sum3 += X[i+3] * Y[i+3];
	/*a = __builtin_ia32_loadups(X+i);
	  b = __builtin_ia32_loadups(Y+i);
	  c = __builtin_ia32_addps(c, __builtin_ia32_mulps(a, b));  */
	/* sum += X_i1[j] * X_i2[j]; */
	/* sum += *(X_i1+j) * *(X_i2+j); */
	/* sum += *X_i1++ * *X_i2++; */
	/* sum += X[i1*ncol+j] * X[i2*ncol+j]; */
	/* sum += (X[i1*ncol+j]-means[i1])*(X[i2*ncol+j]-means[i2]); */
	/* sum += (*X_i1++-mean_i1) * (*X_i2++-mean_i2); */
	/* Use SSE library.  Really slow. */
	/*a = (union f4vector *)(X_i1+j);
	  b = (union f4vector *)(X_i2+j);
	  c.v = c.v + a->v * b->v;*/
    }
    sum += (sum0+sum1) + (sum2+sum3);
    for(; i<n; i++)
	sum += X[i] * Y[i];
    /* sum0 = c.f[0]; sum1 = c.f[1]; sum2 = c.f[2]; sum3 = c.f[3]; */
    /*__builtin_ia32_storeups(temp, c);
      sum += temp[0] + temp[1] + temp[2] + temp[3]; */

    return sum;
}

double dot_d(double *X, double *Y, int n) {
    //return dot_intrin(X, Y, n);
    return dot_C_d(X, Y, n);
}

float dot_f(float *X, float *Y, int n) {
    int i;
    float sum;
    float sum0, sum1, sum2, sum3;

    sum = 0;
    sum0 = sum1 = sum2 = sum3 = 0;
    for(i=0; i<=n-4; i+=4) {
	sum0 += X[i]   * Y[i];
	sum1 += X[i+1] * Y[i+1];
	sum2 += X[i+2] * Y[i+2];
	sum3 += X[i+3] * Y[i+3];
    }
    sum += (sum0+sum1) + (sum2+sum3);
    for(; i<n; i++)
	sum += X[i] * Y[i];
    return sum;
}

#define align_ptr(data) (void *)(((int)data + 15) &~ 0x0F)

double *cor_byrow_d(double *X, int nrow, int ncol, int safe)
{
    int i, j;
    int i1, i2;
    double sum, mean, sd, x;
    double *cor;
    double *X_i, *X_i1, *X_i2;
    double EPS;
    //float *X_aligned;
    //float *stddevs;

    // dot_intrin is not helped by aligning the pointer.
    // BUG: Will leak memory.  Just for test.
    // X_aligned = (float *)malloc(nrow*ncol*sizeof(*X) + 2);
    // X_aligned = (float *)align_ptr(X_aligned);
    // memcpy(X_aligned, X, nrow*ncol*sizeof(*X));
    // X = X_aligned; 

    //stddevs = (float *)malloc(nrow*sizeof(*stddevs));



    /* Rounding error.  See below. */
    EPS = 1e-5;

    if(!X || nrow < 1 || ncol < 1)
	return NULL;

    if(!(cor = (double *)malloc(nrow*nrow*sizeof(*cor))))
	return NULL;

    /* If there's only 1 sample, there will be a divide-by-zero. */ 
    if(ncol == 1) {
	if(safe) {
	    /* Set correlations to default of 0. */
	    memset(cor, 0, nrow*nrow*sizeof(*cor));
	    return cor;
	}
	free(cor);
	return NULL;
    }

    /* Covariance = E((X-mean(X))(Y-mean(Y))) */
    /* Correlation = Cov(X, Y)/(SD(X) * SD(Y)) */

    /* Calculate the means and standard deviations of each row. */
    for(i=0; i<nrow; i++) {
	X_i = X + i*ncol;
	/* Calculate the mean of the row. */
    	sum = 0;
	for(j=0; j<ncol; j++)
	    sum += X_i[j];
	mean = sum / ncol;

	/* Calculate the standard deviation of the row. */
	sum = 0;
	for(j=0; j<ncol; j++)
	    sum += (X_i[j]-mean)*(X_i[j]-mean);
	sd = sqrt(sum / (ncol-1));

	if(!sd && !safe) {
	    /* If the variance is 0, then return an error. */
	    free(cor);
	    return NULL;
	}

	/* Normalize each row to N(0, 1) to make the the correlation
	   coefficient calculation easier. */
	for(j=0; j<ncol; j++) {
	    if(sd == 0) {
		X_i[j] = 0;
		continue;
	    } 
	    //X_i[j] -= mean;
	    X_i[j] = (X_i[j]-mean)/sd;
	}
	
	//stddevs[i] = sd;
    }

    /* Calculate the correlations. */
    for(i1=0; i1<nrow; i1++) {
	X_i1 = X+i1*ncol;
    	for(i2=i1; i2<nrow; i2++) {
	    X_i2 = X+i2*ncol;

	    sum = dot_d(X_i1, X_i2, ncol);
	    //	    sum = 0;
	    //	    sum0 = sum1 = sum2 = sum3 = 0;
	    //	    /*c.f[0] = c.f[1] = c.f[2] = c.f[3] = 0;*/
	    //	    for(j=0; j<=ncol-4; j+=4) {
	    //		/* Implementations in order of fastest to slowest. */
	    //		sum0 += X_i1[j]   * X_i2[j];
	    //		sum1 += X_i1[j+1] * X_i2[j+1];
	    //		sum2 += X_i1[j+2] * X_i2[j+2];
	    //		sum3 += X_i1[j+3] * X_i2[j+3];
	    //		/*a = __builtin_ia32_loadups(X_i1+j);
	    //		b = __builtin_ia32_loadups(X_i2+j);
	    //		c = __builtin_ia32_addps(c, __builtin_ia32_mulps(a, b)); */
	    //		/* sum += X_i1[j] * X_i2[j]; */
	    //		/* sum += *(X_i1+j) * *(X_i2+j); */
	    //		/* sum += *X_i1++ * *X_i2++; */
	    //		/* sum += X[i1*ncol+j] * X[i2*ncol+j]; */
	    //		/* sum += (X[i1*ncol+j]-means[i1])*(X[i2*ncol+j]-means[i2]); */
	    //		/* sum += (*X_i1++-mean_i1) * (*X_i2++-mean_i2); */
	    //		/* Use SSE library.  Really slow. */
	    //		/*a = (union f4vector *)(X_i1+j);
	    //		b = (union f4vector *)(X_i2+j);
	    //		c.v = c.v + a->v * b->v;*/
	    //	    }
	    //	    for(; j<ncol; j++) {
	    //		sum += X_i1[j] * X_i2[j];
	    //	    }
	    //	    /* sum0 = c.f[0]; sum1 = c.f[1]; sum2 = c.f[2]; sum3 = c.f[3]; */
	    //	    sum += (sum0+sum1) + (sum2+sum3);
	    //	    /*__builtin_ia32_storeups(tmp, c);*/

	    /* numpy and R divide by (ncol-1). */
	    x = sum/(ncol-1);
	    /* Calculate the correlation.  Could optimize this by
	       putting the division in the X matrix. */

	    //if(!stddevs[i1] || !stddevs[i2])
	    //x = 0;
	    //else
	    //x = x/(stddevs[i1]*stddevs[i2]);

	    /* Some values might be greater than 1 or -1 due to
	       rounding error.  Fix this. */
	    if(x > 1.0 && x < 1.0+EPS)
		x = 1.0;
	    else if(x < -1.0 && x > -1.0-EPS)
		x = -1.0;
	    cor[i1*nrow+i2] = cor[i2*nrow+i1] = x;
	}
    }

    return cor;
}

float *cor_byrow_f(float *X, int nrow, int ncol, int safe)
{
    int i, j;
    int i1, i2;
    float sum, mean, sd, x;
    float *cor;
    float *X_i, *X_i1, *X_i2;
    float EPS;

    /* Rounding error.  See below. */
    EPS = 1e-5;

    if(!X || nrow < 1 || ncol < 1)
	return NULL;

    if(!(cor = (float *)malloc(nrow*nrow*sizeof(*cor))))
	return NULL;

    /* If there's only 1 sample, there will be a divide-by-zero. */ 
    if(ncol == 1) {
	if(safe) {
	    /* Set correlations to default of 0. */
	    memset(cor, 0, nrow*nrow*sizeof(*cor));
	    return cor;
	}
	free(cor);
	return NULL;
    }

    /* Covariance = E((X-mean(X))(Y-mean(Y))) */
    /* Correlation = Cov(X, Y)/(SD(X) * SD(Y)) */

    /* Calculate the means and standard deviations of each row. */
    for(i=0; i<nrow; i++) {
	X_i = X + i*ncol;
	/* Calculate the mean of the row. */
    	sum = 0;
	for(j=0; j<ncol; j++)
	    sum += X_i[j];
	mean = sum / ncol;

	/* Calculate the standard deviation of the row. */
	sum = 0;
	for(j=0; j<ncol; j++)
	    sum += (X_i[j]-mean)*(X_i[j]-mean);
	sd = sqrt(sum / (ncol-1));

	if(!sd && !safe) {
	    /* If the variance is 0, then return an error. */
	    free(cor);
	    return NULL;
	}

	/* Normalize each row to N(0, 1) to make the the correlation
	   coefficient calculation easier. */
	for(j=0; j<ncol; j++) {
	    if(sd == 0) {
		X_i[j] = 0;
		continue;
	    } 
	    X_i[j] = (X_i[j]-mean)/sd;
	}
    }

    /* Calculate the correlations. */
    for(i1=0; i1<nrow; i1++) {
	X_i1 = X+i1*ncol;
    	for(i2=i1; i2<nrow; i2++) {
	    X_i2 = X+i2*ncol;

	    sum = dot_f(X_i1, X_i2, ncol);

	    /* numpy and R divide by (ncol-1). */
	    x = sum/(ncol-1);
	    /* Calculate the correlation.  Could optimize this by
	       putting the division in the X matrix. */

	    /* Some values might be greater than 1 or -1 due to
	       rounding error.  Fix this. */
	    if(x > 1.0 && x < 1.0+EPS)
		x = 1.0;
	    else if(x < -1.0 && x > -1.0-EPS)
		x = -1.0;
	    cor[i1*nrow+i2] = cor[i2*nrow+i1] = x;
	}
    }

    return cor;
}

double fisher_z(double R, int N)
{
    double W, var_W;
    double EPS;

    /* R cannot be exactly 1 or -1. */
    EPS = 1E-10;
    R = (R <= 1.0-EPS) ? R : 1.0-EPS;
    R = (R >= -1.0+EPS) ? R : -1.0+EPS;

    if(N <= 25)
	var_W = 1.0 / (N-1);
    else 
	var_W = 1.0 / (N-3);

    W = 0.5 * log((1+R)/(1-R));
    if(N <= 25)
	W = W - (3*W+tanh(W)) / (4*(N-1));
    return W / sqrt(var_W);
}

/* Calculate density of N(0, 1) distribution. */
double dnorm(double x, double mean, double sd)
{
    double a, b, c;
    a = sd * sqrt(2*PI);
    b = x - mean;
    c = -(b*b)/(2.0*sd*sd);
    return(1.0/a*exp(c));
}

float dparzen(float x, float *X_obs, int *X_count, int length, float h)
{
    int i;
    int N;
    float total;
    float y, z;

    N = 0;
    total = 0.0;
    for(i=0; i<length; i++) {
	N += X_count[i];
	y = (x-X_obs[i])/h;
	z = -(y*y)/2.0;
	total += (float)X_count[i] * exp(z);
	/*total += dnorm((x-X_obs[i])/h) * X_count[i];*/
    }
    total = total / sqrt(2*PI);
    return total/(N*h);
}

float dparzen_fast(float x, float *X_obs, int *X_count, int length, 
  float h, int N, int max_count, float EPS)
{
    /* N is the sum of X_count.  max_count is max(X_count).  EPS is
       the error we are willing to tolerate for this estimate. Cached
       for faster access.  Also reduces the number of calls to exp. */
    int i;
    float total;
    float y;

    /* y = -[(x-X_obs[i])/h]^2 / 2
     * EPS = [X_count[i] * exp(y)] / [N*h*sqrt(2*PI)]
     * EPS*N*h*sqrt(2*PI) = X_count[i] * exp(y)
     * [EPS*N*h*sqrt(2*PI)] / max_count <= exp(y)
     * Expanding further (e.g. the y) does not make things much faster.
     */

    EPS = log((EPS*N*h*sqrt(2*PI))/max_count);
    total = 0.0;
    for(i=0; i<length; i++) {
	y = (x-X_obs[i])/h;
	y = -(y*y)/2.0;
	if(y < EPS) continue;
	total += (float)X_count[i] * exp(y);
	/*total += dnorm((x-X_obs[i])/h) * X_count[i];*/
    }
    total = total / sqrt(2*PI);
    return total/(N*h);
}

PyObject *c2py_ivector(int *X, int length)
{
    int i;
    PyObject *py_X, *py_item;

    if(!(py_X = PyList_New(length)))
	return NULL;
    for(i=0; i<length; i++) {
	if(!(py_item = PyInt_FromLong((long)X[i]))) {
	    Py_DECREF(py_X); return NULL;
	}
	/* SET_ITEM steals reference to py_item--no DECREF necessary. */
	PyList_SET_ITEM(py_X, i, py_item);
    }

    return py_X;
}

/* Convert Python object to a vector. */
int py2c_ivector(PyObject *py_X, int **X, int *length)
{
    PyObject *py_item;
    int i;

    py_item = NULL;
    *X = NULL;
    if((*length = PySequence_Size(py_X)) == -1)
	return 0;
    if(!(*X = (int *)malloc((*length) * sizeof(**X)))) {
	PyErr_SetString(PyExc_MemoryError, "out of memory");
	return 0;
    }

    for(i=0; i<*length; i++) {
	if(!(py_item = PySequence_GetItem(py_X, i))) {  /* new ref */
	    goto py2c_ivector_err;
	}
	(*X)[i] = int_FromPyObject(py_item);
	if(PyErr_Occurred()) {
	    goto py2c_ivector_err;
	}
	Py_DECREF(py_item); py_item = NULL;
    }

    return 1;
  py2c_ivector_err:
    if(py_item) { Py_DECREF(py_item); }
    if(*X) { free(*X); *X = NULL; }
    return 0;
}

int py2c_fvector(PyObject *py_X, float **X, int *length)
{
    PyObject *py_item;
    int i;

    py_item = NULL;
    *X = NULL;
    if((*length = PySequence_Size(py_X)) == -1)
	return 0;
    if(!(*X = (float *)malloc((*length) * sizeof(**X)))) {
	PyErr_SetString(PyExc_MemoryError, "out of memory");
	return 0;
    }

    for(i=0; i<*length; i++) {
	if(!(py_item = PySequence_GetItem(py_X, i))) {  /* new ref */
	    goto py2c_vector_err;
	}
	(*X)[i] = (float)double_FromPyObject(py_item);
	if(PyErr_Occurred()) {
	    goto py2c_vector_err;
	}
	Py_DECREF(py_item); py_item = NULL;
    }

    return 1;
  py2c_vector_err:
    if(py_item) { Py_DECREF(py_item); }
    if(*X) { free(*X); *X = NULL; }
    return 0;
}

/* Convert a matrix to a Python object. */
PyObject *c2py_fmatrix(float *X, int nrow, int ncol)
{
    int i, j;
    PyObject *py_X, *py_X_row, *py_item;
    float *X_i;

    if(!(py_X = PyList_New(nrow)))
	return NULL;
    for(i=0; i<nrow; i++) {
	if(!(py_X_row = PyList_New(ncol))) {
	    Py_DECREF(py_X); return NULL;
	}
	X_i = X + i*ncol;
	for(j=0; j<ncol; j++) {
	    if(!(py_item = PyFloat_FromDouble((double)X_i[j]))) {
		/* Bug: The items in py_X will not be DECREF'd. */
		Py_DECREF(py_X); return NULL;
	    }
	    /* SET_ITEM steals reference to py_item--no DECREF necessary. */
	    PyList_SET_ITEM(py_X_row, j, py_item);
	}
	/* SET_ITEM steals reference to py_X_row--no DECREF necessary. */
	PyList_SET_ITEM(py_X, i, py_X_row);
    }

    return py_X;
}

PyObject *c2py_dmatrix(double *X, int nrow, int ncol)
{
    int i, j;
    PyObject *py_X, *py_X_row, *py_item;
    double *X_i;

    if(!(py_X = PyList_New(nrow)))
	return NULL;
    for(i=0; i<nrow; i++) {
	if(!(py_X_row = PyList_New(ncol))) {
	    Py_DECREF(py_X); return NULL;
	}
	X_i = X + i*ncol;
	for(j=0; j<ncol; j++) {
	    if(!(py_item = PyFloat_FromDouble(X_i[j]))) {
		/* Bug: The items in py_X will not be DECREF'd. */
		Py_DECREF(py_X); return NULL;
	    }
	    /* SET_ITEM steals reference to py_item--no DECREF necessary. */
	    PyList_SET_ITEM(py_X_row, j, py_item);
	}
	/* SET_ITEM steals reference to py_X_row--no DECREF necessary. */
	PyList_SET_ITEM(py_X, i, py_X_row);
    }

    return py_X;
}

/* Convert Python object to a matrix. */
int py2c_fmatrix(PyObject *py_X, float **X, int *nrow, int *ncol)
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

    if(!(*X = (float *)malloc(*nrow * *ncol * sizeof(**X)))) {
	PyErr_SetString(PyExc_MemoryError, "out of memory");
	goto py2c_matrix_err;
    }

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
	    (*X)[index] = (float)double_FromPyObject(py_item);
	    index++;
	}
	py_item = NULL;
	Py_DECREF(py_Fast_X_row); py_Fast_X_row = NULL;
	if(PyErr_Occurred()) {
	    /* Optimization: This might be set in double_FromPyObject,
	       but check here instead of inside loop. */
	    goto py2c_matrix_err;
	}
    }

    return 1;

  py2c_matrix_err:
    if(py_item) { Py_DECREF(py_item); }
    if(py_X_row) { Py_DECREF(py_X_row); }
    if(py_Fast_X_row) { Py_DECREF(py_Fast_X_row); }
    if(*X) { free(*X); *X = NULL; }
    *nrow = *ncol = 0;
    return 0;
}

int py2c_dmatrix(PyObject *py_X, double **X, int *nrow, int *ncol)
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

    if(!(*X = (double *)malloc(*nrow * *ncol * sizeof(**X)))) {
	PyErr_SetString(PyExc_MemoryError, "out of memory");
	goto py2c_matrix_err;
    }

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
	    (*X)[index] = double_FromPyObject(py_item);
	    index++;
	}
	py_item = NULL;
	Py_DECREF(py_Fast_X_row); py_Fast_X_row = NULL;
	if(PyErr_Occurred()) {
	    /* Optimization: This might be set in double_FromPyObject,
	       but check here instead of inside loop. */
	    goto py2c_matrix_err;
	}
    }

    return 1;

  py2c_matrix_err:
    if(py_item) { Py_DECREF(py_item); }
    if(py_X_row) { Py_DECREF(py_X_row); }
    if(py_Fast_X_row) { Py_DECREF(py_Fast_X_row); }
    if(*X) { free(*X); *X = NULL; }
    *nrow = *ncol = 0;
    return 0;
}

double double_FromPyObject(PyObject *py_obj)
{
    PyObject *py_float;
    double c_float;

    if(!(py_float = PyNumber_Float(py_obj))) {
	return 0;
    }
    c_float = PyFloat_AsDouble(py_float);
    Py_DECREF(py_float);
    return c_float;
}

int int_FromPyObject(PyObject *py_obj)
{
    PyObject *py_int;
    int c_int;

    if(!(py_int = PyNumber_Int(py_obj))) {
	return 0;
    }
    c_int = (int)PyInt_AsLong(py_int);
    Py_DECREF(py_int);
    return c_int;
}

