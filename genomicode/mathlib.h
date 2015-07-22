float mean(float *X, int length);
float var(float *X, int length);
float *cov_byrow(float *X, int nrow, int ncol);
double *cor_byrow_d(double *X, int nrow, int ncol, int safe);
float *cor_byrow_f(float *X, int nrow, int ncol, int safe);
float *cor_bycol_f(float *X, int nrow, int ncol, int safe);

double fisher_z(double R, int N);

double dnorm(double x, double mean, double sd);
float dparzen(float x, float *X_obs, int *X_count, int length, float h);
float dparzen_fast(float x, float *X_obs, int *X_count, int length, 
		   float h, int N, int max_count, float EPS);

PyObject *c2py_ivector(int *X, int length);
PyObject *c2py_fmatrix(float *X, int nrow, int ncol);
PyObject *c2py_dmatrix(double *X, int nrow, int ncol);
int py2c_ivector(PyObject *py_X, int **X, int *length);
int py2c_fvector(PyObject *py_X, float **X, int *length);
int py2c_fmatrix(PyObject *py_X, float **X, int *nrow, int *ncol);
int py2c_dmatrix(PyObject *py_X, double **X, int *nrow, int *ncol);

double double_FromPyObject(PyObject *py_obj);
int int_FromPyObject(PyObject *py_obj);
