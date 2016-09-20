#define NPY_NO_DEPRECATED_API NPY_1_6_API_VERSION
#include </usr/include/python2.7/Python.h>
#include </usr/local/lib/python2.7/dist-packages/numpy/core/include/numpy/arrayobject.h>
#include <stdio.h>
#include <complex.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>

void t0(double complex *x,int64_t NN);

static PyObject *t0_think(PyObject *self, PyObject *args);

static char module_docstring[] =
  "doc stuff ";
static char t0_docstring[] =
  "Rearrange elements so each row's elements share the same irreducibles of whom who share the same non-zero elements."; 

static PyMethodDef module_methods[] = {
  {"t0", t0_think, METH_VARARGS, t0_docstring},
  {NULL, NULL, 0, NULL} // tells the compiler there are no more methods
};

PyMODINIT_FUNC init_think(void)
{
  PyObject *mm = Py_InitModule3("_think", module_methods, module_docstring);
  if (mm == NULL)
    { // i think this checks if the initialization of our function fgt failed
      return; // exit python function call?
    }
  /* Load `numpy` functionality. */
  import_array(); // has an implicit declaration warning, probabliy resolved when python compiles
}

static PyObject *t0_think(PyObject *self, PyObject *args)
{
  PyObject *x_obj,*x_array,*y_array;
  int64_t *x,*y;
  int64_t N0,N1,N2;

  if (!PyArg_ParseTuple(args, "O", &x_obj))
    {
      printf("_fgt error : only one input allowed\n");
      return 1;
    }
  x_array = PyArray_FROM_OTF(x_obj, NPY_UBYTE, NPY_IN_ARRAY);

  // int64_t *x_dims = PyArray_DIMS(x_array);

  x = (uint8_t *)PyArray_DATA(x_array);
  N0 = (int64_t)PyArray_DIM(*x_array, 0);
  N1 = (int64_t)PyArray_DIM(*x_array, 1);
  N2 = (int64_t)PyArray_DIM(*x_array, 2);
  // npy_intp dims[2];
  // dims[0] = NN;
  // dims[1] = NN;
  // y_array = PyArray_SimpleNew(2, dims, NPY_INT64);
  // y = (double complex *)PyArray_DATA(*y_array);

  //  return y_array;
  return;
}

void t0(uint8_t * memory, uint8_t* raw_thought, int64_t m0, int64_t m1, int64_t m2, int64_t rt0)
  {

  }
