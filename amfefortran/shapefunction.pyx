# cython: language_level=3
import numpy as np
cimport numpy as np

cdef extern:
    int c_line2_no_of_nodes()
    int c_line2_no_of_local_coordinates()
    double c_line2_xi_lower()
    double c_line2_xi_upper()
    void c_line2_evaluate(double *xi, double *out)


class Line2:
    def __init__(self):
        return

    def no_of_nodes(self):
        return c_line2_no_of_nodes()

    def no_of_local_coordinates(self):
        return c_line2_no_of_local_coordinates()

    def evaluate(self, double xi):
        cdef np.ndarray[double, ndim=1, mode="fortran"] out
        out = np.zeros(2, np.double, order='F')
        c_line2_evaluate(&xi, &out[0])
        return out
