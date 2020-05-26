cdef extern:
    int c_line2_no_of_nodes()

def line2_no_of_nodes():
    return c_line2_no_of_nodes()
