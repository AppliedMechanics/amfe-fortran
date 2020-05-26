module shapefunction_wrapper

    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use shapefunction_fun

    implicit none

    private
    public c_line2_no_of_nodes

    contains

        function c_line2_no_of_nodes() bind(c, name='c_line2_no_of_nodes')
            integer(c_int) :: c_line2_no_of_nodes
            c_line2_no_of_nodes = line2_no_of_nodes()
        end function c_line2_no_of_nodes

end module shapefunction_wrapper