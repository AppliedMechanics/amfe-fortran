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

        function c_line2_no_of_local_coordinates() bind(c, name='c_line2_no_of_local_coordinates')
            integer(c_int) :: c_line2_no_of_local_coordinates
            c_line2_no_of_local_coordinates = line2_no_of_local_coordinates()
        end function c_line2_no_of_local_coordinates

        function c_line2_xi_lower() bind(c, name='c_line2_xi_lower')
            real(c_double) :: c_line2_xi_lower
            c_line2_xi_lower = line2_xi_lower()
        end function c_line2_xi_lower

        function c_line2_xi_upper() bind(c, name='c_line2_xi_upper')
            real(c_double) :: c_line2_xi_upper
            c_line2_xi_upper = line2_xi_upper()
        end function c_line2_xi_upper

        subroutine c_line2_evaluate(xi, out) bind(c, name='c_line2_evaluate')
            real(c_double), intent(in) :: xi
            real(c_double), intent(out) :: out(2)
            call line2_evaluate(xi, out)
        end subroutine c_line2_evaluate

        subroutine c_line2_jacobian(xi, out) bind(c, name='c_line2_jacobian')
            real(c_double), intent(in) :: xi
            real(c_double), intent(out) :: out(2,1)
            call line2_jacobian(xi, out)
        end subroutine c_line2_jacobian

end module shapefunction_wrapper