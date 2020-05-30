!
! Copyright (c) 2020 TECHNICAL UNIVERSITY OF MUNICH, DEPARTMENT OF MECHANICAL ENGINEERING, CHAIR OF APPLIED MECHANICS,
! BOLTZMANNSTRASSE 15, 85748 GARCHING/MUNICH, GERMANY, RIXEN@TUM.DE.
!
! Distributed under 3-Clause BSD license. See LICENSE file for more information.
!

module integrand_structural2d_obj
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    use integrand_structural2d_nonlinear_fun

    implicit none
    private

    public integrand2d, integrand_structural2d_nonlinear

    type, abstract :: integrand2d
        contains
            procedure(fields_interface), nopass, deferred :: fields
            procedure(no_of_fields_interface), nopass, deferred :: no_of_fields
            procedure(no_of_dofs_interface), nopass, deferred :: no_of_dofs
            procedure(spatial_dimension_interface), nopass, deferred :: spatial_dimension
            procedure(f_int_and_jac_interface), nopass, deferred :: f_int_and_jac
            procedure(m_int_interface), nopass, deferred :: m_int
    end type integrand2d

    abstract interface
        integer function no_of_fields_interface()
        end function no_of_fields_interface

        subroutine fields_interface(fields)
            integer, intent(out) :: fields(:)
        end subroutine fields_interface

        integer function no_of_dofs_interface(no_of_nodes)
            integer, intent(in) :: no_of_nodes
        end function no_of_dofs_interface

        integer function spatial_dimension_interface()
        end function spatial_dimension_interface

        subroutine f_int_and_jac_interface(n, dn_dxi, xn, d, dd, t, matfunc_2d, s0, t0, s1, &
                f_int_k_dk_out, no_of_fields, spatial_dimension, no_of_nodes, no_of_states)
            import dp
            integer, intent(in) :: spatial_dimension
            integer, intent(in) :: no_of_nodes
            integer, intent(in) :: no_of_states
            integer, intent(in) :: no_of_fields
            real(dp), intent(in) :: n(no_of_nodes) ! size(no_of_nodes)
            real(dp), intent(in) :: dn_dxi(no_of_nodes,spatial_dimension) ! size(no_of_nodes, spatial_dimension)
            real(dp), intent(in) :: xn(no_of_nodes*spatial_dimension) ! size(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: d(no_of_nodes*spatial_dimension) ! size(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: dd(no_of_nodes*spatial_dimension) ! size(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: t
            interface
                subroutine sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out)
                    import :: dp
                    real(dp), intent(in) :: f(2,2)
                    real(dp), intent(in) :: df(2,2)
                    real(dp), intent(in) :: t
                    real(dp), intent(in) :: s0(:)
                    real(dp), intent(in) :: t0
                    real(dp), intent(out) :: s1(:)
                    real(dp), intent(out) :: sv_out(3)
                    real(dp), intent(out) :: c_out(3,3)
                    real(dp), intent(out) :: dc_out(3,3)
                end subroutine sv_c_dc_2d
            end interface
            procedure(sv_c_dc_2d) :: matfunc_2d
            real(dp), intent(in) :: s0(no_of_states) ! size(no_of_states)
            real(dp), intent(in) :: t0
            real(dp), intent(inout) :: s1(no_of_states) ! size(no_of_states)
            real(dp), intent(inout) :: f_int_k_dk_out(no_of_nodes*spatial_dimension,no_of_nodes*spatial_dimension*2+1)
        end subroutine f_int_and_jac_interface

        subroutine m_int_interface(n, dn_dxi, xn, d, dd, t, density_times_thickness, m_out, no_of_fields, &
        spatial_dimension, no_of_nodes)
            import dp
            integer, intent(in) :: no_of_fields
            integer, intent(in) :: spatial_dimension
            integer, intent(in) :: no_of_nodes
            real(dp), intent(in) :: n(no_of_nodes)
            real(dp), intent(in) :: dn_dxi(no_of_nodes,spatial_dimension)
            real(dp), intent(in) :: xn(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: d(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: dd(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: density_times_thickness
            real(dp), intent(inout) :: m_out(no_of_fields*no_of_nodes,no_of_fields*no_of_nodes)
        end subroutine m_int_interface

    end interface

    type, extends(integrand2d) :: integrand_structural2d_nonlinear
        contains
            procedure, nopass :: fields => igstrc2dnl_fields
            procedure, nopass :: no_of_fields => igstrc2dnl_no_of_fields
            procedure, nopass :: no_of_dofs => igstrc2dnl_no_of_dofs
            procedure, nopass :: spatial_dimension => igstrc2dnl_spatial_dimension
            procedure, nopass :: f_int_and_jac => igstrc2dnl_f_int_and_jac
            procedure, nopass :: m_int => igstrc2dnl_m_int
    end type integrand_structural2d_nonlinear

end module integrand_structural2d_obj