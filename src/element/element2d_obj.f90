! Copyright (c) 2020, Lehrstuhl fuer Angewandte Mechanik, Technische
! Universitaet Muenchen.
!
! Distributed under BSD-3-Clause License. See LICENSE-File for more information
!
! Author: Christian Meyer
!

module element2d_obj
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    use shapefunction_obj
    use integrand_structural2d_obj
    use material_obj
    implicit none

    private
    public element2d

    ! A PDT implementation was not successful due to following bug
    ! CF. BUG https://gcc.gnu.org/bugzilla/show_bug.cgi?id=82943
    ! Thus, an implementation with allocatables is used
    type element2d
        ! fixed types:
        integer :: no_of_dofs
        type(shapefunction2d) :: shapefunction
        ! exact type unknown before only abstract type known, and thus must be allocatable:
        class(integrand2d), allocatable :: integrand
        class(material), allocatable :: material
        real(dp), allocatable :: gausspoints_kf_points(:,:)
        real(dp), allocatable :: gausspoints_kf_weights(:)
        real(dp), allocatable :: gausspoints_m_points(:,:)
        real(dp), allocatable :: gausspoints_m_weights(:)

        contains
            procedure :: f_int_and_jac => element2d_f_int_and_jac

    end type element2d

    contains

        subroutine element2d_f_int_and_jac(this, out, xn, d, dd, t, s0, t0, s1)
            implicit none
            class(element2d) :: this ! type bound procedure (not object bounded), thus 'class' instead of 'type'
            real(dp), intent(inout) :: out(this%no_of_dofs, this%no_of_dofs*2+1) ! size no_of_dofs, no_of_dofs*2+1
            real(dp), intent(in) :: xn(this%shapefunction%no_of_nodes * 2) ! size 2*no_of_nodes
            real(dp), intent(in) :: d(this%no_of_dofs)
            real(dp), intent(in) :: dd(this%no_of_dofs)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:) !s0(this%material%no_of_states() * this%gausspoints_kf%no_of_gps)
            real(dp), intent(in) :: t0
            real(dp), intent(inout) :: s1(:) !s1(this%material%no_of_states() * this%gausspoints_kf%no_of_gps)

            integer :: current_gp, no_of_gps, no_of_states_per_gp
            integer :: g

            real(dp) :: n(this%shapefunction%no_of_nodes)
            real(dp) :: dn_dxi(this%shapefunction%no_of_nodes,2)
            real(dp) :: f_int_k_dk_gp(this%no_of_dofs, this%no_of_dofs*2+1)

            ! allocatable variables
            real(dp), allocatable :: s0_local(:)
            real(dp), allocatable :: s1_local(:)

            ! allocate:
            allocate(s0_local(this%material%no_of_states()))
            allocate(s1_local(this%material%no_of_states()))

            ! initialize
            current_gp = 1
            s1(:) = 0.0_dp
            out(:,:) = 0.0_dp
            f_int_k_dk_gp(:,:) = 0.0_dp

            no_of_gps = size(this%gausspoints_kf_points, 2)
            no_of_states_per_gp = this%material%no_of_states()

            ! run loop over gauss points
            do g=1,no_of_gps
                ! ---------------------------- Initializations
                s0_local = s0((current_gp-1)*no_of_states_per_gp+1:current_gp*no_of_states_per_gp+1)
                s1_local(:) = 0.0_dp
                f_int_k_dk_gp(:,:) = 0.0_dp


                ! ---------------------------- Call shapefunctions --------------------------------------------
                call this%shapefunction%evaluate(this%gausspoints_kf_points(1,g), this%gausspoints_kf_points(2,g), n)


                ! ---------------------------- Call weak form -------------------------------------------------
                call this%integrand%f_int_and_jac(n, dn_dxi, xn, d, dd, t, matfunc2d_internal, s0_local, t0, s1_local, &
                f_int_k_dk_gp, this%integrand%no_of_fields(), 2, this%shapefunction%no_of_nodes, no_of_states_per_gp)

                ! ---------------------------- Update s1 -------------------------------------------------------
                s1(current_gp-1*no_of_states_per_gp+1:current_gp*no_of_states_per_gp+1) = s1_local
                out = out + this%gausspoints_kf_weights(g) * f_int_k_dk_gp
                ! ---------------------------- Increase current_gp ---------------------------------------------
                current_gp = current_gp + 1
            end do

            contains
                subroutine matfunc2d_internal(f, df, t, s0, t0, s1, sv_out, c_out, dc_out)
                real(dp), intent(in) :: f(2,2)
                real(dp), intent(in) :: df(2,2)
                real(dp), intent(in) :: t
                real(dp), intent(in) :: s0(:)
                real(dp), intent(in) :: t0
                real(dp), intent(out) :: s1(:)
                real(dp), intent(out) :: sv_out(3)
                real(dp), intent(out) :: c_out(3,3)
                real(dp), intent(out) :: dc_out(3,3)

                call this%material%sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out)

                end subroutine matfunc2d_internal

        end subroutine element2d_f_int_and_jac


end module element2d_obj
