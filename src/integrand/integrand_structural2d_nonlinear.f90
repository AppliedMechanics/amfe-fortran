!
! Copyright (c) 2020 TECHNICAL UNIVERSITY OF MUNICH, DEPARTMENT OF MECHANICAL ENGINEERING, CHAIR OF APPLIED MECHANICS,
! BOLTZMANNSTRASSE 15, 85748 GARCHING/MUNICH, GERMANY, RIXEN@TUM.DE.
!
! Distributed under 3-Clause BSD license. See LICENSE file for more information.
!

module integrand_structural2d
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    use tools
    use fields

    implicit none
    private

    public igstrc2dnl_fields, igstrc2dnl_no_of_dofs, igstrc2dnl_spatial_dimension, igstrc2dnl_f_int_and_jac, &
            igstrc2dnl_m_int

    contains

        subroutine igstrc2dnl_fields(fields)
            integer, intent(out) :: fields(2)
            fields = (/ UX, UY /)
        end subroutine igstrc2dnl_fields

        integer function igstrc2dnl_no_of_dofs(no_of_nodes)
            integer, intent(in) :: no_of_nodes
            igstrc2dnl_no_of_dofs = 2*no_of_nodes
        end function igstrc2dnl_no_of_dofs

        integer function igstrc2dnl_spatial_dimension()
            igstrc2dnl_spatial_dimension = 2
        end function igstrc2dnl_spatial_dimension

        subroutine igstrc2dnl_f_int_and_jac(n, dn_dxi, xn, d, dd, t, matfunc_2d, s0, t0, s1, &
                f_int_k_dk_out, no_of_fields, spatial_dimension, no_of_nodes, no_of_states)
            integer, intent(in) :: spatial_dimension
            integer, intent(in) :: no_of_nodes
            integer, intent(in) :: no_of_states
            integer, intent(in) :: no_of_fields
            real(dp), intent(in) :: n(no_of_nodes)
            real(dp), intent(in) :: dn_dxi(no_of_nodes,spatial_dimension)
            real(dp), intent(in) :: xn(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: d(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: dd(no_of_nodes*spatial_dimension)
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
            real(dp), intent(in) :: s0(no_of_states)
            real(dp), intent(in) :: t0
            real(dp), intent(inout) :: s1(no_of_states)
            real(dp), intent(inout) :: f_int_k_dk_out(no_of_nodes*spatial_dimension,no_of_nodes*spatial_dimension*2+1)

            integer :: no_of_dofs
            real(dp) :: dn_dx(no_of_nodes,spatial_dimension)
            real(dp) :: du_dx(no_of_fields,spatial_dimension)
            real(dp) :: ddu_dx(no_of_fields,spatial_dimension)
            real(dp) :: detj
            real(dp) :: fval(2,2)
            real(dp) :: dfval(2,2)
            real(dp) :: bval(3,spatial_dimension*no_of_nodes)
            real(dp) :: k_geo(no_of_nodes*spatial_dimension,no_of_nodes*spatial_dimension)
            real(dp) :: s_v(3)
            real(dp) :: c(3,3)
            real(dp) :: dc(3,3)

            no_of_dofs = no_of_nodes*spatial_dimension

            call f_material_gradients_structural_2d(n, dn_dxi, xn, d, dd, t, du_dx, ddu_dx, dn_dx, detj, &
            no_of_nodes, spatial_dimension, no_of_fields)

            call f_deformation_gradient_2d(du_dx, fval)
            call f_deformation_gradient_2d(ddu_dx, dfval)

            call matfunc_2d(fval,dfval,t,s0,t0,s1,s_v,c,dc)

            call f_b_l_2d(fval, dn_dx, bval, no_of_nodes, spatial_dimension)

            call f_f_int_and_jac_structural_2d(bval, s_v, c, dc, detj, dn_dx, f_int_k_dk_out, k_geo, no_of_nodes, &
                spatial_dimension, no_of_dofs)

        end subroutine igstrc2dnl_f_int_and_jac

        subroutine igstrc2dnl_m_int(n, dn_dxi, xn, d, dd, t, density_times_thickness, m_out, no_of_fields, &
        spatial_dimension, no_of_nodes)

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

            real(dp) :: x_mat(spatial_dimension, no_of_nodes)
            real(dp) :: j(spatial_dimension, spatial_dimension)
            real(dp) :: detj
            real(dp) :: m_small(no_of_nodes, no_of_nodes)
            integer :: k, m

            x_mat = reshape(xn, (/ spatial_dimension, no_of_nodes /))
            j = matmul(x_mat, dn_dxi)

            detj = determinant_2d(j)

            do m=1,no_of_nodes
                do k=1,no_of_nodes
                    m_small(k,m) = n(k)*n(m)
                end do
            end do
            m_small = m_small * detJ * density_times_thickness

            call scatter_matrix(m_small, m_out, spatial_dimension, no_of_nodes, no_of_nodes)

        end subroutine igstrc2dnl_m_int

        subroutine f_material_gradients_structural_2d(n, dn_dxi, xn, d, dd, t, du_dx, ddu_dx, dn_dx, detj, &
        no_of_nodes, spatial_dimension, no_of_fields)
            integer, intent(in) :: no_of_nodes
            integer, intent(in) :: spatial_dimension
            integer, intent(in) :: no_of_fields
            real(dp), intent(in) :: n(no_of_nodes)
            real(dp), intent(in) :: dn_dxi(no_of_nodes, spatial_dimension)
            real(dp), intent(in) :: xn(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: d(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: dd(no_of_nodes*spatial_dimension)
            real(dp), intent(in) :: t
            real(dp), intent(out) :: dn_dx(no_of_nodes, spatial_dimension)
            real(dp), intent(out) :: du_dx(no_of_fields, spatial_dimension)
            real(dp), intent(out) :: ddu_dx(no_of_fields, spatial_dimension)
            real(dp), intent(out) :: detj

            real(dp) :: jinv(spatial_dimension, spatial_dimension)
            real(dp) :: d_mat(no_of_fields, no_of_nodes)
            real(dp) :: dd_mat(no_of_fields, no_of_nodes)
            real(dp) :: x_mat(spatial_dimension, no_of_nodes)
            real(dp) :: j(spatial_dimension, spatial_dimension)


            x_mat = reshape(xn, (/ spatial_dimension, no_of_nodes /))
            j = matmul(x_mat, dn_dxi)
            detj = determinant_2d(j)
            call f_inverse_matrix_2d(j, detj, jinv)

            dn_dX = matmul(dn_dxi, jinv)

            d_mat = reshape(d, (/ no_of_fields, no_of_nodes /))
            dd_mat = reshape(dd, (/ no_of_fields, no_of_nodes /))

            du_dx = matmul(d_mat, dn_dx)
            ddu_dx = matmul(dd_mat, dn_dx)

        end subroutine

        subroutine f_f_int_and_jac_structural_2d(B_L, S_v, C, dC, detj, dN_dX, f_int_K_dK_out, K_geo, no_of_nodes, &
        spatial_dimension, no_of_dofs)
            integer, intent(in) :: no_of_nodes
            integer, intent(in) :: no_of_dofs
            integer, intent(in) :: spatial_dimension
            real(dp), intent(in) :: B_L(3,no_of_dofs)
            real(dp), intent(in) :: dN_dX(no_of_nodes, spatial_dimension)
            real(dp), intent(in) :: S_v(3)
            real(dp), intent(in) :: detj
            real(dp), intent(in) :: C(3,3)
            real(dp), intent(in) :: dC(3,3)
            real(dp), intent(out) :: f_int_K_dK_out(no_of_dofs, 1+2*no_of_dofs)
            real(dp), intent(out) :: K_geo(no_of_dofs, no_of_dofs)

            real(dp) :: S_mat(2,2)
            real(dp) :: K_mat(no_of_dofs, no_of_dofs)
            real(dp) :: K_mat_velocity(no_of_dofs, no_of_dofs)
            real(dp) :: K_geo_small(no_of_nodes, no_of_nodes)

            f_int_K_dK_out(:, 1) = matmul(transpose(B_L),  S_v) * detj
            ! Compute Jacobian
            K_mat = matmul(matmul(transpose(B_L), C), B_L)
            K_mat_velocity = matmul(matmul(transpose(B_L), dC), B_L)
            S_mat(1, 1) = S_v(1)
            S_mat(1, 2) = S_v(3)
            S_mat(2, 1) = S_v(3)
            S_mat(2, 2) = S_v(2)
            K_geo(:, :) = 0.0
            K_geo_small = matmul(matmul(dN_dX, S_mat), transpose(dN_dX))

            call f_scatter_matrix(K_geo_small, K_geo, spatial_dimension, no_of_nodes, no_of_nodes)

            f_int_K_dK_out(:, 2:no_of_dofs+2) = (K_geo + K_mat) * detj
            f_int_K_dK_out(:, 2+no_of_dofs:2+2*no_of_dofs+1) = K_mat_velocity * detj
        end subroutine


        subroutine f_b_l_2d(F, dN_dX, out, no_of_nodes, spatial_dimension)
            integer, intent(in) :: no_of_nodes
            integer, intent(in) :: spatial_dimension
            real(dp), intent(in) :: F(2,2)
            real(dp), intent(in) :: dN_dX(no_of_nodes, spatial_dimension)
            real(dp), intent(inout) :: out(3, spatial_dimension*no_of_nodes)

            integer :: a

            do a=0,no_of_nodes-1
                out(1, 2 * a + 1) = F(1, 1) * dN_dX(a+1, 1)
                out(1, 2 * a + 2) = F(2, 1) * dN_dX(a+1, 1)
                out(2, 2 * a + 1) = F(1, 2) * dN_dX(a+1, 2)
                out(2, 2 * a + 2) = F(2, 2) * dN_dX(a+1, 2)
                out(3, 2 * a + 1) = F(1, 1) * dN_dX(a+1, 2) + F(1, 2) * dN_dX(a+1, 1)
                out(3, 2 * a + 2) = F(2, 1) * dN_dX(a+1, 2) + F(2, 2) * dN_dX(a+1, 1)
            end do
        end subroutine

        subroutine f_deformation_gradient_2d(du_dX, Fout)
            real(dp), intent(in) :: du_dX(2,2)
            real(dp), intent(inout) :: Fout(2,2)

            Fout = du_dX
            Fout(1, 1) = Fout(1,1) + 1.0_dp
            Fout(2, 2) = Fout(2,2) + 1.0_dp
        end subroutine

end module integrand_structural2d