!
! Copyright (c) 2020 TECHNICAL UNIVERSITY OF MUNICH, DEPARTMENT OF MECHANICAL ENGINEERING, CHAIR OF APPLIED MECHANICS,
! BOLTZMANNSTRASSE 15, 85748 GARCHING/MUNICH, GERMANY, RIXEN@TUM.DE.
!
! Distributed under 3-Clause BSD license. See LICENSE file for more information.
!

module integrand_structural3d_nonlinear_fun
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    use tools
    use fields

    implicit none
    private

    public igstrc3dnl_fields, igstrc3dnl_no_of_dofs, igstrc3dnl_spatial_dimension, igstrc3dnl_f_int_and_jac, &
            igstrc3dnl_m_int, igstrc3dnl_no_of_fields

    contains
        
        subroutine igstrc3dnl_fields(fields)
            integer, intent(out) :: fields(:)
            fields = (/ UX, UY, UZ /)
        end subroutine igstrc3dnl_fields

        integer function igstrc3dnl_no_of_fields() result(no_of_fields)
            no_of_fields = 3
        end function igstrc3dnl_no_of_fields

        integer function igstrc3dnl_no_of_dofs(no_of_nodes)
            integer, intent(in) :: no_of_nodes
            igstrc3dnl_no_of_dofs = 3*no_of_nodes
        end function igstrc3dnl_no_of_dofs

        integer function igstrc3dnl_spatial_dimension()
            igstrc3dnl_spatial_dimension = 3
        end function igstrc3dnl_spatial_dimension
        
        subroutine igstrc3dnl_f_int_and_jac(n, dn_dxi, xn, d, dd, t, matfunc_3d, s0, t0, s1, f_int_k_dk_out, &
        no_of_fields, spatial_dimension, no_of_nodes, no_of_states)
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
                subroutine sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out)
                    import :: dp
                    real(dp), intent(in) :: f(3,3)
                    real(dp), intent(in) :: df(3,3)
                    real(dp), intent(in) :: t
                    real(dp), intent(in) :: s0(:)
                    real(dp), intent(in) :: t0
                    real(dp), intent(out) :: s1(:)
                    real(dp), intent(out) :: sv_out(6)
                    real(dp), intent(out) :: c_out(6,6)
                    real(dp), intent(out) :: dc_out(6,6)
                end subroutine sv_c_dc_3d
            end interface
            procedure(sv_c_dc_3d) :: matfunc_3d
            real(dp), intent(in) :: s0(no_of_states)
            real(dp), intent(in) :: t0
            real(dp), intent(inout) :: s1(no_of_states)
            real(dp), intent(inout) :: f_int_k_dk_out(no_of_nodes*spatial_dimension,no_of_nodes*spatial_dimension*2+1)

            integer :: no_of_dofs
            real(dp) :: dn_dx(no_of_nodes,spatial_dimension)
            real(dp) :: du_dx(no_of_fields,spatial_dimension)
            real(dp) :: ddu_dx(no_of_fields,spatial_dimension)
            real(dp) :: detj(1)
            real(dp) :: fval(3,3)
            real(dp) :: dfval(3,3)
            real(dp) :: bval(6,spatial_dimension*no_of_nodes)
            real(dp) :: k_geo(no_of_nodes*spatial_dimension,no_of_nodes*spatial_dimension)
            real(dp) :: s_v(6)
            real(dp) :: c(6,6)
            real(dp) :: dc(6,6)

            no_of_dofs = no_of_nodes*spatial_dimension

            call f_material_gradients_structural_3d(n, dn_dxi, xn, d, dd, t, du_dx, ddu_dx, dn_dx, detj, &
            no_of_nodes, spatial_dimension, no_of_fields)

            call f_deformation_gradient_3d(du_dx, fval)
            call f_deformation_gradient_3d(ddu_dx, dfval)

            call matfunc_3d(fval,dfval,t,s0,t0,s1,s_v,c,dc)

            call f_b_l_3d(fval, dn_dx, bval, no_of_nodes, spatial_dimension)

            call f_f_int_and_jac_structural_3d(bval, s_v, c, dc, detj(1), dn_dx, f_int_k_dk_out, k_geo, no_of_nodes, &
                spatial_dimension, no_of_dofs)

        end subroutine igstrc3dnl_f_int_and_jac

        subroutine igstrc3dnl_m_int(n, dn_dxi, xn, d, dd, t, density, m_out, no_of_fields, &
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
            real(dp), intent(in) :: density
            real(dp), intent(inout) :: m_out(no_of_fields*no_of_nodes,no_of_fields*no_of_nodes)

            real(dp) :: x_mat(spatial_dimension, no_of_nodes)
            real(dp) :: j(spatial_dimension, spatial_dimension)
            real(dp) :: detj
            real(dp) :: m_small(no_of_nodes, no_of_nodes)
            integer :: k, m

            x_mat = reshape(xn, (/ spatial_dimension, no_of_nodes /))
            j = matmul(x_mat, dn_dxi)


            call f_determinant_3d(j, detj)

            do m=1,no_of_nodes
                do k=1,no_of_nodes
                    m_small(k,m) = n(k)*n(m)
                end do
            end do
            m_small = m_small * detJ * density

            call scatter_matrix(m_small, m_out, spatial_dimension, no_of_nodes, no_of_nodes)

        end subroutine igstrc3dnl_m_int

        subroutine f_material_gradients_structural_3d(n, dn_dxi, xn, d, dd, t, du_dx, ddu_dx, dn_dx, detj, &
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
            real(dp), intent(inout) :: dn_dx(no_of_nodes, spatial_dimension)
            real(dp), intent(inout) :: du_dx(no_of_fields, spatial_dimension)
            real(dp), intent(inout) :: ddu_dx(no_of_fields, spatial_dimension)
            real(dp), intent(inout) :: detj(1)

            real(dp) :: jinv(spatial_dimension, spatial_dimension)
            real(dp) :: d_mat(no_of_fields, no_of_nodes)
            real(dp) :: dd_mat(no_of_fields, no_of_nodes)
            real(dp) :: x_mat(spatial_dimension, no_of_nodes)
            real(dp) :: j(spatial_dimension, spatial_dimension)
            real(dp) :: detjval

            x_mat = reshape(xn, (/ spatial_dimension, no_of_nodes /))
            j = matmul(x_mat, dn_dxi)
            call f_determinant_3d(j, detjval)
            detj(1) = detjval
            call f_inverse_matrix_3d(j, detjval, jinv)

            dn_dX = matmul(dn_dxi, jinv)

            d_mat = reshape(d, (/ no_of_fields, no_of_nodes /))
            dd_mat = reshape(dd, (/ no_of_fields, no_of_nodes /))

            du_dx = matmul(d_mat, dn_dx)
            ddu_dx = matmul(dd_mat, dn_dx)

        end subroutine

        subroutine f_f_int_and_jac_structural_3d(B_L, S_v, C, dC, detJval, dN_dX, f_int_K_dK_out, K_geo, no_of_nodes, &
        spatial_dimension, no_of_dofs)
            integer, intent(in) :: no_of_nodes
            integer, intent(in) :: no_of_dofs
            integer, intent(in) :: spatial_dimension
            real(dp), intent(in) :: B_L(6,no_of_dofs)
            real(dp), intent(in) :: dN_dX(no_of_nodes, spatial_dimension)
            real(dp), intent(in) :: S_v(6)
            real(dp), intent(in) :: detJval
            real(dp), intent(in) :: C(6,6)
            real(dp), intent(in) :: dC(6,6)
            real(dp), intent(inout) :: f_int_K_dK_out(no_of_dofs, 1+2*no_of_dofs)
            real(dp), intent(inout) :: K_geo(no_of_dofs, no_of_dofs)

            real(dp) :: S_mat(3,3)
            real(dp) :: K_mat(no_of_dofs, no_of_dofs)
            real(dp) :: K_mat_velocity(no_of_dofs, no_of_dofs)
            real(dp) :: K_geo_small(no_of_nodes, no_of_nodes)

            f_int_K_dK_out(:, 1) = matmul(transpose(B_L),  S_v) * detJval
            ! Compute Jacobian
            K_mat = matmul(matmul(transpose(B_L), C), B_L)
            K_mat_velocity = matmul(matmul(transpose(B_L), dC), B_L)
            S_mat(1, 1) = S_v(1)
            S_mat(2, 2) = S_v(2)
            S_mat(3, 3) = S_v(3)
            S_mat(2, 3) = S_v(4)
            S_mat(3, 2) = S_v(4)
            S_mat(1, 3) = S_v(5)
            S_mat(3, 1) = S_v(5)
            S_mat(1, 2) = S_v(6)
            S_mat(2, 1) = S_v(6)

            K_geo(:, :) = 0.0
            K_geo_small = matmul(matmul(dN_dX, S_mat), transpose(dN_dX))

            call scatter_matrix(K_geo_small, K_geo, spatial_dimension, no_of_nodes, no_of_nodes)

            f_int_K_dK_out(:, 2:no_of_dofs+2) = (K_geo + K_mat) * detJval
            f_int_K_dK_out(:, 2+no_of_dofs:2+2*no_of_dofs+1) = K_mat_velocity * detJval
        end subroutine


        subroutine f_b_l_3d(F, dN_dX, out, no_of_nodes, spatial_dimension)
            integer, intent(in) :: no_of_nodes
            integer, intent(in) :: spatial_dimension
            real(dp), intent(in) :: F(3,3)
            real(dp), intent(in) :: dN_dX(no_of_nodes, spatial_dimension)
            real(dp), intent(inout) :: out(6, spatial_dimension*no_of_nodes)

            integer :: a

            do a=0,no_of_nodes-1
                out(1, 3 * a + 1) = F(1, 1) * dN_dX(a+1, 1)
                out(1, 3 * a + 2) = F(2, 1) * dN_dX(a+1, 1)
                out(1, 3 * a + 3) = F(3, 1) * dN_dX(a+1, 1)
                out(2, 3 * a + 1) = F(1, 2) * dN_dX(a+1, 2)
                out(2, 3 * a + 2) = F(2, 2) * dN_dX(a+1, 2)
                out(2, 3 * a + 3) = F(3, 2) * dN_dX(a+1, 2)
                out(3, 3 * a + 1) = F(1, 3) * dN_dX(a+1, 3)
                out(3, 3 * a + 2) = F(2, 3) * dN_dX(a+1, 3)
                out(3, 3 * a + 3) = F(3, 3) * dN_dX(a+1, 3)
                out(4, 3 * a + 1) = F(1, 2) * dN_dX(a+1, 3) + F(1, 3) * dN_dX(a+1, 2)
                out(4, 3 * a + 2) = F(2, 2) * dN_dX(a+1, 3) + F(2, 3) * dN_dX(a+1, 2)
                out(4, 3 * a + 3) = F(3, 2) * dN_dX(a+1, 3) + F(3, 3) * dN_dX(a+1, 2)
                out(5, 3 * a + 1) = F(1, 3) * dN_dX(a+1, 1) + F(1, 1) * dN_dX(a+1, 3)
                out(5, 3 * a + 2) = F(2, 3) * dN_dX(a+1, 1) + F(2, 1) * dN_dX(a+1, 3)
                out(5, 3 * a + 3) = F(3, 3) * dN_dX(a+1, 1) + F(3, 1) * dN_dX(a+1, 3)
                out(6, 3 * a + 1) = F(1, 1) * dN_dX(a+1, 2) + F(1, 2) * dN_dX(a+1, 1)
                out(6, 3 * a + 2) = F(2, 1) * dN_dX(a+1, 2) + F(2, 2) * dN_dX(a+1, 1)
                out(6, 3 * a + 3) = F(3, 1) * dN_dX(a+1, 2) + F(3, 2) * dN_dX(a+1, 1)
            end do

        end subroutine

        subroutine f_determinant_3d(mat, det)
            real(dp), intent(in) :: mat(3, 3)
            real(dp), intent(out) :: det

            real(dp) :: a11, a12, a13, a21, a22, a23, a31, a32, a33

            a11 = mat(1,1)
            a12 = mat(1,2)
            a13 = mat(1,3)
            a21 = mat(2,1)
            a22 = mat(2,2)
            a23 = mat(2,3)
            a31 = mat(3,1)
            a32 = mat(3,2)
            a33 = mat(3,3)
            det = a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31
        end subroutine

        subroutine f_inverse_matrix_3d(mat, det, matout)
            real(dp), intent(in) :: mat(3, 3)
            real(dp), intent(in) :: det
            real(dp), intent(inout) :: matout(3, 3)

            real(dp) :: a11, a12, a13, a21, a22, a23, a31, a32, a33
            a11 = mat(1,1)
            a12 = mat(1,2)
            a13 = mat(1,3)
            a21 = mat(2,1)
            a22 = mat(2,2)
            a23 = mat(2,3)
            a31 = mat(3,1)
            a32 = mat(3,2)
            a33 = mat(3,3)
            matout = reshape((/ a22*a33 - a23*a32, -a21*a33 + a23*a31, &
                               a21*a32 - a22*a31, -a12*a33 + a13*a32, &
                               a11*a33 - a13*a31, -a11*a32 + a12*a31, &
                               a12*a23 - a13*a22, -a11*a23 + a13*a21, &
                               a11*a22 - a12*a21 /), (/3,3/)) / det
        end subroutine

        subroutine f_deformation_gradient_3d(du_dX, Fout)
            real(dp), intent(in) :: du_dX(3,3)
            real(dp), intent(inout) :: Fout(3,3)

            Fout = du_dX
            Fout(1, 1) = Fout(1,1) + 1.0_dp
            Fout(2, 2) = Fout(2,2) + 1.0_dp
            Fout(3, 3) = Fout(3,3) + 1.0_dp

        end subroutine

end module integrand_structural3d_nonlinear_fun