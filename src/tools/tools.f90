!
! Copyright (c) 2020 TECHNICAL UNIVERSITY OF MUNICH, DEPARTMENT OF MECHANICAL ENGINEERING, CHAIR OF APPLIED MECHANICS,
! BOLTZMANNSTRASSE 15, 85748 GARCHING/MUNICH, GERMANY, RIXEN@TUM.DE.
!
! Distributed under 3-Clause BSD license. See LICENSE file for more information.
!


module tools
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    implicit none


    contains
        subroutine green_lagrange_strain_2d(F, Evout)
            real(dp), intent(in) :: F(2, 2)
            real(dp), intent(out) :: Evout(3)

            real(dp) :: eye(2, 2)
            real(dp) :: E(2, 2)

            eye = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp /), shape=(/2,2/))
            E = 0.5_dp*(matmul(transpose(F), F) - eye)
            Evout(1) = E(1,1)
            Evout(2) = E(2,2)
            Evout(3) = 2.0_dp * E(1,2)
        end subroutine green_lagrange_strain_2d

        subroutine green_lagrange_strain_3d(F, Evout)
            real(dp), intent(in) :: F(3,3)
            real(dp), intent(out) :: Evout(6)

            real(dp) :: eye(3,3)
            real(dp) :: E(3,3)

            eye = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp /), shape=(/3,3/))
            E = 0.5_dp*(matmul(transpose(F), F) - eye)

            Evout(1) = E(1, 1)
            Evout(2) = E(2, 2)
            Evout(3) = E(3, 3)
            Evout(4) = 2.0_dp * E(2, 3)
            Evout(5) = 2.0_dp * E(1, 3)
            Evout(6) = 2.0_dp * E(1, 2)
        end subroutine green_lagrange_strain_3d

        subroutine scatter_matrix(A, B, ndim, m, n)
            integer, intent(in) :: m
            integer, intent(in) :: n
            integer, intent(in) :: ndim
            real(dp), intent(in) :: A(m, n)
            real(dp), intent(inout) :: B(ndim*m, ndim*n)

            integer i, j, k

            !    Initialization to zero; otherwise random stuff is there...
            B(:,:) = 0.0_dp
            ! looping style like in Python:
            do j=0,n-1
                do i=0,m-1
                    do k=0,ndim-1
                        B(ndim*i+k+1,ndim*j+k+1) = A(i+1,j+1)
                    end do
                end do
            end do

        end subroutine scatter_matrix

        real(dp) function determinant_2d(mat)
            real(dp), intent(in) :: mat(2, 2)
            determinant_2d = mat(1, 1) * mat(2, 2) - mat(1, 2) * mat(2, 1)
        end function determinant_2d

        subroutine inverse_matrix_2d(mat, det, matout)
            real(dp), intent(in) :: mat(2, 2)
            real(dp), intent(in) :: det
            real(dp), intent(out) :: matout(2, 2)

            matout = reshape((/ mat(2, 2) / det, &
                                -mat(2, 1) / det, &
                                -mat(1, 2) / det, &
                                mat(1, 1) / det /), shape=(/2,2/))
        end subroutine inverse_matrix_2d


end module tools