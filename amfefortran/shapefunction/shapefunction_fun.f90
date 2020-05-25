!
! Copyright (c) 2020 TECHNICAL UNIVERSITY OF MUNICH, DEPARTMENT OF MECHANICAL ENGINEERING, CHAIR OF APPLIED MECHANICS,
! BOLTZMANNSTRASSE 15, 85748 GARCHING/MUNICH, GERMANY, RIXEN@TUM.DE.
!
! Distributed under 3-Clause BSD license. See LICENSE file for more information.
!


module shapefunc
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    implicit none

    contains

        ! -------------------------------------------------- Line 2 ---------------------------------------------------
        integer function line2_no_of_nodes()
            line2_no_of_nodes = 2
        end function line2_no_of_nodes

        integer function line2_no_of_local_coordinates()
            line2_no_of_local_coordinates = 1
        end function line2_no_of_local_coordinates

        real(dp) function line2_xi_lower()
            line2_xi_lower = -1.0_dp
        end function line2_xi_lower

        real(dp) function line2_xi_upper()
            line2_xi_upper = 1.0_dp
        end function line2_xi_upper

        subroutine line2_evaluate(xi, out)
            real(dp), intent(in) :: xi
            real(dp), intent(out) :: out(2)
            !real(dp), pointer, intent(out) :: out(:)

            out = (/ 0.5_dp * (1.0_dp - xi), &
                     0.5_dp * (1.0_dp + xi) /)
        end subroutine line2_evaluate

        subroutine line2_jacobian(xi, out)
            real(dp), intent(in) :: xi
            real(dp), intent(out) :: out(2,1)
            out = reshape((/-0.5_dp, 0.5_dp/), shape=(/2,1/))
        end subroutine line2_jacobian

        !------------------------------------------------ Line 3 ------------------------------------------------------

        integer function line3_no_of_nodes()
            line3_no_of_nodes = 3
        end function line3_no_of_nodes

        integer function line3_no_of_local_coordinates()
            line3_no_of_local_coordinates = 1
        end function line3_no_of_local_coordinates

        real(dp) function line3_xi_lower()
            line3_xi_lower = -1.0_dp
        end function line3_xi_lower

        real(dp) function line3_xi_upper()
            line3_xi_upper = 1.0_dp
        end function line3_xi_upper

        subroutine line3_evaluate(xi, out)
            real(dp), intent(in) :: xi
            real(dp), intent(out) :: out(3)

            out = (/ -0.5_dp * xi * (1.0_dp - xi), &
                             0.5_dp * xi * (1.0_dp + xi), &
                             (1.0_dp - xi**2) /)
        end subroutine line3_evaluate

        subroutine line3_jacobian(xi, out)
            real(dp), intent(in) :: xi
            real(dp), intent(out) :: out(3,1)

            out = reshape((/1.0_dp*xi - 0.5_dp, 1.0_dp*xi + 0.5_dp, -2.0_dp*xi/), shape=(/3,1/))
        end subroutine line3_jacobian

        !------------------------------------------ Tri3 --------------------------------------------------------------
        integer function tri3_no_of_nodes()
            tri3_no_of_nodes = 3
        end function tri3_no_of_nodes

        integer function tri3_no_of_local_coordinates()
            tri3_no_of_local_coordinates = 2
        end function tri3_no_of_local_coordinates

        real(dp) function tri3_xi_lower()
            tri3_xi_lower = 0.0_dp
        end function tri3_xi_lower

        real(dp) function tri3_xi_upper()
            tri3_xi_upper = 1.0_dp
        end function tri3_xi_upper

        real(dp) function tri3_eta_lower(xi)
            real(dp), intent(in) :: xi
            tri3_eta_lower = 0.0_dp
        end function tri3_eta_lower

        real(dp) function tri3_eta_upper(xi)
            real(dp), intent(in) :: xi
            tri3_eta_upper = 1.0_dp - xi
        end function tri3_eta_upper

        subroutine tri3_evaluate(xi, eta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: out(3)

            out = (/ 1.0_dp - xi - eta, &
                              xi, &
                              eta /)
        end subroutine tri3_evaluate

        subroutine tri3_jacobian(xi, eta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: out(3,2)

            out = reshape([ -1.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp ], [3, 2])
        end subroutine tri3_jacobian

        !------------------------------------------------ Quad4 -------------------------------------------------------
        integer function quad4_no_of_nodes()
            quad4_no_of_nodes = 4
        end function quad4_no_of_nodes

        integer function quad4_no_of_local_coordinates()
            quad4_no_of_local_coordinates = 2
        end function quad4_no_of_local_coordinates

        real(dp) function quad4_xi_lower()
            quad4_xi_lower = -1.0_dp
        end function quad4_xi_lower

        real(dp) function quad4_xi_upper()
            quad4_xi_upper = 1.0_dp
        end function quad4_xi_upper

        real(dp) function quad4_eta_lower(xi)
            real(dp), intent(in) :: xi
            quad4_eta_lower = -1.0_dp
        end function quad4_eta_lower

        real(dp) function quad4_eta_upper(xi)
            real(dp), intent(in) :: xi
            quad4_eta_upper = 1.0_dp
        end function quad4_eta_upper

        subroutine quad4_evaluate(xi, eta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: out(4)

            out = (/ 0.25_dp * (1.0_dp - xi) * (1.0_dp - eta), &
                            0.25_dp * (1.0_dp + xi) * (1.0_dp - eta), &
                            0.25_dp * (1.0_dp + xi) * (1.0_dp + eta), &
                            0.25_dp * (1.0_dp - xi) * (1.0_dp + eta) /)
        end subroutine quad4_evaluate

        subroutine quad4_jacobian(xi, eta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: out(4,2)

            out = reshape((/ 0.25_dp * (eta - 1.0_dp), &
                            0.25_dp * (1.0_dp - eta), &
                            0.25_dp * (1.0_dp + eta), &
                            0.25_dp * (-1.0_dp - eta), &
                            0.25_dp * (xi - 1.0_dp), &
                            0.25_dp * (-xi - 1.0_dp), &
                            0.25_dp * (1.0_dp + xi), &
                            0.25_dp * (1.0_dp - xi)/), &
                            shape=(/4,2/))
        end subroutine quad4_jacobian

        !---------------------------------------------------------- Quad8 ---------------------------------------------
        integer function quad8_no_of_nodes()
            quad8_no_of_nodes = 8
        end function quad8_no_of_nodes

        integer function quad8_no_of_local_coordinates()
            quad8_no_of_local_coordinates = 2
        end function quad8_no_of_local_coordinates

        real(dp) function quad8_xi_lower()
            quad8_xi_lower = -1.0_dp
        end function quad8_xi_lower

        real(dp) function quad8_xi_upper()
            quad8_xi_upper = 1.0_dp
        end function quad8_xi_upper

        real(dp) function quad8_eta_lower(xi)
            real(dp), intent(in) :: xi
            quad8_eta_lower = -1.0_dp
        end function quad8_eta_lower

        real(dp) function quad8_eta_upper(xi)
            real(dp), intent(in) :: xi
            quad8_eta_upper = 1.0_dp
        end function quad8_eta_upper

        subroutine quad8_evaluate(xi, eta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: out(8)

            out = (/0.25_dp * ((1.0_dp + xi) * (1.0_dp + eta) - (1.0_dp - xi**2) * (1.0_dp + eta) - &
            (1.0_dp - eta**2) * (1.0_dp + xi)), &
                          0.25_dp *((1.0_dp - xi) * (1.0_dp + eta) - (1.0_dp - xi ** 2) * (1.0_dp + eta) &
                                  - (1.0_dp - eta ** 2) * (1.0_dp - xi)),&
                          0.25_dp *((1.0_dp - xi) * (1.0_dp - eta) - (1.0_dp - eta ** 2) * (1.0_dp - xi) &
                                  - (1.0_dp - xi ** 2) * (1.0_dp - eta)),&
                          0.25_dp *((1.0_dp + xi) * (1.0_dp - eta) - (1.0_dp - xi ** 2) * (1.0_dp - eta) &
                                  - (1.0_dp - eta ** 2) * (1.0_dp + xi)),&
                          0.5_dp * (1.0_dp - xi ** 2) * (1.0_dp + eta),&
                          0.5_dp * (1.0_dp - eta ** 2) * (1.0_dp - xi),&
                          0.5_dp * (1.0_dp - xi ** 2) * (1.0_dp - eta),&
                          0.5_dp * (1.0_dp - eta ** 2) * (1.0_dp + xi) /)
        end subroutine quad8_evaluate

        subroutine quad8_jacobian(xi, eta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: out(8,2)

            out = reshape((/ 0.25_dp * eta ** 2 + 0.25_dp * eta + 0.5_dp * xi * (eta + 1.0_dp), &
            -0.25_dp * eta ** 2 - 0.25_dp * eta + 0.5_dp * xi * (eta + 1.0_dp), &
            -0.25_dp * eta ** 2 + 0.25_dp * eta + 0.5_dp * xi * (1.0_dp - eta), &
            0.25_dp * eta ** 2 - 0.25_dp * eta + 0.5_dp * xi * (1.0_dp - eta), &
            -1.0_dp * xi * (eta + 1.0_dp), &
            0.5_dp * eta ** 2 - 0.5_dp, &
            -1.0_dp * xi * (1.0_dp - eta), &
            0.5_dp - 0.5_dp * eta ** 2, &
            0.5_dp * eta * (xi + 1.0_dp) + 0.25_dp * xi ** 2 + 0.25_dp * xi, &
            0.5_dp * eta * (1.0_dp - xi) + 0.25_dp * xi ** 2 - 0.25_dp * xi, &
            0.5_dp * eta * (1.0_dp - xi) - 0.25_dp * xi ** 2 + 0.25_dp * xi, &
            0.5_dp * eta * (xi + 1.0_dp) - 0.25_dp * xi ** 2 - 0.25_dp * xi, &
            0.5_dp - 0.5_dp * xi ** 2, &
            -1.0_dp * eta * (1.0_dp - xi), &
            0.5_dp * xi ** 2 - 0.5_dp, &
            -1.0_dp * eta * (xi + 1.0_dp) /), &
            shape=(/8, 2/))
        end subroutine quad8_jacobian

        !---------------------------------------------------------- Tri6 ---------------------------------------------
        integer function tri6_no_of_nodes()
            tri6_no_of_nodes = 6
        end function tri6_no_of_nodes

        integer function tri6_no_of_local_coordinates()
            tri6_no_of_local_coordinates = 2
        end function tri6_no_of_local_coordinates

        real(dp) function tri6_xi_lower()
            tri6_xi_lower = 0.0_dp
        end function tri6_xi_lower

        real(dp) function tri6_xi_upper()
            tri6_xi_upper = 1.0_dp
        end function tri6_xi_upper

        real(dp) function tri6_eta_lower(xi)
            real(dp), intent(in) :: xi
            tri6_eta_lower = 0.0_dp
        end function tri6_eta_lower

        real(dp) function tri6_eta_upper(xi)
            real(dp), intent(in) :: xi
            tri6_eta_upper = 1.0_dp - xi
        end function tri6_eta_upper

        subroutine tri6_evaluate(xi, eta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: out(6)

            out = (/ (1.0_dp - xi - eta) * (1.0_dp - 2.0_dp * (xi + eta)), &
                           - 2.0_dp * xi * (0.5_dp - xi), &
                           - 2.0_dp * eta * (0.5_dp - eta), &
                             4.0_dp * xi * (1.0_dp - xi - eta), &
                             4.0_dp * xi * eta, &
                             4.0_dp * eta * (1.0_dp - xi - eta) /)
        end subroutine tri6_evaluate

        subroutine tri6_jacobian(xi, eta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: out(6,2)

            out = reshape((/ 4.0_dp*eta + 4.0_dp*xi - 3.0_dp, &
                            4.0_dp*xi - 1.0_dp, &
                            0.0_dp, &
                            -4.0_dp*eta - 8.0_dp*xi + 4.0_dp, &
                            4.0_dp * eta, &
                            -4.0_dp*eta, &
                            4.0_dp*eta + 4.0_dp*xi - 3.0_dp, &
                            0.0_dp, &
                            4.0_dp*eta - 1.0_dp, &
                            -4.0_dp*xi, &
                            4.0_dp * xi, &
                            -8.0_dp*eta - 4.0_dp*xi + 4.0_dp /), &
                            shape = (/6, 2/))
        end subroutine tri6_jacobian

        !---------------------------------------------------------- Hexa8 ---------------------------------------------
        integer function hexa8_no_of_nodes()
            hexa8_no_of_nodes = 8
        end function hexa8_no_of_nodes

        integer function hexa8_no_of_local_coordinates()
            hexa8_no_of_local_coordinates = 3
        end function hexa8_no_of_local_coordinates

        real(dp) function hexa8_xi_lower()
            hexa8_xi_lower = -1.0_dp
        end function hexa8_xi_lower

        real(dp) function hexa8_xi_upper()
            hexa8_xi_upper = 1.0_dp
        end function hexa8_xi_upper

        real(dp) function hexa8_eta_lower(xi)
            real(dp), intent(in) :: xi
            hexa8_eta_lower = -1.0_dp
        end function hexa8_eta_lower

        real(dp) function hexa8_eta_upper(xi)
            real(dp), intent(in) :: xi
            hexa8_eta_upper = 1.0_dp
        end function hexa8_eta_upper

        real(dp) function hexa8_zeta_lower(xi, eta)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            hexa8_zeta_lower = -1.0_dp
        end function hexa8_zeta_lower

        real(dp) function hexa8_zeta_upper(xi, eta)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            hexa8_zeta_upper = 1.0_dp
        end function hexa8_zeta_upper

        subroutine hexa8_evaluate(xi, eta, zeta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(in) :: zeta
            real(dp), intent(out) :: out(8)

            out = (/ (-eta + 1.0_dp) * (-xi + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
                            (-eta + 1.0_dp) * (xi + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
                            (eta + 1.0_dp) * (xi + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
                            (eta + 1.0_dp) * (-xi + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
                            (-eta + 1.0_dp) * (-xi + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
                            (-eta + 1.0_dp) * (xi + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
                            (eta + 1.0_dp) * (xi + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
                            (eta + 1.0_dp) * (-xi + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp /)
        end subroutine hexa8_evaluate

        subroutine hexa8_jacobian(xi, eta, zeta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(in) :: zeta
            real(dp), intent(out) :: out(8,3)

            out = reshape((/ -(-eta + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
            (-eta + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
            (eta + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
            -(eta + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
            -(-eta + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
            (-eta + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
            (eta + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
            -(eta + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
            -(-xi + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
            -(xi + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
            (xi + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
            (-xi + 1.0_dp) * (-zeta + 1.0_dp) / 8.0_dp, &
            -(-xi + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
            -(xi + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
            (xi + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
            (-xi + 1.0_dp) * (zeta + 1.0_dp) / 8.0_dp, &
            -(-eta + 1.0_dp) * (-xi + 1.0_dp) / 8.0_dp, &
            -(-eta + 1.0_dp) * (xi + 1.0_dp) / 8.0_dp, &
            -(eta + 1.0_dp) * (xi + 1.0_dp) / 8.0_dp, &
            -(eta + 1.0_dp) * (-xi + 1.0_dp) / 8.0_dp, &
            (-eta + 1.0_dp) * (-xi + 1.0_dp) / 8.0_dp, &
            (-eta + 1.0_dp) * (xi + 1.0_dp) / 8.0_dp, &
            (eta + 1.0_dp) * (xi + 1.0_dp) / 8.0_dp, &
            (eta + 1.0_dp) * (-xi + 1.0_dp) / 8.0_dp /), &
            shape=(/8,3/))
        end subroutine hexa8_jacobian

        !---------------------------------------------------------- Tet4 ---------------------------------------------
        integer function tet4_no_of_nodes()
            tet4_no_of_nodes = 4
        end function tet4_no_of_nodes

        integer function tet4_no_of_local_coordinates()
            tet4_no_of_local_coordinates = 3
        end function tet4_no_of_local_coordinates

        real(dp) function tet4_xi_lower()
            tet4_xi_lower = 0.0_dp
        end function tet4_xi_lower

        real(dp) function tet4_xi_upper()
            tet4_xi_upper = 1.0_dp
        end function tet4_xi_upper

        real(dp) function tet4_eta_lower(xi)
            real(dp), intent(in) :: xi
            tet4_eta_lower = 0.0_dp
        end function tet4_eta_lower

        real(dp) function tet4_eta_upper(xi)
            real(dp), intent(in) :: xi
            tet4_eta_upper = 1.0_dp - xi
        end function tet4_eta_upper

        real(dp) function tet4_zeta_lower(xi, eta)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            tet4_zeta_lower = 0.0_dp
        end function tet4_zeta_lower

        real(dp) function tet4_zeta_upper(xi, eta)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            tet4_zeta_upper = 1.0_dp - xi - eta
        end function tet4_zeta_upper

        subroutine tet4_evaluate(xi, eta, zeta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(in) :: zeta
            real(dp), intent(out) :: out(4)

            out = (/ 1.0_dp - xi - eta - zeta, &
                           xi, &
                           eta, &
                           zeta /)
        end subroutine tet4_evaluate

        subroutine tet4_jacobian(xi, eta, zeta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(in) :: zeta
            real(dp), intent(out) :: out(4,3)

            out = reshape((/-1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp,&
                    1.0_dp/),(/4,3/))
        end subroutine tet4_jacobian


        !---------------------------------------------------------- Hexa20 --------------------------------------------
        integer function hexa20_no_of_nodes()
            hexa20_no_of_nodes = 20
        end function hexa20_no_of_nodes

        integer function hexa20_no_of_local_coordinates()
            hexa20_no_of_local_coordinates = 3
        end function hexa20_no_of_local_coordinates

        real(dp) function hexa20_xi_lower()
            hexa20_xi_lower = -1.0_dp
        end function hexa20_xi_lower

        real(dp) function hexa20_xi_upper()
            hexa20_xi_upper = 1.0_dp
        end function hexa20_xi_upper

        real(dp) function hexa20_eta_lower(xi)
            real(dp), intent(in) :: xi
            hexa20_eta_lower = -1.0_dp
        end function hexa20_eta_lower

        real(dp) function hexa20_eta_upper(xi)
            real(dp), intent(in) :: xi
            hexa20_eta_upper = 1.0_dp
        end function hexa20_eta_upper

        real(dp) function hexa20_zeta_lower(xi, eta)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            hexa20_zeta_lower = -1.0_dp
        end function hexa20_zeta_lower

        real(dp) function hexa20_zeta_upper(xi, eta)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            hexa20_zeta_upper = 1.0_dp
        end function hexa20_zeta_upper

        subroutine hexa20_evaluate(xi, eta, zeta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(in) :: zeta
            real(dp), intent(out) :: out(20)

            out = (/ (1.0_dp - eta) * (1.0_dp - xi) * (1.0_dp - zeta) * (-eta - xi - zeta - 2.0_dp) / 8.0_dp, &
            (1.0_dp - eta) * (1.0_dp - zeta) * (xi + 1.0_dp) * (-eta + xi - zeta - 2.0_dp) / 8.0_dp, &
            (1.0_dp - zeta) * (eta + 1.0_dp) * (xi + 1.0_dp) * (eta + xi - zeta - 2.0_dp) / 8.0_dp, &
            (1.0_dp - xi) * (1.0_dp - zeta) * (eta + 1.0_dp) * (eta - xi - zeta - 2.0_dp) / 8.0_dp, &
            (1.0_dp - eta) * (1.0_dp - xi) * (zeta + 1.0_dp) * (-eta - xi + zeta - 2.0_dp) / 8.0_dp, &
            (1.0_dp - eta) * (xi + 1.0_dp) * (zeta + 1.0_dp) * (-eta + xi + zeta - 2.0_dp) / 8.0_dp, &
            (eta + 1.0_dp) * (xi + 1.0_dp) * (zeta + 1.0_dp) * (eta + xi + zeta - 2.0_dp) / 8.0_dp, &
            (1.0_dp - xi) * (eta + 1.0_dp) * (zeta + 1.0_dp) * (eta - xi + zeta - 2.0_dp) / 8.0_dp, &
            (0.25_dp - 0.25_dp * xi ** 2) * (1.0_dp - eta) * (1.0_dp - zeta), &
            (1.0_dp - eta ** 2) * (1.0_dp - zeta) * (0.25_dp * xi + 0.25), &
            (0.25_dp - 0.25_dp * xi ** 2) * (1.0_dp - zeta) * (eta + 1.0_dp), &
            (0.25_dp - 0.25_dp * xi) * (1.0_dp - eta ** 2) * (1.0_dp - zeta), &
            (0.25_dp - 0.25_dp * xi ** 2) * (1.0_dp - eta) * (zeta + 1.0_dp), &
            (1.0_dp - eta ** 2) * (0.25_dp * xi + 0.25) * (zeta + 1.0_dp), &
            (0.25_dp - 0.25_dp * xi ** 2) * (eta + 1.0_dp) * (zeta + 1.0_dp), &
            (0.25_dp - 0.25_dp * xi) * (1.0_dp - eta ** 2) * (zeta + 1.0_dp), &
            (0.25_dp - 0.25_dp * xi) * (1.0_dp - eta) * (1.0_dp - zeta ** 2), &
            (1.0_dp - eta) * (1.0_dp - zeta ** 2) * (0.25_dp * xi + 0.25), &
            (1.0_dp - zeta ** 2) * (eta + 1.0_dp) * (0.25_dp * xi + 0.25), &
            (0.25_dp - 0.25_dp * xi) * (1.0_dp - zeta ** 2) * (eta + 1.0_dp) /)
        end subroutine hexa20_evaluate

        subroutine hexa20_jacobian(xi, eta, zeta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(in) :: zeta
            real(dp), intent(out) :: out(20,3)

            out = reshape((/ -(1.0_dp / 8.0_dp - eta / 8.0_dp) * (1.0_dp - xi) * (1.0_dp - zeta) + (1.0_dp - zeta) * ( &
            eta / 8.0_dp - 1.0_dp / 8.0_dp) * (-eta - xi - zeta - 2.0_dp), &
            (1.0_dp / 8.0_dp - eta / 8.0_dp) * (1.0_dp - zeta) * (-eta + xi - zeta - 2.0_dp) + (1.0_dp - eta) * ( &
            1.0_dp - zeta) * (xi / 8.0_dp + 1.0_dp / 8.0_dp), &
            (1.0_dp - zeta) * (eta / 8.0_dp + 1.0_dp / 8.0_dp) * (xi + 1.0_dp) + (1.0_dp - zeta) * &
            (eta / 8.0_dp + 1.0_dp / 8.0_dp) * ( eta + xi - zeta - 2.0_dp), &
            -(1.0_dp - xi) * (1.0_dp - zeta) * (eta / 8.0_dp + 1.0_dp / 8.0_dp) + (1.0_dp - zeta) * ( &
            -eta / 8.0_dp - 1.0_dp / 8.0_dp) * (eta - xi - zeta - 2.0_dp), &
            -(1.0_dp - eta) * (1.0_dp - xi) * (zeta / 8.0_dp + 1.0_dp / 8.0_dp) - (1.0_dp - eta) * ( &
            zeta / 8.0_dp + 1.0_dp / 8.0_dp) * (-eta - xi + zeta - 2.0_dp), &
            (1.0_dp - eta) * (xi / 8.0_dp + 1.0_dp / 8.0_dp) * (zeta + 1.0_dp) + (1.0_dp - eta) * &
            (zeta / 8.0_dp + 1.0_dp / 8.0_dp) * (-eta + xi + zeta - 2.0_dp), &
            (eta / 8.0_dp + 1.0_dp / 8.0_dp) * (xi + 1) * (zeta + 1.0_dp) + (eta / 8.0_dp + 1.0_dp / 8.0_dp) * &
            (zeta + 1.0_dp) * ( eta + xi + zeta - 2.0_dp), &
            -(1.0_dp - xi) * (eta / 8.0_dp + 1.0_dp / 8.0_dp) * (zeta + 1.0_dp) - (eta / 8.0_dp + 1.0_dp / 8.0_dp) * &
            ( zeta + 1.0_dp) * (eta - xi + zeta - 2.0_dp), &
            -xi * (1.0_dp - eta) * (1.0_dp - zeta) / 2.0_dp, &
            (1.0_dp - eta ** 2) * (1.0_dp - zeta) / 4.0_dp, &
            -xi * (1.0_dp - zeta) * (eta + 1) / 2.0_dp, &
            -(1.0_dp - eta ** 2) * (1.0_dp - zeta) / 4.0_dp, &
            -xi * (1.0_dp - eta) * (zeta + 1.0_dp) / 2.0_dp, &
            (1.0_dp - eta ** 2) * (zeta + 1.0_dp) / 4.0_dp, &
            -xi * (eta + 1) * (zeta + 1) / 2.0_dp, &
            -(1.0_dp - eta ** 2) * (zeta + 1.0_dp) / 4.0_dp, &
            -(1.0_dp - eta) * (1.0_dp - zeta ** 2) / 4.0_dp, &
            (1.0_dp - eta) * (1.0_dp - zeta ** 2) / 4.0_dp, &
            (1.0_dp - zeta ** 2) * (eta + 1.0_dp) / 4.0_dp, &
            -(1.0_dp - zeta ** 2) * (eta + 1.0_dp) / 4.0_dp, &
            -(1.0_dp / 8.0_dp - eta / 8.0_dp) * (1.0_dp - xi) * (1.0_dp - zeta) + (1.0_dp - zeta) * &
            (xi / 8.0_dp - 1.0_dp / 8.0_dp) * ( -eta - xi - zeta - 2.0_dp), &
            -(1.0_dp - eta) * (1.0_dp - zeta) * (xi / 8.0_dp + 1.0_dp / 8.0_dp) + (1.0_dp - zeta) * ( &
            -xi / 8.0_dp - 1.0_dp / 8.0_dp) * (-eta + xi - zeta - 2.0_dp), &
            (1.0_dp - zeta) * (eta / 8.0_dp + 1.0_dp / 8.0_dp) * (xi + 1.0_dp) + (1.0_dp - zeta) * &
            (xi / 8.0_dp + 1.0_dp / 8.0_dp) * ( eta + xi - zeta - 2.0_dp), &
            (1.0_dp / 8.0_dp - xi / 8) * (1.0_dp - zeta) * (eta - xi - zeta - 2.0_dp) + (1.0_dp - xi) * &
            (1.0_dp - zeta) * ( eta / 8.0_dp + 1.0_dp / 8.0_dp), &
            -(1.0_dp - eta) * (1.0_dp - xi) * (zeta / 8.0_dp + 1.0_dp / 8.0_dp) + (1.0_dp - xi) * ( &
            -zeta / 8.0_dp - 1.0_dp / 8.0_dp) * (-eta - xi + zeta - 2.0_dp), &
            -(1.0_dp - eta) * (xi / 8.0_dp + 1.0_dp / 8.0_dp) * (zeta + 1.0_dp) - (xi / 8.0_dp + 1.0_dp / 8.0_dp) * &
            (zeta + 1.0_dp) * ( -eta + xi + zeta - 2.0_dp), &
            (eta / 8.0_dp + 1.0_dp / 8.0_dp) * (xi + 1) * (zeta + 1.0_dp) + (xi / 8.0_dp + 1.0_dp / 8.0_dp) * &
            (zeta + 1.0_dp) * ( eta + xi + zeta - 2.0_dp), &
            (1.0_dp - xi) * (eta / 8.0_dp + 1.0_dp / 8.0_dp) * (zeta + 1.0_dp) + (1.0_dp - xi) * &
            (zeta / 8.0_dp + 1.0_dp / 8.0_dp) * ( eta - xi + zeta - 2.0_dp), &
            (0.25_dp - xi ** 2 / 4.0_dp) * (zeta - 1.0_dp), &
            -2.0_dp * eta * (1.0_dp - zeta) * (xi / 4.0_dp + 0.25), &
            (0.25_dp - xi ** 2 / 4.0_dp) * (1.0_dp - zeta), &
            -2.0_dp * eta * (0.25_dp - xi / 4) * (1.0_dp - zeta), &
            (0.25_dp - xi ** 2 / 4.0_dp) * (-zeta - 1.0_dp), &
            -2.0_dp * eta * (xi / 4 + 0.25) * (zeta + 1.0_dp), &
            (0.25_dp - xi ** 2 / 4.0_dp) * (zeta + 1.0_dp), &
            -2.0_dp * eta * (0.25_dp - xi / 4.0_dp) * (zeta + 1.0_dp), &
            (0.25_dp - xi / 4.0_dp) * (zeta ** 2 - 1.0_dp), &
            (xi / 4.0_dp + 0.25) * (zeta ** 2 - 1.0_dp), &
            (1.0_dp - zeta ** 2) * (xi / 4.0_dp + 0.25), &
            (0.25_dp - xi / 4) * (1.0_dp - zeta ** 2), &
            -(1.0_dp / 8.0_dp - eta / 8.0_dp) * (1.0_dp - xi) * (1.0_dp - zeta) - (1.0_dp / 8.0_dp - eta / 8.0_dp) * &
                    (1.0_dp - xi) * ( -eta - xi - zeta - 2.0_dp), &
            -(1.0_dp - eta) * (1.0_dp - zeta) * (xi / 8.0_dp + 1.0_dp / 8.0_dp) - (1.0_dp - eta) * &
                    (xi / 8.0_dp + 1.0_dp / 8.0_dp) * ( -eta + xi - zeta - 2.0_dp), &
            -(1.0_dp - zeta) * (eta / 8.0_dp + 1.0_dp / 8.0_dp) * (xi + 1.0_dp) - (eta / 8.0_dp + 1.0_dp / 8.0_dp) * &
                    (xi + 1) * ( eta + xi - zeta - 2.0_dp), &
            -(1.0_dp - xi) * (1.0_dp - zeta) * (eta / 8.0_dp + 1.0_dp / 8.0_dp) - (1.0_dp - xi) * &
                    (eta / 8.0_dp + 1.0_dp / 8.0_dp) * ( eta - xi - zeta - 2.0_dp), &
            (1.0_dp / 8.0_dp - eta / 8.0_dp) * (1.0_dp - xi) * (-eta - xi + zeta - 2.0_dp) + (1.0_dp - eta) * &
                    (1.0_dp - xi) * ( zeta / 8.0_dp + 1.0_dp / 8.0_dp), &
            (1.0_dp - eta) * (xi / 8.0_dp + 1.0_dp / 8.0_dp) * (zeta + 1.0_dp) + (1.0_dp - eta) * &
                    (xi / 8.0_dp + 1.0_dp / 8.0_dp) * ( -eta + xi + zeta - 2.0_dp), &
            (eta / 8.0_dp + 1.0_dp / 8.0_dp) * (xi + 1) * (zeta + 1.0_dp) + (eta / 8.0_dp + 1.0_dp / 8.0_dp) * &
                    (xi + 1.0_dp) * ( eta + xi + zeta - 2.0_dp), &
            (1.0_dp - xi) * (eta / 8.0_dp + 1.0_dp / 8.0_dp) * (zeta + 1.0_dp) + (1.0_dp - xi) * &
                    (eta / 8.0_dp + 1.0_dp / 8.0_dp) * ( eta - xi + zeta - 2.0_dp), &
            (0.25_dp - xi ** 2 / 4.0_dp) * (eta - 1.0_dp), &
            (eta ** 2 - 1.0_dp) * (xi / 4.0_dp + 0.25), &
            (0.25_dp - xi ** 2 / 4.0_dp) * (-eta - 1), &
            (0.25_dp - xi / 4.0_dp) * (eta ** 2 - 1.0_dp), &
            (0.25_dp - xi ** 2 / 4.0_dp) * (1.0_dp - eta), &
            (1.0_dp - eta ** 2) * (xi / 4.0_dp + 0.25), &
            (0.25_dp - xi ** 2 / 4.0_dp) * (eta + 1.0_dp), &
            (0.25_dp - xi / 4.0_dp) * (1.0_dp - eta ** 2), &
            -2.0_dp * zeta * (0.25_dp - xi / 4.0_dp) * (1.0_dp - eta), &
            -2 * zeta * (1.0_dp - eta) * (xi / 4.0_dp + 0.25), &
            -2.0_dp * zeta * (eta + 1.0_dp) * (xi / 4.0_dp + 0.25), &
            -2.0_dp * zeta * (0.25_dp - xi / 4.0_dp) * (eta + 1.0_dp) /), &
            shape=(/20,3/))
        end subroutine hexa20_jacobian

        !---------------------------------------------------------- Tet10 ---------------------------------------------
        integer function tet10_no_of_nodes()
            tet10_no_of_nodes = 10
        end function tet10_no_of_nodes

        integer function tet10_no_of_local_coordinates()
            tet10_no_of_local_coordinates = 3
        end function tet10_no_of_local_coordinates

        real(dp) function tet10_xi_lower()
            tet10_xi_lower = 0.0_dp
        end function tet10_xi_lower

        real(dp) function tet10_xi_upper()
            tet10_xi_upper = 1.0_dp
        end function tet10_xi_upper

        real(dp) function tet10_eta_lower(xi)
            real(dp), intent(in) :: xi
            tet10_eta_lower = 0.0_dp
        end function tet10_eta_lower

        real(dp) function tet10_eta_upper(xi)
            real(dp), intent(in) :: xi
            tet10_eta_upper = 1.0_dp - xi
        end function tet10_eta_upper

        real(dp) function tet10_zeta_lower(xi, eta)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            tet10_zeta_lower = 0.0_dp
        end function tet10_zeta_lower

        real(dp) function tet10_zeta_upper(xi, eta)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            tet10_zeta_upper = 1.0_dp - xi - eta
        end function tet10_zeta_upper

        subroutine tet10_evaluate(xi, eta, zeta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(in) :: zeta
            real(dp), intent(out) :: out(10)

            out = (/ 2.0_dp * eta * (eta + xi + zeta - 1.0_dp) - eta + 2.0_dp * xi * ( &
                    eta + xi + zeta - 1.0_dp) - xi + 2.0_dp * zeta * ( &
                    eta + xi + zeta - 1.0_dp) - zeta + 1.0_dp, &
                    xi * (2.0_dp * xi - 1.0_dp), &
                    eta * (2.0_dp * eta - 1.0_dp), &
                    zeta * (2.0_dp * zeta - 1.0_dp), &
                    4.0_dp * xi * (-eta - xi - zeta + 1.0_dp), &
                    4.0_dp * eta * xi, &
                    4.0_dp * eta * (-eta - xi - zeta + 1.0_dp), &
                    4.0_dp * xi * zeta, &
                    4.0_dp * eta * zeta, &
                    4.0_dp * zeta * (-eta - xi - zeta + 1.0_dp) /)
        end subroutine tet10_evaluate

        subroutine tet10_jacobian(xi, eta, zeta, out)
            real(dp), intent(in) :: xi
            real(dp), intent(in) :: eta
            real(dp), intent(in) :: zeta
            real(dp), intent(out) :: out(10,3)

            out = reshape((/ 4.0_dp * eta + 4.0_dp * xi + 4.0_dp * zeta - 3.0_dp, &
            4.0_dp * xi - 1.0_dp, &
            0.0_dp, &
            0.0_dp, &
            -4.0_dp * eta - 8.0_dp * xi - 4.0_dp * zeta + 4.0_dp, &
            4.0_dp * eta, &
            -4.0_dp * eta, &
            4.0_dp * zeta, &
            0.0_dp, &
            -4.0_dp * zeta, &
            4.0_dp * eta + 4.0_dp * xi + 4.0_dp * zeta - 3.0_dp, &
            0.0_dp, &
            4.0_dp * eta - 1.0_dp, &
            0.0_dp, &
            -4.0_dp * xi, &
            4.0_dp * xi, &
            -8.0_dp * eta - 4.0_dp * xi - 4.0_dp * zeta + 4.0_dp, &
            0.0_dp, &
            4.0_dp * zeta, &
            -4.0_dp * zeta, &
            4.0_dp * eta + 4.0_dp * xi + 4.0_dp * zeta - 3.0_dp, &
            0.0_dp, &
            0.0_dp, &
            4.0_dp * zeta - 1.0_dp, &
            -4.0_dp * xi, &
            0.0_dp, &
            -4.0_dp * eta, &
            4.0_dp * xi, &
            4.0_dp * eta, &
            -4.0_dp * eta - 4.0_dp * xi - 8.0_dp * zeta + 4.0_dp /), &
            shape=(/10,3/))
        end subroutine tet10_jacobian

end module shapefunc