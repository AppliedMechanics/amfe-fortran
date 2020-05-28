! Copyright (c) 2020, Lehrstuhl fuer Angewandte Mechanik, Technische
! Universitaet Muenchen.
!
! Distributed under BSD-3-Clause License. See LICENSE-File for more information
!
! Author: Christian Meyer
!

module shapefunction_obj
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    use shapefunction_fun
    implicit none
    private

    public :: shapefunction1d, shapefunction2d, shapefunction3d, &
        create_shapefunction_1d, create_shapefunction_2d, create_shapefunction_3d

    integer, public, parameter :: &
            LINE2 = 12, &
            LINE3 = 13, &
            TRI3 = 23, &
            QUAD4 = 24, &
            TRI6 = 26, &
            QUAD8 = 28, &
            TET4 = 34, &
            HEXA8 = 38, &
            TET10 = 310, &
            HEXA20 = 320

    type :: shapefunction
        integer :: no_of_nodes
        integer :: no_of_local_coordinates
        procedure(), nopass, pointer :: evaluate
        procedure(), nopass, pointer :: jacobian
    end type shapefunction

    type, extends(shapefunction) :: shapefunction1d
        procedure(xi_func), nopass, pointer :: xi_lower
        procedure(xi_func), nopass, pointer :: xi_upper
    end type shapefunction1d

    type, extends(shapefunction) :: shapefunction2d
        procedure(xi_func), nopass, pointer :: xi_lower
        procedure(xi_func), nopass, pointer :: xi_upper
        procedure(eta_func), nopass, pointer :: eta_lower
        procedure(eta_func), nopass, pointer :: eta_upper
    end type shapefunction2d

    type, extends(shapefunction) :: shapefunction3d
        procedure(xi_func), nopass, pointer :: xi_lower
        procedure(xi_func), nopass, pointer :: xi_upper
        procedure(eta_func), nopass, pointer :: eta_lower
        procedure(eta_func), nopass, pointer :: eta_upper
        procedure(zeta_func), nopass, pointer :: zeta_lower
        procedure(zeta_func), nopass, pointer :: zeta_upper
    end type shapefunction3d

    interface
        function xi_func() result(out)
            import :: dp
            real(dp) :: out
        end function xi_func

        function eta_func(xi) result(out)
            import :: dp
            real(dp), intent(in) :: xi
            real(dp) :: out
        end function eta_func

        function zeta_func(xi, eta) result(out)
            import :: dp
            real(dp), intent(in) :: xi, eta
            real(dp) :: out
        end function zeta_func
    end interface

    contains
        function create_line2() result(obj)
            type(shapefunction1d) :: obj
            obj%no_of_nodes = line2_no_of_nodes()
            obj%no_of_local_coordinates = line2_no_of_local_coordinates()
            obj%evaluate => line2_evaluate
            obj%jacobian => line2_jacobian
            obj%xi_lower => line2_xi_lower
            obj%xi_upper => line2_xi_upper
        end function create_line2

        function create_line3() result(obj)
            type(shapefunction1d) :: obj
            obj%no_of_nodes = line3_no_of_nodes()
            obj%no_of_local_coordinates = line3_no_of_local_coordinates()
            obj%evaluate => line3_evaluate
            obj%jacobian => line3_jacobian
            obj%xi_lower => line3_xi_lower
            obj%xi_upper => line3_xi_upper
        end function create_line3

        function create_tri3() result(obj)
            type(shapefunction2d) :: obj
            obj%no_of_nodes = tri3_no_of_nodes()
            obj%no_of_local_coordinates = tri3_no_of_local_coordinates()
            obj%evaluate => tri3_evaluate
            obj%jacobian => tri3_jacobian
            obj%xi_lower => tri3_xi_lower
            obj%xi_upper => tri3_xi_upper
            obj%eta_lower => tri3_eta_lower
            obj%eta_upper => tri3_eta_upper
        end function create_tri3

        function create_quad4() result(obj)
            type(shapefunction2d) :: obj
            obj%no_of_nodes = quad4_no_of_nodes()
            obj%no_of_local_coordinates = quad4_no_of_local_coordinates()
            obj%evaluate => quad4_evaluate
            obj%jacobian => quad4_jacobian
            obj%xi_lower => quad4_xi_lower
            obj%xi_upper => quad4_xi_upper
            obj%eta_lower => quad4_eta_lower
            obj%eta_upper => quad4_eta_upper
        end function create_quad4

        function create_quad8() result(obj)
            type(shapefunction2d) :: obj
            obj%no_of_nodes = quad8_no_of_nodes()
            obj%no_of_local_coordinates = quad8_no_of_local_coordinates()
            obj%evaluate => quad8_evaluate
            obj%jacobian => quad8_jacobian
            obj%xi_lower => quad8_xi_lower
            obj%xi_upper => quad8_xi_upper
            obj%eta_lower => quad8_eta_lower
            obj%eta_upper => quad8_eta_upper
        end function create_quad8

        function create_tri6() result(obj)
            type(shapefunction2d) :: obj
            obj%no_of_nodes = tri6_no_of_nodes()
            obj%no_of_local_coordinates = tri6_no_of_local_coordinates()
            obj%evaluate => tri6_evaluate
            obj%jacobian => tri6_jacobian
            obj%xi_lower => tri6_xi_lower
            obj%xi_upper => tri6_xi_upper
            obj%eta_lower => tri6_eta_lower
            obj%eta_upper => tri6_eta_upper
        end function create_tri6

        function create_hexa8() result(obj)
            type(shapefunction3d) :: obj
            obj%no_of_nodes = hexa8_no_of_nodes()
            obj%no_of_local_coordinates = hexa8_no_of_local_coordinates()
            obj%evaluate => hexa8_evaluate
            obj%jacobian => hexa8_jacobian
            obj%xi_lower => hexa8_xi_lower
            obj%xi_upper => hexa8_xi_upper
            obj%eta_lower => hexa8_eta_lower
            obj%eta_upper => hexa8_eta_upper
            obj%zeta_lower => hexa8_zeta_lower
            obj%zeta_upper => hexa8_zeta_upper
        end function create_hexa8

        function create_tet4() result(obj)
            type(shapefunction3d) :: obj
            obj%no_of_nodes = tet4_no_of_nodes()
            obj%no_of_local_coordinates = tet4_no_of_local_coordinates()
            obj%evaluate => tet4_evaluate
            obj%jacobian => tet4_jacobian
            obj%xi_lower => tet4_xi_lower
            obj%xi_upper => tet4_xi_upper
            obj%eta_lower => tet4_eta_lower
            obj%eta_upper => tet4_eta_upper
        end function create_tet4

        function create_hexa20() result(obj)
            type(shapefunction3d) :: obj
            obj%no_of_nodes = hexa20_no_of_nodes()
            obj%no_of_local_coordinates = hexa20_no_of_local_coordinates()
            obj%evaluate => hexa20_evaluate
            obj%jacobian => hexa20_jacobian
            obj%xi_lower => hexa20_xi_lower
            obj%xi_upper => hexa20_xi_upper
            obj%eta_lower => hexa20_eta_lower
            obj%eta_upper => hexa20_eta_upper
        end function create_hexa20

        function create_tet10() result(obj)
            type(shapefunction3d) :: obj
            obj%no_of_nodes = tet10_no_of_nodes()
            obj%no_of_local_coordinates = tet10_no_of_local_coordinates()
            obj%evaluate => tet10_evaluate
            obj%jacobian => tet10_jacobian
            obj%xi_lower => tet10_xi_lower
            obj%xi_upper => tet10_xi_upper
            obj%eta_lower => tet10_eta_lower
            obj%eta_upper => tet10_eta_upper
        end function create_tet10

        function create_shapefunction_1d(identifier) result(obj)
            integer, intent(in) :: identifier
            type(shapefunction1d) :: obj

            select case(identifier)
                case(LINE2)
                    obj = create_line2()
                case(LINE3)
                    obj = create_line3()
                case default
                    stop 'Shapefunction not available or wrong dimension'
            end select
        end function create_shapefunction_1d

        function create_shapefunction_2d(identifier) result(obj)
            integer, intent(in) :: identifier
            type(shapefunction2d) :: obj

            select case(identifier)
                case(TRI3 )
                    obj = create_tri3()
                case(QUAD4)
                    obj = create_quad4()
                case(TRI6)
                    obj = create_tri6()
                case(QUAD8)
                    obj = create_quad8()
                case default
                    stop 'Shapefunction not available or wrong dimension'
            end select

        end function create_shapefunction_2d

        function create_shapefunction_3d(identifier) result(obj)
            integer, intent(in) :: identifier
            type(shapefunction3d) :: obj

            select case(identifier)
                case(TET4)
                    obj = create_tet4()
                case(HEXA8)
                    obj = create_hexa8()
                case(TET10)
                    obj = create_tet10()
                case(HEXA20)
                    obj = create_hexa20()
                case default
                    stop 'Shapefunction not available or wrong dimension'
            end select
        end function create_shapefunction_3d

end module shapefunction_obj
