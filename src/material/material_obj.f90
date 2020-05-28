! Copyright (c) 2020, Lehrstuhl fuer Angewandte Mechanik, Technische
! Universitaet Muenchen.
!
! Distributed under BSD-3-Clause License. See LICENSE-File for more information
!
! Author: Christian Meyer
!

module material_obj
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    use material_fun
    implicit none
    private

    public kirchhoffmaterial

    integer, public, parameter :: HypothesisTridimensional = 0
    integer, public, parameter :: HypothesisPlanestrain = 1
    integer, public, parameter :: HypothesisPlanestress = 2


    type, abstract :: material
        contains
            procedure(sv_c_dc_3d_interface), deferred :: sv_c_dc_3d
            procedure(sv_c_dc_2d_interface), deferred :: sv_c_dc_2d
    end type material

    interface
        subroutine sv_c_dc_3d_interface(this, f, df, t, s0, t0, s1, sv_out, c_out, dc_out)
            import material, dp
            class(material) :: this
            real(dp), intent(in) :: f(3,3)
            real(dp), intent(in) :: df(3,3)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(6)
            real(dp), intent(out) :: c_out(6,6)
            real(dp), intent(out) :: dc_out(6,6)
        end subroutine sv_c_dc_3d_interface
    end interface

    interface
        subroutine sv_c_dc_2d_interface(this, f, df, t, s0, t0, s1, sv_out, c_out, dc_out)
            import material, dp
            class(material) :: this
            real(dp), intent(in) :: f(2,2)
            real(dp), intent(in) :: df(2,2)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(3)
            real(dp), intent(out) :: c_out(3,3)
            real(dp), intent(out) :: dc_out(3,3)
        end subroutine sv_c_dc_2d_interface
    end interface

    type, extends(material) :: kirchhoffmaterial
        real(dp) :: youngs_modulus
        real(dp) :: poissons_ratio
        real(dp) :: density
        integer :: hypothesis
        contains
            procedure :: sv_c_dc_3d => kirchhoff_obj_sv_c_dc_3d
            procedure :: sv_c_dc_2d => kirchhoff_obj_sv_c_dc_2d
    end type kirchhoffmaterial

    contains
        subroutine kirchhoff_obj_sv_c_dc_3d(this, f, df, t, s0, t0, s1, sv_out, c_out, dc_out)
            class(kirchhoffmaterial) :: this
            real(dp), intent(in) :: f(3,3)
            real(dp), intent(in) :: df(3,3)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(6)
            real(dp), intent(out) :: c_out(6,6)
            real(dp), intent(out) :: dc_out(6,6)
            call kirchhoff_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                    this%youngs_modulus, this%poissons_ratio, this%density, this%hypothesis)
        end subroutine kirchhoff_obj_sv_c_dc_3d

        subroutine kirchhoff_obj_sv_c_dc_2d(this, f, df, t, s0, t0, s1, sv_out, c_out, dc_out)
            class(kirchhoffmaterial) :: this
            real(dp), intent(in) :: f(2,2)
            real(dp), intent(in) :: df(2,2)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(3)
            real(dp), intent(out) :: c_out(3,3)
            real(dp), intent(out) :: dc_out(3,3)

            call kirchhoff_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                    this%youngs_modulus, this%poissons_ratio, this%density, this%hypothesis)
        end subroutine kirchhoff_obj_sv_c_dc_2d

end module material_obj
