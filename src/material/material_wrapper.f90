module material_wrapper

    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use material_fun

    implicit none

    private

    public c_kirchhoff_sv_c_dc_3d, c_kirchhoff_sv_c_dc_2d, &
            c_mooney_rivlin_sv_c_dc_3d, c_mooney_rivlin_sv_c_dc_2d, &
            c_neo_hookean_sv_c_dc_3d, c_neo_hookean_sv_c_dc_2d

    contains

        subroutine c_kirchhoff_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                e_modul, nu, density, hypothesis) bind(c, name='c_kirchhoff_sv_c_dc_3d')
            real(c_double), intent(in) :: f(3,3)
            real(c_double), intent(in) :: df(3,3)
            real(c_double), intent(in) :: t
            real(c_double), intent(in) :: s0(:)
            real(c_double), intent(in) :: t0
            real(c_double), intent(out) :: s1(:)
            real(c_double), intent(out) :: sv_out(6)
            real(c_double), intent(out) :: c_out(6,6)
            real(c_double), intent(out) :: dc_out(6,6)
            real(c_double), intent(in) :: e_modul
            real(c_double), intent(in) :: nu
            real(c_double), intent(in) :: density
            integer(c_int), intent(in) :: hypothesis

            call kirchhoff_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                e_modul, nu, density, hypothesis)
        end subroutine c_kirchhoff_sv_c_dc_3d

        subroutine c_kirchhoff_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, e_modul, nu, density, &
                hypothesis) bind(c, name='c_kirchhoff_sv_c_dc_2d')
            real(c_double), intent(in) :: f(2,2)
            real(c_double), intent(in) :: df(2,2)
            real(c_double), intent(in) :: t
            real(c_double), intent(in) :: s0(:)
            real(c_double), intent(in) :: t0
            real(c_double), intent(out) :: s1(:)
            real(c_double), intent(out) :: sv_out(3)
            real(c_double), intent(out) :: c_out(3,3)
            real(c_double), intent(out) :: dc_out(3,3)
            real(c_double), intent(in) :: e_modul
            real(c_double), intent(in) :: nu
            real(c_double), intent(in) :: density
            integer(c_int), intent(in) :: hypothesis

            call kirchhoff_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                e_modul, nu, density, hypothesis)

        end subroutine c_kirchhoff_sv_c_dc_2d

        subroutine c_mooney_rivlin_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                a01, a10, kappa, hypothesis) bind(c, name='c_mooney_rivlin_sv_c_dc_3d')
            real(c_double), intent(in) :: f(3,3)
            real(c_double), intent(in) :: df(3,3)
            real(c_double), intent(in) :: t
            real(c_double), intent(in) :: s0(:)
            real(c_double), intent(in) :: t0
            real(c_double), intent(out) :: s1(:)
            real(c_double), intent(out) :: sv_out(6)
            real(c_double), intent(out) :: c_out(6,6)
            real(c_double), intent(out) :: dc_out(6,6)
            real(c_double), intent(in) :: a01
            real(c_double), intent(in) :: a10
            real(c_double), intent(in) :: kappa
            integer(c_int), intent(in) :: hypothesis

            call mooney_rivlin_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                a01, a10, kappa, hypothesis)
        end subroutine c_mooney_rivlin_sv_c_dc_3d

        subroutine c_mooney_rivlin_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, a01, a10, kappa, &
                hypothesis) bind(c, name='c_mooney_rivlin_sv_c_dc_2d')
            real(c_double), intent(in) :: f(2,2)
            real(c_double), intent(in) :: df(2,2)
            real(c_double), intent(in) :: t
            real(c_double), intent(in) :: s0(:)
            real(c_double), intent(in) :: t0
            real(c_double), intent(out) :: s1(:)
            real(c_double), intent(out) :: sv_out(3)
            real(c_double), intent(out) :: c_out(3,3)
            real(c_double), intent(out) :: dc_out(3,3)
            real(c_double), intent(in) :: a01
            real(c_double), intent(in) :: a10
            real(c_double), intent(in) :: kappa
            integer(c_int), intent(in) :: hypothesis

            call mooney_rivlin_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                a01, a10, kappa, hypothesis)
        end subroutine c_mooney_rivlin_sv_c_dc_2d

        subroutine c_neo_hookean_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                kappa, mu, hypothesis) bind(c, name='c_neo_hookean_sv_c_dc_3d')
            real(c_double), intent(in) :: f(3,3)
            real(c_double), intent(in) :: df(3,3)
            real(c_double), intent(in) :: t
            real(c_double), intent(in) :: s0(:)
            real(c_double), intent(in) :: t0
            real(c_double), intent(out) :: s1(:)
            real(c_double), intent(out) :: sv_out(6)
            real(c_double), intent(out) :: c_out(6,6)
            real(c_double), intent(out) :: dc_out(6,6)
            real(c_double), intent(in) :: kappa
            real(c_double), intent(in) :: mu
            integer(c_int), intent(in) :: hypothesis

            call neo_hookean_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                kappa, mu, hypothesis)
        end subroutine c_neo_hookean_sv_c_dc_3d

        subroutine c_neo_hookean_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, kappa, mu, &
                hypothesis) bind(c, name='c_neo_hookean_sv_c_dc_2d')
            real(c_double), intent(in) :: f(2,2)
            real(c_double), intent(in) :: df(2,2)
            real(c_double), intent(in) :: t
            real(c_double), intent(in) :: s0(:)
            real(c_double), intent(in) :: t0
            real(c_double), intent(out) :: s1(:)
            real(c_double), intent(out) :: sv_out(3)
            real(c_double), intent(out) :: c_out(3,3)
            real(c_double), intent(out) :: dc_out(3,3)
            real(c_double), intent(in) :: kappa
            real(c_double), intent(in) :: mu
            integer(c_int), intent(in) :: hypothesis

            call neo_hookean_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, kappa, mu, &
                hypothesis)

        end subroutine c_neo_hookean_sv_c_dc_2d

end module material_wrapper
