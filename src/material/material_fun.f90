!
! Copyright (c) 2020 TECHNICAL UNIVERSITY OF MUNICH, DEPARTMENT OF MECHANICAL ENGINEERING, CHAIR OF APPLIED MECHANICS,
! BOLTZMANNSTRASSE 15, 85748 GARCHING/MUNICH, GERMANY, RIXEN@TUM.DE.
!
! Distributed under 3-Clause BSD license. See LICENSE file for more information.
!


module material_fun
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    use tools
    implicit none
    private

    public kirchhoff_sv_c_dc_3d, kirchhoff_sv_c_dc_2d, &
            mooney_rivlin_sv_c_dc_3d, mooney_rivlin_sv_c_dc_2d, &
            neo_hookean_sv_c_dc_3d, neo_hookean_sv_c_dc_2d

    contains
    
        subroutine kirchhoff_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, e_modul, nu, density, hypothesis)
            real(dp), intent(in) :: f(3,3)
            real(dp), intent(in) :: df(3,3)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(6)
            real(dp), intent(out) :: c_out(6,6)
            real(dp), intent(out) :: dc_out(6,6)
            real(dp), intent(in) :: e_modul
            real(dp), intent(in) :: nu
            real(dp), intent(in) :: density
            integer, intent(in) :: hypothesis

            real(dp) :: E(3,3)
            real(dp) :: lam, mu, Ev(6)

            lam = nu * e_modul / ((1.0_dp + nu) * (1.0_dp - 2*nu))
            mu  = e_modul / (2.0_dp*(1 + nu))

            C_out = reshape((/lam + 2*mu, lam, lam, 0.0_dp, 0.0_dp, 0.0_dp, &
                    lam, lam + 2*mu, lam, 0.0_dp, 0.0_dp, 0.0_dp, &
                    lam, lam, lam + 2*mu, 0.0_dp, 0.0_dp, 0.0_dp, &
                    0.0_dp, 0.0_dp, 0.0_dp, mu, 0.0_dp, 0.0_dp, &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, mu, 0.0_dp, &
                    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, mu/), (/6,6/))

            call green_lagrange_strain_3d(f, Ev)

            Sv_out = matmul(C_out, Ev)
            dC_out = dC_out * 0.0_dp
        end subroutine kirchhoff_sv_c_dc_3d

        subroutine kirchhoff_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, e_modul, nu, density, &
                hypothesis)
            real(dp), intent(in) :: f(2,2)
            real(dp), intent(in) :: df(2,2)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(3)
            real(dp), intent(out) :: c_out(3,3)
            real(dp), intent(out) :: dc_out(3,3)
            real(dp), intent(in) :: e_modul
            real(dp), intent(in) :: nu
            real(dp), intent(in) :: density
            integer, intent(in) :: hypothesis

            real(dp) :: lam, mu, Ev(3)
            real(dp) :: E(2,2)

            lam = nu * E_modul / ((1.0_dp + nu) * (1.0_dp - 2*nu))
            mu  = E_modul / (2.0_dp*(1 + nu))

            !     This is slow, maybe one wants to change that in the future
            if (hypothesis == 2) then
                C_out = E_modul/(1.0_dp- nu**2)*reshape((/1.0_dp, nu, 0.0_dp, &
                        nu, 1.0_dp, 0.0_dp, &
                        0.0_dp, 0.0_dp, (1.0_dp-nu)/2.0_dp/), shape(C_out))
            elseif (hypothesis == 1) then
                C_out = reshape((/lam + 2*mu, lam, 0.0_dp, &
                        lam, lam + 2*mu, 0.0_dp, &
                        0.0_dp, 0.0_dp, mu/), shape(C_out))
            else
                stop 'f_kirchhoff_sv_c_dc may only be called for plane strain/stress (hypothesis = 1 or 2)'
            end if

            dC_out = dC_out * 0.0_dp
            call green_lagrange_strain_2d(F, Ev)
            Sv_out = matmul(C_out, Ev)

        end subroutine kirchhoff_sv_c_dc_2d

        subroutine mooney_rivlin_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                a01, a10, kappa, hypothesis)
            real(dp), intent(in) :: f(3,3)
            real(dp), intent(in) :: df(3,3)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(6)
            real(dp), intent(out) :: c_out(6,6)
            real(dp), intent(out) :: dc_out(6,6)
            real(dp), intent(in) :: a01
            real(dp), intent(in) :: a10
            real(dp), intent(in) :: kappa
            integer, intent(in) :: hypothesis

            real(dp) :: C(3,3), I1E(6,1), I2E(6,1), I3E(6,1), I2EE(6,6), I3EE(6,6), &
                    J1E(6,1), J2E(6,1), J3E(6,1), EYE(3,3), J1EE(6,6), J2EE(6,6), J3EE(6,6)
            real(dp) :: C11, C22, C33, C23, C13, C12, I1, I2, I3, J3, &
                    J1I1, J1I3, J2I2, J2I3, J3I3, J1I1I3, J1I3I3, J2I2I3, J2I3I3, J3I3I3, E(3,3)

            if (hypothesis == 0) then
                EYE = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp /), shape(EYE))
                C = matmul(transpose(F),F)

                E = F
                E = 0.5 * (matmul(transpose(F),F) - EYE)

                C11 = C(1,1)
                C22 = C(2,2)
                C33 = C(3,3)
                C23 = C(2,3)
                C13 = C(1,3)
                C12 = C(1,2)
                !   invariants and reduced invariants
                I1  = C11 + C22 + C33
                I2  = C11*C22 + C11*C33 - C12**2 - C13**2 + C22*C33 - C23**2
                I3  = C11*C22*C33 - C11*C23**2 - C12**2*C33 + 2*C12*C13*C23 - C13**2*C22

                J3  = sqrt(I3)
                !   derivatives
                J1I1 = I3**(-1.0_dp/3.0_dp)
                J1I3 = -I1/(3.0_dp*I3**(4.0_dp/3.00D0))
                J2I2 = I3**(-2.0_dp/3)
                J2I3 = -2.0_dp*I2/(3.0_dp*I3**(5.0_dp/3))
                J3I3 = 1./(2.0_dp*sqrt(I3))

                I1E = reshape((/2.0_dp, 2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), shape(I1E))
                I2E = reshape((/2*C22 + 2*C33, 2*C11 + 2*C33, 2*C11 + 2*C22, -2*C23, -2*C13, -2*C12/), shape(I2E))
                I3E = reshape((/2*C22*C33 - 2*C23**2, 2*C11*C33 - 2*C13**2, 2*C11*C22 - 2*C12**2, &
                        -2*C11*C23 + 2*C12*C13, 2*C12*C23 - 2*C13*C22, -2*C12*C33 + 2*C13*C23/), shape(I3E))

                J1E = J1I1*I1E + J1I3*I3E
                J2E = J2I2*I2E + J2I3*I3E
                J3E = J3I3*I3E

                sv_out(:) = a10 * reshape(J1E, (/6/)) + a01 * reshape(J2E, (/6/)) + &
                        kappa*(J3 - 1)*reshape(J3E, (/6/))

                I2EE = reshape((/0.0_dp, 4.0_dp, 4.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
                        4.0_dp, 0.0_dp, 4.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
                        4.0_dp, 4.0_dp, 0.0_dp,  0.0_dp,  0.0_dp,  0.0_dp, &
                        0.0_dp, 0.0_dp, 0.0_dp, -2.0_dp,  0.0_dp,  0.0_dp, &
                        0.0_dp, 0.0_dp, 0.0_dp,  0.0_dp, -2.0_dp,  0.0_dp, &
                        0.0_dp, 0.0_dp, 0.0_dp,  0.0_dp,  0.0_dp, -2.0_dp/), shape(I2EE))

                I3EE = reshape((/    0.0_dp,  4.0_dp*C33,  4.0_dp*C22, -4.0_dp*C23,      0.0_dp,      0.0_dp, &
                        4.0_dp*C33,      0.0_dp,  4.0_dp*C11,      0.0_dp, -4.0_dp*C13,      0.0_dp, &
                        4.0_dp*C22,  4.0_dp*C11,      0.0_dp,      0.0_dp,      0.0_dp, -4.0_dp*C12, &
                        -4.0_dp*C23,      0.0_dp,      0.0_dp, -2.0_dp*C11,  2.0_dp*C12,  2.0_dp*C13, &
                        0.0_dp, -4.0_dp*C13,      0.0_dp,  2.0_dp*C12, -2.0_dp*C22,  2.0_dp*C23, &
                        0.0_dp,      0.0_dp, -4.0_dp*C12,  2.0_dp*C13,  2.0_dp*C23, -2.0_dp*C33/), shape(I3EE))

                !   second derivatives
                J1I1I3 = -1.0_dp/(3.0_dp*I3**(4.0_dp/3))
                J1I3I3 = 4*I1/(9*I3**(7.0_dp/3))
                J2I2I3 = -2.0_dp/(3*I3**(5.0_dp/3))
                J2I3I3 = 10*I2/(9*I3**(8.0_dp/3))
                J3I3I3 = -1.0_dp/(4*I3**(3.0_dp/2))

                J1EE = J1I1I3*(matmul(I1E, transpose(I3E)) + matmul(I3E, transpose(I1E))) &
                        + J1I3I3*matmul(I3E, transpose(I3E)) + J1I3*I3EE
                J2EE = J2I2I3*(matmul(I2E, transpose(I3E)) + matmul(I3E, transpose(I2E))) &
                        + J2I3I3*matmul(I3E, transpose(I3E)) + J2I2*I2EE + J2I3*I3EE
                J3EE = J3I3I3*(matmul(I3E, transpose(I3E))) + J3I3*I3EE

                c_out(:,:) = a10 * J1EE + a01*J2EE + kappa * (matmul(J3E, transpose(J3E))) + &
                        kappa *(J3-1.0)*J3EE
                dc_out(:,:) = 0.0_dp
            else
                stop 'mooney_rivlin_sv_c_dc_3d may only be called for tridimensional hypothesis (hypothesis = 0)'
            end if
        end subroutine mooney_rivlin_sv_c_dc_3d

        subroutine mooney_rivlin_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, a01, a10, kappa, &
                hypothesis)
            real(dp), intent(in) :: f(2,2)
            real(dp), intent(in) :: df(2,2)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(3)
            real(dp), intent(out) :: c_out(3,3)
            real(dp), intent(out) :: dc_out(3,3)
            real(dp), intent(in) :: a01
            real(dp), intent(in) :: a10
            real(dp), intent(in) :: kappa
            integer, intent(in) :: hypothesis

            stop 'mooney rivlin is not implemented for 2d'
        end subroutine mooney_rivlin_sv_c_dc_2d

        subroutine neo_hookean_sv_c_dc_3d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, &
                kappa, mu, hypothesis)
            real(dp), intent(in) :: f(3,3)
            real(dp), intent(in) :: df(3,3)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(6)
            real(dp), intent(out) :: c_out(6,6)
            real(dp), intent(out) :: dc_out(6,6)
            real(dp), intent(in) :: kappa
            real(dp), intent(in) :: mu
            integer, intent(in) :: hypothesis

            real(dp) :: C(3,3), I1E(6,1), I3E(6,1), I3EE(6,6), J1E(6,1), J3E(6,1), EYE(3,3), J1EE(6,6), J3EE(6,6), E(3,3)
            real(dp) :: C11, C22, C33, C23, C13, C12, I1, I3, J3, J1I1, J1I3, J3I3, J1I1I3, J1I3I3, J3I3I3

            if (hypothesis == 0) then
                EYE = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), shape(EYE))
                C = matmul(transpose(F),F)

                E = 0.5_dp * (matmul(transpose(F),F) - EYE)

                C11 = C(1,1)
                C22 = C(2,2)
                C33 = C(3,3)
                C23 = C(2,3)
                C13 = C(1,3)
                C12 = C(1,2)
                !   invariants and reduced invariants
                I1  = C11 + C22 + C33
                I3  = C11*C22*C33 - C11*C23**2 - C12**2*C33 + 2*C12*C13*C23 - C13**2*C22

                J3  = sqrt(I3)
                !   derivatives
                J1I1 = I3**(-1.0/3.0_dp)
                J1I3 = -I1/(3*I3**(4.0_dp/3.0_dp))
                J3I3 = 1.0_dp/(2.0_dp*sqrt(I3))

                I1E = reshape((/2.0_dp, 2.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), shape(I1E))
                I3E = reshape((/2*C22*C33 - 2*C23**2, 2*C11*C33 - 2*C13**2, 2*C11*C22 - 2*C12**2, &
                        -2*C11*C23 + 2*C12*C13, 2*C12*C23 - 2*C13*C22, -2*C12*C33 + 2*C13*C23/), shape(I3E))

                J1E = J1I1*I1E + J1I3*I3E
                J3E = J3I3*I3E

                sv_out(:) = mu /2.*reshape(J1E, (/6/)) + kappa * (J3 - 1)*reshape(J3E, (/6/))

                I3EE = reshape((/    0.0_dp,  4.0_dp*C33,  4.0_dp*C22, -4.0_dp*C23,      0.0_dp,      0.0_dp, &
                        4.0_dp*C33,      0.0_dp,  4.0_dp*C11,      0.0_dp, -4.0_dp*C13,      0.0_dp, &
                        4.0_dp*C22,  4.0_dp*C11,      0.0_dp,      0.0_dp,      0.0_dp, -4.0_dp*C12, &
                        -4.0_dp*C23,      0.0_dp,      0.0_dp, -2.0_dp*C11,  2.0_dp*C12,  2.0_dp*C13, &
                        0.0_dp, -4.0_dp*C13,      0.0_dp,  2.0_dp*C12, -2.0_dp*C22,  2.0_dp*C23, &
                        0.0_dp,      0.0_dp, -4.0_dp*C12,  2.0_dp*C13,  2.0_dp*C23, -2.0_dp*C33/), shape(I3EE))

                !     second derivatives
                J1I1I3 = -1.0_dp/(3.0_dp*I3**(4.0_dp/3))
                J1I3I3 = 4*I1/(9*I3**(7.0_dp/3.0_dp))
                J3I3I3 = -1.0_dp/(4*I3**(3.0_dp/2.0_dp))

                J1EE = J1I1I3*(matmul(I1E, transpose(I3E)) + matmul(I3E, transpose(I1E))) &
                        + J1I3I3*matmul(I3E, transpose(I3E)) + J1I3*I3EE
                J3EE = J3I3I3*(matmul(I3E, transpose(I3E))) + J3I3*I3EE

                c_out(:,:) = mu /2*J1EE + kappa * (matmul(J3E, transpose(J3E))) + kappa * (J3-1.0_dp)*J3EE
                dc_out(:,:) = 0.0_dp
            else
                stop 'neo_hookean_sv_c_dc_3d is only allowed for tridimensional hypothesis (hypothesis = 0)'

            end if
        end subroutine neo_hookean_sv_c_dc_3d

        subroutine neo_hookean_sv_c_dc_2d(f, df, t, s0, t0, s1, sv_out, c_out, dc_out, kappa, mu, hypothesis)
            real(dp), intent(in) :: f(2,2)
            real(dp), intent(in) :: df(2,2)
            real(dp), intent(in) :: t
            real(dp), intent(in) :: s0(:)
            real(dp), intent(in) :: t0
            real(dp), intent(out) :: s1(:)
            real(dp), intent(out) :: sv_out(3)
            real(dp), intent(out) :: c_out(3,3)
            real(dp), intent(out) :: dc_out(3,3)
            real(dp), intent(in) :: kappa
            real(dp), intent(in) :: mu
            integer, intent(in) :: hypothesis

            real(dp) :: C(2,2), I1E(3,1), I3E(3,1), I3EE(3,3), J1E(3,1), J3E(3,1), EYE(2,2), J1EE(3,3), J3EE(3,3),E(2,2)
            real(dp) :: C11, C22, C33, C12, I1, I3, J3, J1I1, J1I3, J3I3, J1I1I3, J1I3I3, J3I3I3

            if (hypothesis == 1) then

                EYE = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), shape(EYE))

                C = matmul(transpose(F),F)
                E = 0.5_dp * (C - EYE)

                C11 = C(1,1)
                C22 = C(2,2)
                C12 = C(1,2)
                C33 = 1.0_dp

                !   invariants and reduced invariants
                I1  = C11 + C22 + C33
                I3  = C11*C22 - C12**2
                J3  = sqrt(I3)

                !   derivatives
                J1I1 = I3**(-1.0/3.0_dp)
                J1I3 = -I1/(3*I3**(4.00D0/3.00D0))
                J3I3 = 1.0_dp/(2.0_dp*sqrt(I3))

                I1E = reshape((/2.0_dp, 2.0_dp, 0.0_dp /), shape(I1E))
                I3E = reshape((/2*C22*C33, 2*C11*C33, -2*C12*C33/), shape(I3E))

                J1E = J1I1*I1E + J1I3*I3E
                J3E = J3I3*I3E

                sv_out(:) = mu /2.0_dp*reshape(J1E, (/3/)) + kappa*(J3 - 1.0_dp)*reshape(J3E, (/3/))

                I3EE = reshape((/ 0.0_dp, 4.0_dp, 0.0_dp, &
                        4.0_dp, 0.0_dp, 0.0_dp, &
                        0.0_dp, 0.0_dp,-2.0_dp /), shape(I3EE))

                !   second derivatives
                J1I1I3 = -1.0_dp/(3.0_dp*I3**(4.0_dp/3))
                J1I3I3 = 4.0_dp*I1/(9.0_dp*I3**(7.0_dp/3))
                J3I3I3 = -1.0_dp/(4.0_dp*I3**(3.0_dp/2))

                J1EE = J1I1I3*(matmul(I1E, transpose(I3E)) + matmul(I3E, transpose(I1E))) &
                        + J1I3I3*matmul(I3E, transpose(I3E)) + J1I3*I3EE
                J3EE = J3I3I3*(matmul(I3E, transpose(I3E))) + J3I3*I3EE

                c_out(:,:) = mu /2.0_dp*J1EE + kappa * (matmul(J3E, transpose(J3E))) + kappa*(J3-1.0_dp)*J3EE
                dc_out(:,:) = 0.0_dp
            else
                stop 'neo_hookean_2d is only available for plane strain (hypothesis = 1)'
            end if

        end subroutine neo_hookean_sv_c_dc_2d

end module material_fun