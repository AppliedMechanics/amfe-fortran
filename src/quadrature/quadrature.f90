! Copyright (c) 2020, Lehrstuhl fuer Angewandte Mechanik, Technische
! Universitaet Muenchen.
!
! Distributed under BSD-3-Clause License. See LICENSE-File for more information
!
! Author: Christian Meyer
!
module quadrature
    use, intrinsic:: iso_fortran_env, only: dp=>real64
    implicit none
    
    private
    public quadrature_1d, quadrature_2d, quadrature_3d
    
    contains
        
        subroutine quadrature_1d(f_1d, out, evaluationpoints, weights, out_no_of_rows, out_no_of_cols, no_of_gps)
            integer, intent(in) :: out_no_of_rows
            integer, intent(in) :: out_no_of_cols
            integer, intent(in) :: no_of_gps
            interface
                subroutine integrand_1d(xi_, out_)
                    import :: out_no_of_rows, out_no_of_cols, dp
                    real(dp), intent(in) :: xi_
                    real(dp), intent(out) :: out_(out_no_of_rows, out_no_of_cols)
                end subroutine integrand_1d
            end interface
            procedure(integrand_1d) :: f_1d
            real(dp), intent(out) :: out(out_no_of_rows, out_no_of_cols)
            real(dp), intent(in) :: evaluationpoints(1,no_of_gps)
            real(dp), intent(in) :: weights(no_of_gps)

            real(dp) :: out_local(out_no_of_rows, out_no_of_cols)
            integer :: g
        
            out(:,:) = 0.0_dp
            do g=1,no_of_gps
                call f_1d(evaluationpoints(1,g),out_local)
                out = out + weights(g) * out_local
            end do
        
        end subroutine quadrature_1d
        
        subroutine quadrature_2d(f_2d, out, evaluationpoints, weights, out_no_of_rows, out_no_of_cols, no_of_gps)
            integer, intent(in) :: out_no_of_rows
            integer, intent(in) :: out_no_of_cols
            integer, intent(in) :: no_of_gps
            interface
                subroutine integrand_2d(xi_, eta_, out_)
                    import :: out_no_of_rows, out_no_of_cols, dp
                    real(dp), intent(in) :: xi_
                    real(dp), intent(in) :: eta_
                    real(dp), intent(out) :: out_(out_no_of_rows, out_no_of_cols)
                end subroutine integrand_2d
            end interface
            procedure(integrand_2d) :: f_2d
            real(dp), intent(inout) :: out(out_no_of_rows, out_no_of_cols)
            real(dp), intent(in) :: evaluationpoints(2,no_of_gps)
            real(dp), intent(in) :: weights(no_of_gps)
        
            real(dp) :: out_local(out_no_of_rows, out_no_of_cols)
            integer :: g
        
            out(:,:) = 0.0_dp
            do g=1,no_of_gps
                call f_2d(evaluationpoints(1,g), evaluationpoints(2,g), out_local)
                out = out + weights(g) * out_local
            end do
        
        end subroutine quadrature_2d
        
        subroutine quadrature_3d(f_3d, out, evaluationpoints, weights, out_no_of_rows, out_no_of_cols, no_of_gps)
            integer, intent(in) :: out_no_of_rows
            integer, intent(in) :: out_no_of_cols
            integer, intent(in) :: no_of_gps
            interface
                subroutine integrand_3d(xi_, eta_, zeta_, out_)
                    import :: out_no_of_rows, out_no_of_cols, dp
                    real(dp), intent(in) :: xi_
                    real(dp), intent(in) :: eta_
                    real(dp), intent(in) :: zeta_
                    real(dp), intent(out) :: out_(out_no_of_rows, out_no_of_cols)
                end subroutine integrand_3d
            end interface
            procedure(integrand_3d) :: f_3d
            real(dp), intent(inout) :: out(out_no_of_rows, out_no_of_cols)
            real(dp), intent(in) :: evaluationpoints(3,no_of_gps)
            real(dp), intent(in) :: weights(no_of_gps)
        
            real(dp) :: out_local(out_no_of_rows, out_no_of_cols)
            integer :: g
        
            out(:,:) = 0.0_dp
            do g=1,no_of_gps
                call f_3d(evaluationpoints(1,g), evaluationpoints(1,g), evaluationpoints(1,g), &
                        out_local)
                out = out + weights(g) * out_local
            end do
        
        end subroutine quadrature_3d

end module quadrature
