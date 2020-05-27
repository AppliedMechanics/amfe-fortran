!
! Copyright (c) 2020 TECHNICAL UNIVERSITY OF MUNICH, DEPARTMENT OF MECHANICAL ENGINEERING, CHAIR OF APPLIED MECHANICS,
! BOLTZMANNSTRASSE 15, 85748 GARCHING/MUNICH, GERMANY, RIXEN@TUM.DE.
!
! Distributed under 3-Clause BSD license. See LICENSE file for more information.
!

module fields
    use, intrinsic:: iso_fortran_env, only: dp=>real64

    implicit none
    private

    public UX, UY, UZ

    integer :: UX = 11
    integer :: UY = 12
    integer :: UZ = 13
end module fields
