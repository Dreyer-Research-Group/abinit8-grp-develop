!!****m* ABINIT/interfaces_28_numeric_noabirule
!! NAME
!! interfaces_28_numeric_noabirule
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/28_numeric_noabirule
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_28_numeric_noabirule

 implicit none

interface
 subroutine GAMMA_FUNCTION(X,GA)
  use defs_basis
  implicit none
  real(dp),intent(out) :: ga
  real(dp),intent(in) :: x
 end subroutine GAMMA_FUNCTION
end interface


end module interfaces_28_numeric_noabirule
!!***
