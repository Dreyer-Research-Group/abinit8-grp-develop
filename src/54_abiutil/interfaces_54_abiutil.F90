!!****m* ABINIT/interfaces_54_abiutil
!! NAME
!! interfaces_54_abiutil
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/54_abiutil
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

module interfaces_54_abiutil

 implicit none

interface
 subroutine calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,nunit,ratsph,rhor,rprimd,typat,ucvol,xred,&  
  &  prtopt,cplex,intgden,dentot)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: cplex
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: ntypat
  integer,intent(in) :: nunit
  integer ,intent(in) :: prtopt
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out),optional :: dentot(nspden)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(out),optional :: intgden(nspden,natom)
  real(dp),intent(in) :: ratsph(ntypat)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine calcdensph
end interface

interface
 subroutine mkfilename(filnam,filnam_out,get,idtset,ird,jdtset_,ndtset,stringfil,stringvar,will_read)
  use defs_basis
  implicit none
  integer,intent(in) :: get
  integer,intent(in) :: idtset
  integer,intent(in) :: ird
  integer,intent(in) :: ndtset
  integer,intent(out) :: will_read
  character(len=fnlen),intent(out) :: filnam_out
  character(len=*),intent(in) :: stringfil
  character(len=*),intent(in) :: stringvar
  character(len=fnlen),intent(in) :: filnam(5)
  integer,intent(in) :: jdtset_(0:ndtset)
 end subroutine mkfilename
end interface

interface
 subroutine printmagvtk(mpi_enreg,cplex,nspden,nfft,ngfft,rhor,rprimd,fname)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: cplex
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  character(len=*),intent(in) :: fname
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rhor(cplex*nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine printmagvtk
end interface

end module interfaces_54_abiutil
!!***
