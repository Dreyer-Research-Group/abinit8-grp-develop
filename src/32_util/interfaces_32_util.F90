!!****m* ABINIT/interfaces_32_util
!! NAME
!! interfaces_32_util
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/32_util
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

module interfaces_32_util

 implicit none

interface
 subroutine appdig(integ,string,strinn)
  implicit none
  integer,intent(in) :: integ
  character(len=*),intent(in) :: string
  character(len=*),intent(out) :: strinn
 end subroutine appdig
end interface

interface
 subroutine fappnd(filapp,filnam,iapp,&  
  &  suff) ! optional argument
  use defs_basis
  implicit none
  integer,intent(in) :: iapp
  character(len=fnlen),intent(out) :: filapp
  character(len=fnlen),intent(in) :: filnam
  character(len=3),optional,intent(in) :: suff
 end subroutine fappnd
end interface

interface
 subroutine isfile(filnam,status)
  use defs_basis
  implicit none
  character(len=fnlen),intent(inout) :: filnam
  character(len=3),intent(in) :: status
 end subroutine isfile
end interface

interface
 subroutine littlegroup_q(nsym,qpt,symq,symrec,symafm,timrev,prtvol,use_sym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(in),optional :: prtvol
  integer,intent(out) :: timrev
  integer,intent(in),optional :: use_sym
  real(dp),intent(in) :: qpt(3)
  integer,intent(in) :: symafm(nsym)
  integer,intent(out) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine littlegroup_q
end interface

interface
 subroutine test_unused_arg(used_arg,unused_arg)
  implicit none
  integer, intent(in) :: unused_arg
  integer, intent(inout) :: used_arg
 end subroutine test_unused_arg
end interface

interface
 subroutine test_same_actual_arg(dummy_out1,dummy_out2,used_arg)
  implicit none
  integer, intent(out) :: dummy_out1
  integer, intent(out) :: dummy_out2
  integer, intent(in) :: used_arg
 end subroutine test_same_actual_arg
end interface

interface
 subroutine mati3det(mm,det)
  implicit none
  integer,intent(out) :: det
  integer,intent(in) :: mm(3,3)
 end subroutine mati3det
end interface

interface
 subroutine mati3inv(mm,mit)
  implicit none
  integer,intent(out) :: mit(3,3)
  integer,intent(in) :: mm(3,3)
 end subroutine mati3inv
end interface

interface
 subroutine matr3inv(aa,ait)
  use defs_basis
  implicit none
  real(dp),intent(in) :: aa(3,3)
  real(dp),intent(out) :: ait(3,3)
 end subroutine matr3inv
end interface

interface
 subroutine mknormpath(nbounds,bounds,gmet,ndiv_small,ndiv,npt_tot,path)
  use defs_basis
  implicit none
  integer,intent(in) :: nbounds
  integer,intent(in) :: ndiv_small
  integer,intent(inout) :: npt_tot
  real(dp),intent(in) :: bounds(3,nbounds)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(inout) :: ndiv(nbounds-1)
  real(dp),intent(out),optional :: path(3,npt_tot)
 end subroutine mknormpath
end interface

interface
 function proc_distrb_cycle(distrb,ikpt,iband1,iband2,isppol,me) 
  implicit none
  integer,intent(in) :: iband1
  integer,intent(in) :: iband2
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: me
  logical :: proc_distrb_cycle
  integer,allocatable,intent(in) :: distrb(:,:,:)
 end function proc_distrb_cycle
end interface

interface
 subroutine radsintr(funr,funq,mqgrid,mrgrid,qgrid,rgrid,yq1,yqn)
  use defs_basis
  implicit none
  integer , intent(in) :: mqgrid
  integer , intent(in) :: mrgrid
  real(dp), intent(out) :: yq1
  real(dp), intent(out) :: yqn
  real(dp), intent(out) :: funq(mqgrid)
  real(dp), intent(in) :: funr(mrgrid)
  real(dp), intent(in) :: qgrid(mqgrid)
  real(dp), intent(in) :: rgrid(mrgrid)
 end subroutine radsintr
end interface

interface
 function radsmear(r, rsph, rsm)
  use defs_basis
  implicit none
  real(dp), intent(in) :: r
  real(dp) :: radsmear
  real(dp), intent(in) :: rsm
  real(dp), intent(in) :: rsph
 end function radsmear
end interface

interface
 subroutine smatrix(cg,cgq,cg1_k,ddkflag,dtm_k,icg,icg1,itrs,job,maxbd,&  
  &  mcg_k,mcg_q,mcg1_k,minbd,mpw,mband_occ,nband_occ,npw_k1,npw_k2,nspinor,&  
  &  pwind_k,pwnsfac_k,sflag_k,shiftbd,smat_inv,smat_k,smat_k_paw,usepaw)
  use defs_basis
  implicit none
  integer,intent(in) :: ddkflag
  integer,intent(in) :: icg
  integer,intent(in) :: icg1
  integer,intent(in) :: itrs
  integer,intent(in) :: job
  integer,intent(in) :: maxbd
  integer,intent(in) :: mband_occ
  integer,intent(in) :: mcg1_k
  integer,intent(in) :: mcg_k
  integer,intent(in) :: mcg_q
  integer,intent(in) :: minbd
  integer,intent(in) :: mpw
  integer,intent(in) :: nband_occ
  integer,intent(in) :: npw_k1
  integer,intent(in) :: npw_k2
  integer,intent(in) :: nspinor
  integer,intent(in) :: shiftbd
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: cg(2,mcg_k)
  real(dp),intent(out) :: cg1_k(2,mcg1_k)
  real(dp),intent(in) :: cgq(2,mcg_q)
  real(dp),intent(out) :: dtm_k(2)
  integer,intent(in) :: pwind_k(mpw)
  real(dp),intent(in) :: pwnsfac_k(4,mpw)
  integer,intent(inout) :: sflag_k(mband_occ)
  real(dp),intent(out) :: smat_inv(2,mband_occ,mband_occ)
  real(dp),intent(inout) :: smat_k(2,mband_occ,mband_occ)
  real(dp),intent(in) :: smat_k_paw(2,usepaw*mband_occ,usepaw*mband_occ)
 end subroutine smatrix
end interface

interface
 subroutine status(counter,filstat,istat,level,routine)
  implicit none
  integer,intent(in) :: counter
  integer,intent(in) :: istat
  integer,intent(in) :: level
  character(len=*),intent(in) :: filstat
  character(len=*),intent(in) :: routine
 end subroutine status
end interface


end module interfaces_32_util
!!***
