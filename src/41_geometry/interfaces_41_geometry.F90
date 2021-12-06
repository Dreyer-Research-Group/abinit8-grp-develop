!!****m* ABINIT/interfaces_41_geometry
!! NAME
!! interfaces_41_geometry
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/41_geometry
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

module interfaces_41_geometry

 implicit none

interface
 subroutine bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(inout) :: nogen
  integer,intent(in) :: nsym
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine bldgrp
end interface

interface
 subroutine chkgrp(nsym,symafm,symrel,ierr)
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: nsym
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine chkgrp
end interface

interface
 subroutine sg_multable(nsym,symafm,symrel,tnons,tnons_tol,ierr,multable,toinv)
  use defs_basis
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: nsym
  real(dp),intent(in) :: tnons_tol
  integer,optional,intent(out) :: multable(4,nsym,nsym)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,optional,intent(out) :: toinv(4,nsym)
 end subroutine sg_multable
end interface

interface
 subroutine chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel)
  use defs_basis
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(in) :: nsym
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine chkorthsy
end interface

interface
 subroutine chkprimit(chkprim,multi,nsym,symafm,symrel)
  implicit none
  integer,intent(in) :: chkprim
  integer,intent(out) :: multi
  integer,intent(in) :: nsym
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine chkprimit
end interface

interface
 subroutine fillcell(natom,natrd,nsym,nucdipmom,spinat,symafm,symrel,tnons,tolsym,typat,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: natrd
  integer,intent(in) :: nsym
  real(dp),intent(in) :: tolsym
  real(dp),intent(inout) :: nucdipmom(3,natom)
  real(dp),intent(inout) :: spinat(3,natom)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(inout) :: typat(natom)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine fillcell
end interface

interface
 subroutine gensymshub(genafm,spgroup,spgroupma,shubnikov)
  use defs_basis
  implicit none
  integer,intent(out) :: shubnikov
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  real(dp),intent(out) :: genafm(3)
 end subroutine gensymshub
end interface

interface
 subroutine gensymshub4(genafm,msym,nsym,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(inout) :: nsym
  real(dp),intent(in) :: genafm(3)
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine gensymshub4
end interface

interface
 subroutine gensymspgr(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,&  
  &  spgroup,spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(inout) :: brvltt
  integer,intent(in) :: msym
  integer,intent(out) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine gensymspgr
end interface

interface
 subroutine getptgroupma(ptgroup,ptgroupha,ptgroupma)
  implicit none
  integer,intent(out) :: ptgroupma
  character(len=5),intent(in) :: ptgroup
  character(len=5),intent(in) :: ptgroupha
 end subroutine getptgroupma
end interface

interface
 subroutine getspinrot(rprimd,spinrot,symrel_conv)
  use defs_basis
  implicit none
  integer,intent(in) :: symrel_conv(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: spinrot(4)
 end subroutine getspinrot
end interface

interface
 subroutine holocell(cell_base,enforce,foundc,iholohedry,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: enforce
  integer,intent(out) :: foundc
  integer,intent(in) :: iholohedry
  real(dp),intent(in) :: tolsym
  real(dp),intent(inout) :: cell_base(3,3)
 end subroutine holocell
end interface

interface
 subroutine irreducible_set_pert(indsym,mpert,natom,nsym,pertsy,rfdir,rfpert,symq,symrec,symrel)
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: rfdir(3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(out) :: pertsy(3,mpert)
  integer,intent(in) :: rfpert(mpert)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine irreducible_set_pert
end interface

interface
 subroutine littlegroup_pert(gprimd,idir,indsym,iout,ipert,natom,nsym,nsym1,&  
  &  rfmeth,symafm,symaf1,symq,symrec,symrel,symrl1,syuse,tnons,tnons1,&  
  &  unit) ! Optional
  use defs_basis
  implicit none
  integer,intent(in) :: idir
  integer,intent(in) :: iout
  integer,intent(in) :: ipert
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(out) :: nsym1
  integer,intent(in) :: rfmeth
  integer,intent(in) :: syuse
  integer,intent(in),optional :: unit
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(out) :: symaf1(nsym)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  integer,intent(out) :: symrl1(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(out) :: tnons1(3,nsym)
 end subroutine littlegroup_pert
end interface

interface
 subroutine mksupercell(xred_org,magv_org,rprimd_org,nat_org,nat_sc,xred_sc,magv_sc,rprimd_sc,ext,prtvol) 
  use defs_basis
  implicit none
  integer,intent(in) :: nat_org
  integer,intent(in) :: nat_sc
  integer,intent(in),optional :: prtvol
  integer,intent(in) :: ext(3)
  integer,intent(in),optional :: magv_org(nat_org)
  real(dp),intent(out) :: magv_sc(nat_sc)
  real(dp),intent(in) :: rprimd_org(3,3)
  real(dp),intent(out) :: rprimd_sc(3,3)
  real(dp),intent(in) :: xred_org(3,nat_org)
  real(dp),intent(out) :: xred_sc(3,nat_sc)
 end subroutine mksupercell
end interface

interface
 subroutine polcart(red_ptot,pel,pel_cart,pelev,pion,pion_cart,polunit,&  
  &  ptot_cart,rprimd,ucvol,unit_out,usepaw)
  use defs_basis
  implicit none
  integer,intent(in) :: polunit
  integer,intent(in) :: unit_out
  integer,intent(in) :: usepaw
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: pel(3)
  real(dp),intent(out) :: pel_cart(3)
  real(dp),intent(in) :: pelev(3)
  real(dp),intent(in) :: pion(3)
  real(dp),intent(out) :: pion_cart(3)
  real(dp),intent(out) :: ptot_cart(3)
  real(dp),intent(in) :: red_ptot(3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine polcart
end interface

interface
 subroutine prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: ptgroupma
  integer,intent(in) :: spgroup
  integer,intent(in) :: bravais(11)
  real(dp),intent(inout) :: genafm(3)
 end subroutine prtspgroup
end interface

interface
 subroutine ptgmadata(ptgroupma,ptgrpmasb)
  implicit none
  integer,intent(in) :: ptgroupma
  character(len=10),intent(out) :: ptgrpmasb
 end subroutine ptgmadata
end interface

interface
 subroutine remove_inversion(nsym,symrel,tnons,nsym_out,symrel_out,tnons_out,pinv)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: nsym_out
  integer,intent(out) :: pinv
  integer,pointer :: symrel_out(:,:,:)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),pointer :: tnons_out(:,:)
 end subroutine remove_inversion
end interface

interface
 subroutine smallprim(metmin,minim,rprimd)
  use defs_basis
  implicit none
  real(dp),intent(out) :: metmin(3,3)
  real(dp),intent(out) :: minim(3,3)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine smallprim
end interface

interface
 subroutine spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,&  
  &  schsb,spgaxor,spgroup,sporder,spgorig)
  implicit none
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(out) :: sporder
  character(len=1),intent(out) :: brvsb
  character(len=15),intent(out) :: intsb
  character(len=35),intent(out) :: intsbl
  character(len=15),intent(out) :: ptintsb
  character(len=15),intent(out) :: ptschsb
  character(len=15),intent(out) :: schsb
 end subroutine spgdata
end interface

interface
 subroutine stresssym(gprimd,nsym,stress,sym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(inout) :: stress(6)
  integer,intent(in) :: sym(3,3,nsym)
 end subroutine stresssym
end interface

interface
 subroutine strconv(frac,gprimd,cart)
  use defs_basis
  implicit none
  real(dp),intent(inout) :: cart(6)
  real(dp),intent(in) :: frac(6)
  real(dp),intent(in) :: gprimd(3,3)
 end subroutine strconv
end interface

interface
 subroutine symanal(bravais,chkprim,genafm,msym,nsym,ptgroupma,rprimd,spgroup,symafm,symrel,tnons,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: chkprim
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(out) :: ptgroupma
  integer,intent(out) :: spgroup
  real(dp),intent(in) :: tolsym
  integer,intent(out) :: bravais(11)
  real(dp),intent(out) :: genafm(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symanal
end interface

interface
 subroutine symatm(indsym,natom,nsym,symrec,tnons,tolsym,typat,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  real(dp), intent(in) :: tolsym
  integer,intent(out) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrec(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine symatm
end interface

interface
 subroutine symaxes(center,iholohedry,isym,isymrelconv,label,ordersym,tnons_order,trialt,type_axis)
  use defs_basis
  implicit none
  integer,intent(in) :: center
  integer,intent(in) :: iholohedry
  integer,intent(in) :: isym
  integer,intent(in) :: ordersym
  integer,intent(in) :: tnons_order
  integer,intent(out) :: type_axis
  character(len=128),intent(out) :: label
  integer,intent(in) :: isymrelconv(3,3)
  real(dp),intent(in) :: trialt(3)
 end subroutine symaxes
end interface

interface
 subroutine symbrav(bravais,msym,nsym,ptgroup,rprimd,symrel,tolsym,axis)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  character(len=5),intent(out) :: ptgroup
  real(dp),intent(in) :: tolsym
  integer,optional,intent(out) :: axis(3)
  integer,intent(out) :: bravais(11)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,msym)
 end subroutine symbrav
end interface

interface
 subroutine symcharac(center, determinant, iholohedry, isym, label, symrel, tnons, type_axis)
  use defs_basis
  implicit none
  integer, intent(in) :: center
  integer, intent(in) :: determinant
  integer, intent(in) :: iholohedry
  integer, intent(in) :: isym
  integer, intent(out) :: type_axis
  character(len = 128), intent(out) :: label
  integer,intent(in) :: symrel(3,3)
  real(dp),intent(in) :: tnons(3)
 end subroutine symcharac
end interface

interface
 subroutine symchk(difmin,eatom,natom,tratom,transl,trtypat,typat,xred)
  use defs_basis
  implicit none
  integer,intent(out) :: eatom
  integer,intent(in) :: natom
  integer,intent(in) :: trtypat
  integer,intent(out) :: transl(3)
  real(dp),intent(out) :: difmin(3)
  real(dp),intent(in) :: tratom(3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine symchk
end interface

interface
 subroutine symdet(determinant,nsym,sym)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: determinant(nsym)
  integer,intent(in) :: sym(3,3,nsym)
 end subroutine symdet
end interface

interface
 subroutine symfind(berryopt,efield,gprimd,jellslab,msym,natom,noncoll,nptsym,nsym,&  
  &  nzchempot,ptsymrel,spinat,symafm,symrel,tnons,tolsym,typat,use_inversion,xred,&  
  &  nucdipmom)
  use defs_basis
  implicit none
  integer,intent(in) :: berryopt
  integer,intent(in) :: jellslab
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  integer,intent(in) :: noncoll
  integer,intent(in) :: nptsym
  integer,intent(out) :: nsym
  integer,intent(in) :: nzchempot
  integer,intent(in) :: use_inversion
  real(dp),intent(in) :: tolsym
  real(dp),intent(in) :: efield(3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),optional :: nucdipmom(3,natom)
  integer,intent(in) :: ptsymrel(3,3,msym)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine symfind
end interface

interface
 subroutine symkpt(chksymbreak,gmet,ibz2bz,iout,kbz,nkbz,nkibz,nsym,&  
  &  symrec,timrev,wtk,wtk_folded)
  use defs_basis
  implicit none
  integer,intent(in) :: chksymbreak
  integer,intent(in) :: iout
  integer,intent(in) :: nkbz
  integer,intent(out) :: nkibz
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(inout) :: ibz2bz(nkbz)
  real(dp),intent(in) :: kbz(3,nkbz)
  integer,intent(in) :: symrec(3,3,nsym)
  real(dp),intent(in) :: wtk(nkbz)
  real(dp),intent(out) :: wtk_folded(nkbz)
 end subroutine symkpt
end interface

interface
 subroutine symlatt(bravais,msym,nptsym,ptsymrel,rprimd,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(out) :: nptsym
  real(dp),intent(in) :: tolsym
  integer,intent(out) :: bravais(11)
  integer,intent(out) :: ptsymrel(3,3,msym)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine symlatt
end interface

interface
 subroutine symlist_bcc(additional_info,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: additional_info
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_bcc
end interface

interface
 subroutine symlist_fcc(nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_fcc
end interface

interface
 subroutine symlist_others(brvltt,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: brvltt
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_others
end interface

interface
 subroutine symlist_prim(additional_info,nsym,n_axes,spgroup)
  implicit none
  integer,intent(in) :: additional_info
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  integer,intent(in) :: n_axes(31)
 end subroutine symlist_prim
end interface

interface
 subroutine symmetrize_rprimd(bravais,nsym,rprimd,symrel,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: tolsym
  integer,intent(in) :: bravais(11)
  real(dp),intent(inout) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symmetrize_rprimd
end interface

interface
 subroutine symmetrize_xred(indsym,natom,nsym,symrel,tnons,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine symmetrize_xred
end interface

interface
 subroutine symplanes(center,iholohedry,isym,isymrelconv,itnonsconv,label,type_axis)
  use defs_basis
  implicit none
  integer,intent(in) :: center
  integer,intent(in) :: iholohedry
  integer,intent(in) :: isym
  integer,intent(out) :: type_axis
  character(len = 128), intent(out) :: label
  integer,intent(in) :: isymrelconv(3,3)
  real(dp),intent(in) :: itnonsconv(3)
 end subroutine symplanes
end interface

interface
 subroutine symptgroup(iholohedry,nsym,ptgroup,symrel)
  implicit none
  integer,intent(out) :: iholohedry
  integer,intent(in) :: nsym
  character(len=5),intent(out) :: ptgroup
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine symptgroup
end interface

interface
 subroutine symredcart(aprim,bprim,symcart,symred)
  use defs_basis
  implicit none
  integer,intent(in) :: symred(3,3)
  real(dp),intent(in) :: aprim(3,3)
  real(dp),intent(in) :: bprim(3,3)
  real(dp),intent(out) :: symcart(3,3)
 end subroutine symredcart
end interface

interface
 subroutine symrelrot(nsym,rprimd,rprimd_new,symrel,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  real(dp),intent(in) :: tolsym
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rprimd_new(3,3)
  integer,intent(inout) :: symrel(3,3,nsym)
 end subroutine symrelrot
end interface

interface
 subroutine symsgcube(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(inout) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(out) :: tnons(3,msym)
 end subroutine symsgcube
end interface

interface
 subroutine symsghexa(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(inout) :: brvltt
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symsghexa
end interface

interface
 subroutine symsgmono(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(inout) :: brvltt
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symsgmono
end interface

interface
 subroutine symsgortho(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symsgortho
end interface

interface
 subroutine symsgtetra(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&  
  &  spgroupma,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: msym
  integer,intent(in) :: nsym
  integer,intent(in) :: shubnikov
  integer,intent(in) :: spgaxor
  integer,intent(in) :: spgorig
  integer,intent(in) :: spgroup
  integer,intent(in) :: spgroupma
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
 end subroutine symsgtetra
end interface

interface
 subroutine symspgr(bravais,nsym,spgroup,symrel,tnons,tolsym)
  use defs_basis
  implicit none
  integer,intent(in) :: nsym
  integer,intent(out) :: spgroup
  real(dp),intent(in) :: tolsym
  integer,intent(in) :: bravais(11)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(inout) :: tnons(3,nsym)
 end subroutine symspgr
end interface

end module interfaces_41_geometry
!!***
