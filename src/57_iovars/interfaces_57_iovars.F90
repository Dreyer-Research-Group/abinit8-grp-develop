!!****m* ABINIT/interfaces_57_iovars
!! NAME
!! interfaces_57_iovars
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/57_iovars
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

module interfaces_57_iovars

 implicit none

interface
 subroutine chkinp(dtsets,iout,mpi_enregs,ndtset,ndtset_alloc,npsp,pspheads)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  type(mpi_type),intent(in) :: mpi_enregs(0:ndtset_alloc)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine chkinp
end interface

interface
 subroutine chkvars (string)
  implicit none
  character(len=*),intent(in) :: string
 end subroutine chkvars
end interface

interface
 subroutine chkvars_in_string(protocol, list_vars, list_logicals, list_strings, string)
  implicit none
  integer,intent(in) :: protocol
  character(len=*),intent(in) :: list_logicals
  character(len=*),intent(in) :: list_strings
  character(len=*),intent(in) :: list_vars
  character(len=*),intent(in) :: string
 end subroutine chkvars_in_string
end interface

interface
 subroutine compute_kgb_indicator(acc_kgb,bandpp,glb_comm,mband,mpw,npband,npfft,npslk,&  
  &  uselinalggpu)
  use defs_basis
  implicit none
  integer,intent(in) :: bandpp
  integer,intent(in) :: glb_comm
  integer,intent(in) :: mband
  integer,intent(in) :: mpw
  integer,intent(in) :: npband
  integer,intent(in) :: npfft
  integer,intent(inout) :: npslk
  integer,intent(inout) :: uselinalggpu
  real(dp),intent(inout) :: acc_kgb
 end subroutine compute_kgb_indicator
end interface

interface
 subroutine finddistrproc(dtsets,filnam,idtset,iexit,mband,mpi_enreg,ndtset_alloc,tread)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: idtset
  integer,intent(inout) :: iexit
  integer,intent(in) :: mband
  integer,intent(in) :: ndtset_alloc
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: tread(11)
  type(dataset_type),intent(inout),target :: dtsets(0:ndtset_alloc)
  character(len=fnlen),intent(in) :: filnam(5)
 end subroutine finddistrproc
end interface

interface
 subroutine getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,mixalch,nimage,npsp,npspalch,&  
  &  ntypat,ntypalch,pspheads)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: lmnmax
  integer,intent(out) :: lmnmaxso
  integer,intent(out) :: lnmax
  integer,intent(out) :: lnmaxso
  integer,intent(in) :: nimage
  integer,intent(in) :: npsp
  integer,intent(in) :: npspalch
  integer,intent(in) :: ntypalch
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: mixalch(npspalch,ntypalch,nimage)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine getdim_nloc
end interface

interface
 subroutine indefo(dtsets,ndtset_alloc,nprocs)
  use defs_abitypes
  implicit none
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: nprocs
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 end subroutine indefo
end interface

interface
 subroutine indefo1(dtset)
  use defs_abitypes
  implicit none
  type(dataset_type),intent(inout) :: dtset
 end subroutine indefo1
end interface

interface
 subroutine ingeo (acell,amu,dtset,bravais,&  
  &  genafm,iatfix,icoulomb,iimage,iout,jdtset,jellslab,lenstr,mixalch,&  
  &  msym,natom,nimage,npsp,npspalch,nspden,nsppol,nsym,ntypalch,ntypat,&  
  &  nucdipmom,nzchempot,pawspnorb,&  
  &  ptgroupma,ratsph,rprim,slabzbeg,slabzend,spgroup,spinat,string,symafm,&  
  &  symmorphi,symrel,tnons,tolsym,typat,vel,vel_cell,xred,znucl)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(out) :: icoulomb
  integer,intent(in) :: iimage
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(out) :: jellslab
  integer,intent(in) :: lenstr
  integer,intent(in) :: msym
  integer,intent(inout) :: natom
  integer,intent(in) :: nimage
  integer,intent(in) :: npsp
  integer,intent(in) :: npspalch
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(inout) :: nsym
  integer,intent(in) :: ntypalch
  integer,intent(in) :: ntypat
  integer,intent(in) :: nzchempot
  integer,intent(in) :: pawspnorb
  integer,intent(out) :: ptgroupma
  integer,intent(out) :: spgroup
  integer,intent(inout) :: symmorphi
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: slabzbeg
  real(dp),intent(out) :: slabzend
  character(len=*),intent(in) :: string
  real(dp),intent(out) :: tolsym
  integer,intent(out) :: bravais(11)
  real(dp),intent(out) :: acell(3)
  real(dp),intent(out) :: amu(ntypat)
  real(dp),intent(out) :: genafm(3)
  integer,intent(out) :: iatfix(3,natom)
  real(dp),intent(out) :: mixalch(npspalch,ntypalch)
  real(dp),intent(inout) :: nucdipmom(3,natom)
  real(dp),intent(in) :: ratsph(ntypat)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(inout) :: spinat(3,natom)
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
  integer,intent(out) :: typat(natom)
  real(dp),intent(out) :: vel(3,natom)
  real(dp),intent(out) :: vel_cell(3,3)
  real(dp),intent(out) :: xred(3,natom)
  real(dp),intent(in) :: znucl(npsp)
 end subroutine ingeo
end interface

interface
 subroutine ingeobld (iout,jdtset,lenstr,natrd,natom,nobj,string,typat,typat_read,xcart,xcart_read)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: lenstr
  integer,intent(in) :: natom
  integer,intent(in) :: natrd
  integer,intent(in) :: nobj
  character(len=*),intent(in) :: string
  integer,intent(out) :: typat(natom)
  integer,intent(in) :: typat_read(natrd)
  real(dp),intent(out) :: xcart(3,natom)
  real(dp),intent(in) :: xcart_read(3,natrd)
 end subroutine ingeobld
end interface

interface
 subroutine inkpts(bravais,chksymbreak,fockdownsampling,iout,iscf,istwfk,jdtset,&  
  &  kpt,kpthf,kptopt,kptnrm,kptrlatt_orig,kptrlatt,kptrlen,lenstr,msym,&  
  &  nkpt,nkpthf,nqpt,ngkpt,nshiftk,nshiftk_orig,shiftk_orig,nsym,&  
  &  occopt,qptn,response,rprimd,shiftk,string,symafm,symrel,vacuum,wtk,&  
  &  impose_istwf_1) ! Optional argument
  use defs_basis
  implicit none
  integer,intent(in) :: chksymbreak
  integer,intent(in),optional :: impose_istwf_1
  integer,intent(in) :: iout
  integer,intent(in) :: iscf
  integer,intent(in) :: jdtset
  integer,intent(in) :: kptopt
  integer,intent(in) :: lenstr
  integer,intent(in) :: msym
  integer,intent(inout) :: nkpt
  integer,intent(inout) :: nkpthf
  integer,intent(in) :: nqpt
  integer,intent(out) :: nshiftk
  integer,intent(out) :: nshiftk_orig
  integer,intent(in) :: nsym
  integer,intent(in) :: occopt
  integer,intent(in) :: response
  real(dp),intent(out) :: kptnrm
  real(dp),intent(out) :: kptrlen
  character(len=*),intent(in) :: string
  integer,intent(in) :: bravais(11)
  integer,intent(out) :: fockdownsampling(3)
  integer,intent(out) :: kptrlatt(3,3)
  integer,intent(out) :: kptrlatt_orig(3,3)
  integer,intent(out) :: ngkpt(3)
  integer,intent(in) :: vacuum(3)
  integer,intent(out) :: istwfk(nkpt)
  real(dp),intent(out) :: kpt(3,nkpt)
  real(dp),intent(out) :: kpthf(3,nkpthf)
  real(dp),intent(in) :: qptn(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: shiftk(3,210)
  real(dp),intent(out) :: shiftk_orig(3,210)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(out) :: wtk(nkpt)
 end subroutine inkpts
end interface

interface
 subroutine inqpt(chksymbreak,iout,jdtset,lenstr,msym,natom,qptn,wtqc,rprimd,spinat,string,typat,vacuum,xred,qptrlatt)
  use defs_basis
  implicit none
  integer,intent(in) :: chksymbreak
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: lenstr
  integer,intent(in) :: msym
  integer,intent(in) :: natom
  character(len=*),intent(in) :: string
  real(dp),intent(inout) :: wtqc
  integer,intent(inout) :: qptrlatt(3,3)
  integer,intent(in) :: vacuum(3)
  real(dp),intent(out) :: qptn(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: spinat(3,natom)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine inqpt
end interface

interface
 subroutine invacuum(jdtset,lenstr,natom,rprimd,string,vacuum,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: jdtset
  integer,intent(in) :: lenstr
  integer,intent(in) :: natom
  character(len=*),intent(in) :: string
  integer,intent(out) :: vacuum(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine invacuum
end interface

interface
 subroutine invars0(dtsets,istatr,istatshft,lenstr,&  
  &  msym,mxnatom,mxnimage,mxntypat,ndtset,ndtset_alloc,npsp,papiopt,timopt,string)
  use defs_abitypes
  implicit none
  integer,intent(out) :: istatr
  integer,intent(out) :: istatshft
  integer,intent(in) :: lenstr
  integer,intent(out) :: msym
  integer,intent(out) :: mxnatom
  integer,intent(out) :: mxnimage
  integer,intent(out) :: mxntypat
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(out) :: npsp
  integer,intent(out) :: papiopt
  integer,intent(inout) :: timopt
  character(len=*),intent(in) :: string
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 end subroutine invars0
end interface

interface
 subroutine invars1(bravais,dtset,iout,jdtset,lenstr,mband_upper,msym,npsp1,&  
  &  string,symafm,symrel,tnons,zionpsp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: lenstr
  integer,intent(out) :: mband_upper
  integer,intent(in) :: msym
  integer,intent(in) :: npsp1
  type(dataset_type),intent(inout) :: dtset
  character(len=*),intent(inout) :: string
  integer,intent(inout) :: bravais(11)
  integer,intent(inout) :: symafm(msym)
  integer,intent(inout) :: symrel(3,3,msym)
  real(dp),intent(inout) :: tnons(3,msym)
  real(dp),intent(in) :: zionpsp(npsp1)
 end subroutine invars1
end interface

interface
 subroutine invars1m(dmatpuflag,dtsets,iout,lenstr,mband_upper_,&  
  &  msym,mxga_n_rules,mxgw_nqlwl,mxlpawu,mxmband_upper,mxnatom,&  
  &  mxnatpawu,mxnatsph,mxnatsph_extra,mxnatvshift,mxnconeq,&  
  &  mxnimage,mxn_efmas_dirs,mxnkpt,mxnkptgw,mxnkpthf,mxnnos,mxnqptdm,mxnspinor,&  
  &  mxnsppol,mxnsym,mxntypat,mxnimfrqs,mxnfreqsp,mxnzchempot,&  
  &  mxn_projection_frequencies,ndtset,ndtset_alloc,string,npsp,zionpsp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(out) :: dmatpuflag
  integer,intent(in) :: iout
  integer,intent(in) :: lenstr
  integer,intent(in) :: msym
  integer,intent(out) :: mxga_n_rules
  integer,intent(out) :: mxgw_nqlwl
  integer,intent(out) :: mxlpawu
  integer,intent(out) :: mxmband_upper
  integer,intent(out) :: mxn_efmas_dirs
  integer,intent(out) :: mxn_projection_frequencies
  integer,intent(in) :: mxnatom
  integer,intent(out) :: mxnatpawu
  integer,intent(out) :: mxnatsph
  integer,intent(out) :: mxnatsph_extra
  integer,intent(out) :: mxnatvshift
  integer,intent(out) :: mxnconeq
  integer,intent(out) :: mxnfreqsp
  integer,intent(in) :: mxnimage
  integer,intent(out) :: mxnimfrqs
  integer,intent(out) :: mxnkpt
  integer,intent(out) :: mxnkptgw
  integer,intent(out) :: mxnkpthf
  integer,intent(out) :: mxnnos
  integer,intent(out) :: mxnqptdm
  integer,intent(out) :: mxnspinor
  integer,intent(out) :: mxnsppol
  integer,intent(out) :: mxnsym
  integer,intent(out) :: mxntypat
  integer,intent(out) :: mxnzchempot
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  character(len=*),intent(inout) :: string
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
  integer,intent(out) :: mband_upper_(0:ndtset_alloc)
  real(dp),intent(in) :: zionpsp(npsp)
 end subroutine invars1m
end interface

interface
 subroutine invars2(bravais,dtset,iout,jdtset,lenstr,&  
  &  mband,msym,npsp,string,usepaw,zionpsp)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: jdtset
  integer,intent(in) :: lenstr
  integer,intent(in) :: mband
  integer,intent(in) :: msym
  integer,intent(in) :: npsp
  integer,intent(in) :: usepaw
  type(dataset_type),intent(inout) :: dtset
  character(len=*),intent(in) :: string
  integer,intent(in) :: bravais(11)
  real(dp),intent(in) :: zionpsp(npsp)
 end subroutine invars2
end interface

interface
 subroutine invars2m(dtsets,iout,lenstr,mband_upper_,msym,ndtset,ndtset_alloc,npsp,pspheads,string)
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: lenstr
  integer,intent(in) :: msym
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  character(len=*),intent(in) :: string
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
  integer,intent(in) :: mband_upper_(0:ndtset_alloc)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine invars2m
end interface

interface
 subroutine macroin(dtsets,ecut_tmp,lenstr,ndtset_alloc,string)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: lenstr
  integer,intent(in) :: ndtset_alloc
  character(len=*),intent(inout) :: string
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
  real(dp),intent(in) :: ecut_tmp(3,2,10)
 end subroutine macroin
end interface

interface
 subroutine macroin2(dtsets,ndtset_alloc)
  use defs_abitypes
  implicit none
  integer,intent(in) :: ndtset_alloc
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
 end subroutine macroin2
end interface

interface
 subroutine mpi_setup(dtsets,filnam,lenstr,mpi_enregs,ndtset,ndtset_alloc,string)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: lenstr
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  character(len=*),intent(inout) :: string
  type(dataset_type),intent(inout) :: dtsets(0:ndtset_alloc)
  character(len=fnlen),intent(in) :: filnam(5)
  type(mpi_type),intent(inout) :: mpi_enregs(0:ndtset_alloc)
 end subroutine mpi_setup
end interface

interface
 subroutine out_acknowl(dtsets,iout,ndtset_alloc,npsp,pspheads) 
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  type(pspheader_type),intent(in) :: pspheads(npsp)
 end subroutine out_acknowl
end interface

interface
 subroutine outvar_a_h (choice,dmatpuflag,dtsets,iout,&  
  &  jdtset_,marr,multivals,mxvals,ncid,ndtset,ndtset_alloc,&  
  &  results_out,strimg)
  use m_results_out
  use defs_abitypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: dmatpuflag
  integer,intent(in) :: iout
  integer,intent(in) :: marr
  integer,intent(in) :: ncid
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  type(ab_dimensions),intent(in) :: multivals
  type(ab_dimensions),intent(in) :: mxvals
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  integer,intent(in) :: jdtset_(0:ndtset_alloc)
  type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
  character(len=8),intent(in) :: strimg(mxvals%nimage)
 end subroutine outvar_a_h
end interface

interface
 subroutine outvar_i_n (dtsets,iout,&  
  &  jdtset_,marr,multivals,mxvals,ncid,ndtset,ndtset_alloc,npsp,prtvol_glob,&  
  &  response_,results_out,strimg)
  use m_results_out
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: marr
  integer,intent(in) :: ncid
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  integer,intent(in) :: prtvol_glob
  type(ab_dimensions),intent(in) :: multivals
  type(ab_dimensions),intent(in) :: mxvals
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  integer,intent(in) :: jdtset_(0:ndtset_alloc)
  integer,intent(in) :: response_(ndtset_alloc)
  type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
  character(len=8),intent(in) :: strimg(mxvals%nimage)
 end subroutine outvar_i_n
end interface

interface
 subroutine outvar_o_z(choice,dtsets,iout,&  
  &  jdtset_,marr,multivals,mxvals,ncid,ndtset,ndtset_alloc,npsp,prtvol_glob,&  
  &  results_out,strimg,timopt)
  use m_results_out
  use defs_abitypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: iout
  integer,intent(in) :: marr
  integer,intent(in) :: ncid
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  integer,intent(in) :: prtvol_glob
  integer,intent(in) :: timopt
  type(ab_dimensions),intent(in) :: multivals
  type(ab_dimensions),intent(in) :: mxvals
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  integer,intent(in) :: jdtset_(0:ndtset_alloc)
  type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
  character(len=8),intent(in) :: strimg(mxvals%nimage)
 end subroutine outvar_o_z
end interface

interface
 subroutine outvars(choice,dmatpuflag,dtsets,filnam4,iout,&  
  &  mxvals,ndtset,ndtset_alloc,npsp,results_out,timopt)
  use m_results_out
  use defs_abitypes
  implicit none
  integer,intent(in) :: choice
  integer,intent(in) :: dmatpuflag
  integer,intent(in) :: iout
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  integer,intent(in) :: timopt
  character(len=*),intent(in) :: filnam4
  type(ab_dimensions),intent(in) :: mxvals
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
 end subroutine outvars
end interface

interface
 subroutine prtocc(dtsets,iout,jdtset_,ndtset_alloc,prtvol_glob,results_out)
  use m_results_out
  use defs_abitypes
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: prtvol_glob
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  integer,intent(in) :: jdtset_(0:ndtset_alloc)
  type(results_out_type),intent(in) :: results_out(0:ndtset_alloc)
 end subroutine prtocc
end interface

interface
 subroutine prttagm(dprarr,intarr,iout,jdtset_,length,&  
  &  marr,narr,narrm,ncid,ndtset_alloc,token,typevarphys,use_narrm,&  
  firstchar,forceprint)  ! optional
  use defs_basis
  implicit none
  integer,intent(in),optional :: forceprint
  integer,intent(in) :: iout
  integer,intent(in) :: length
  integer,intent(in) :: marr
  integer,intent(in) :: narr
  integer,intent(in) :: ncid
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: use_narrm
  character(len=1),intent(in),optional :: firstchar
  character(len=*),intent(in) :: token
  character(len=3),intent(in) :: typevarphys
  real(dp),intent(in) :: dprarr(marr,0:ndtset_alloc)
  integer,intent(in) :: intarr(marr,0:ndtset_alloc)
  integer,intent(in) :: jdtset_(0:ndtset_alloc)
  integer,intent(in) :: narrm(0:ndtset_alloc)
 end subroutine prttagm
end interface

interface
 subroutine prttagm_images(dprarr_images,iout,jdtset_,length,&  
  &  marr,narrm,ncid,ndtset_alloc,token,typevarphys,&  
  &  mxnimage,nimage,ndtset,prtimg,strimg,firstchar,forceprint)
  use defs_basis
  implicit none
  integer,intent(in),optional :: forceprint
  integer,intent(in) :: iout
  integer,intent(in) :: length
  integer,intent(in) :: marr
  integer,intent(in) :: mxnimage
  integer,intent(in) :: ncid
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  character(len=1),intent(in),optional :: firstchar
  character(len=*),intent(in) :: token
  character(len=3),intent(in) :: typevarphys
  real(dp),intent(in) :: dprarr_images(marr,mxnimage,0:ndtset_alloc)
  integer,intent(in) :: jdtset_(0:ndtset_alloc)
  integer,intent(in) :: narrm(0:ndtset_alloc)
  integer,intent(in) :: nimage(0:ndtset_alloc)
  integer,intent(in) :: prtimg(mxnimage,0:ndtset_alloc)
  character(len=8),intent(in) :: strimg(mxnimage)
 end subroutine prttagm_images
end interface

end module interfaces_57_iovars
!!***
