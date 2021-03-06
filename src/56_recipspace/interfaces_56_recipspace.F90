!!****m* ABINIT/interfaces_56_recipspace
!! NAME
!! interfaces_56_recipspace
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/56_recipspace
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

module interfaces_56_recipspace

 implicit none

interface
 subroutine cg_rotate(cryst,kpt1,isym,itimrev,shiftg,nspinor,ndat,&  
  npw1,kg1,npw2,kg2,istwf1,istwf2,cg1,cg2,work_ngfft,work)
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: istwf1
  integer,intent(in) :: istwf2
  integer,intent(in) :: isym
  integer,intent(in) :: itimrev
  integer,intent(in) :: ndat
  integer,intent(in) :: npw1
  integer,intent(in) :: npw2
  integer,intent(in) :: nspinor
  type(crystal_t),intent(in) :: cryst
  integer,intent(in) :: shiftg(3)
  integer,intent(in) :: work_ngfft(18)
  real(dp),intent(in) :: cg1(2,npw1,nspinor,ndat)
  real(dp),intent(out) :: cg2(2,npw2,nspinor,ndat)
  integer,intent(in) :: kg1(3,npw1)
  integer,intent(in) :: kg2(3,npw2)
  real(dp),intent(in) :: kpt1(3)
  real(dp),intent(out) :: work(2,work_ngfft(4),work_ngfft(5),work_ngfft(6))
 end subroutine cg_rotate
end interface

interface
 subroutine cinitylmgi3(gprimd,kg,kptns,mkmem,mpi_enreg,mpsang,mpw,&  
  &  nband,nkpt,npwarr,nsppol,optder,rprimd,ylm,ylm_gr)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: optder
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: rprimd(3,3)
  complex(dp),intent(out) :: ylm(mpw*mkmem,mpsang*mpsang)
  complex(dp),intent(out) :: ylm_gr(mpw*mkmem,3+6*(optder/2),mpsang*mpsang)
 end subroutine cinitylmgi3
end interface

interface
 subroutine get_full_kgrid(indkpt,kpt,kpt_fullbz,kptrlatt,nkpt,&  
  &  nkpt_fullbz,nshiftk,nsym,shiftk,symrel)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt_fullbz
  integer,intent(in) :: nshiftk
  integer,intent(in) :: nsym
  integer,intent(in) :: kptrlatt(3,3)
  integer,intent(out) :: indkpt(nkpt_fullbz)
  real(dp),intent(in) :: kpt(3,nkpt)
  real(dp),intent(out) :: kpt_fullbz(3,nkpt_fullbz)
  real(dp),intent(in) :: shiftk(3,nshiftk)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine get_full_kgrid
end interface

interface
 subroutine get_kpt_fullbz(kpt_fullbz,kptrlatt,nkpt_fullbz,nshiftk,shiftk)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt_fullbz
  integer,intent(in) :: nshiftk
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(out) :: kpt_fullbz(3,nkpt_fullbz)
  real(dp),intent(in) :: shiftk(3,nshiftk)
 end subroutine get_kpt_fullbz
end interface

interface
 subroutine getkgrid(chksymbreak,iout,iscf,kpt,kptopt,kptrlatt,kptrlen,&  
  &  msym,nkpt,nkpt_computed,nshiftk,nsym,rprimd,shiftk,symafm,symrel,vacuum,wtk,&  
  &  fullbz,nkpthf,kpthf,downsampling) ! optional
  use defs_basis
  implicit none
  integer,intent(in) :: chksymbreak
  integer,intent(in) :: iout
  integer,intent(in) :: iscf
  integer,intent(in) :: kptopt
  integer,intent(in) :: msym
  integer,intent(in) :: nkpt
  integer,intent(inout) :: nkpt_computed
  integer,intent(inout),optional :: nkpthf
  integer,intent(inout) :: nshiftk
  integer,intent(in) :: nsym
  real(dp),intent(out) :: kptrlen
  integer,optional,intent(in) :: downsampling(3)
  integer,intent(inout) :: kptrlatt(3,3)
  integer,intent(in) :: vacuum(3)
  real(dp),optional,allocatable,intent(out) :: fullbz(:,:)
  real(dp),intent(inout) :: kpt(3,nkpt)
  real(dp),optional,intent(out) :: kpthf(:,:)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: shiftk(3,210)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
  real(dp),intent(inout) :: wtk(nkpt)
 end subroutine getkgrid
end interface

interface
 subroutine initylmg(gprimd,kg,kptns,mkmem,mpi_enreg,mpsang,mpw,&  
  &  nband,nkpt,npwarr,nsppol,optder,rprimd,ylm,ylm_gr)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: optder
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),intent(in) :: kptns(3,nkpt)
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: ylm(mpw*mkmem,mpsang*mpsang)
  real(dp),intent(out) :: ylm_gr(mpw*mkmem,3+6*(optder/2),mpsang*mpsang)
 end subroutine initylmg
end interface

interface
 subroutine irrzg(irrzon,nspden,nsppol,nsym,n1,n2,n3,phnons,symafm,symrel,tnons)
  use defs_basis
  implicit none
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: n3
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(out) :: irrzon(n1*n2*n3,2,(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: phnons(2,n1*n2*n3,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine irrzg
end interface

interface
 subroutine laplacian(gprimd,mpi_enreg,nfft,nfunc,ngfft,paral_kgb,rdfuncr,&  
  &  laplacerdfuncr,rdfuncg_out,laplacerdfuncg_out,g2cart_out,rdfuncg_in,g2cart_in)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: nfunc
  integer,intent(in) :: paral_kgb
  type(mpi_type),intent(in) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in),optional,target :: g2cart_in(nfft)
  real(dp),intent(out),optional,target :: g2cart_out(nfft)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(out),optional,target :: laplacerdfuncg_out(2,nfft,nfunc)
  real(dp),intent(inout),optional :: laplacerdfuncr(nfft,nfunc)
  real(dp),intent(in),optional,target :: rdfuncg_in(2,nfft,nfunc)
  real(dp),intent(out),optional,target :: rdfuncg_out(2,nfft,nfunc)
  real(dp),intent(inout),optional,target :: rdfuncr(nfft,nfunc)
 end subroutine laplacian
end interface

interface
 subroutine listkk(dksqmax,gmet,indkk,kptns1,kptns2,nkpt1,nkpt2,nsym,&  
  &  sppoldbl,symafm,symmat,timrev,use_symrec)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt1
  integer,intent(in) :: nkpt2
  integer,intent(in) :: nsym
  integer,intent(in) :: sppoldbl
  integer,intent(in) :: timrev
  real(dp),intent(out) :: dksqmax
  logical,optional,intent(in) :: use_symrec
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(out) :: indkk(nkpt2*sppoldbl,6)
  real(dp),intent(in) :: kptns1(3,nkpt1)
  real(dp),intent(in) :: kptns2(3,nkpt2)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symmat(3,3,nsym)
 end subroutine listkk
end interface

interface
 subroutine mknucdipmom_k(gmet,kg,kpt,natom,nucdipmom,nucdipmom_k,npw,rprimd,ucvol,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: npw
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: kg(3,npw)
  real(dp),intent(in) :: kpt(3)
  real(dp),intent(in) :: nucdipmom(3,natom)
  complex(dpc),intent(out) :: nucdipmom_k(npw*(npw+1)/2)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine mknucdipmom_k
end interface

interface
 subroutine mkpwind_k(dk,dtset,fnkpt,fkptns,gmet,indkk_f2ibz,ikpt,ikpt1,&  
  &  kg,kgindex,mpi_enreg,npw_k,pwind_k1,symrec)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: fnkpt
  integer,intent(in) :: ikpt
  integer,intent(in) :: ikpt1
  integer,intent(in) :: npw_k
  type(dataset_type),intent(in) :: dtset
  type(mpi_type), intent(inout) :: mpi_enreg
  real(dp),intent(in) :: dk(3)
  real(dp),intent(in) :: fkptns(3,fnkpt)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: indkk_f2ibz(fnkpt,6)
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: kgindex(dtset%nkpt)
  integer,intent(out) :: pwind_k1(dtset%mpw)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
 end subroutine mkpwind_k
end interface

interface
 subroutine setmqgrid(mqgrid,mqgriddg,ecut,ecutdg,gprimd,nptsgvec,usepaw)
  use defs_basis
  implicit none
  integer , intent(inout) :: mqgrid
  integer , intent(inout) :: mqgriddg
  integer , intent(in) :: nptsgvec
  integer , intent(in) :: usepaw
  real(dp), intent(in) :: ecut
  real(dp), intent(in) :: ecutdg
  real(dp), intent(in) :: gprimd(3,3)
 end subroutine setmqgrid
end interface

interface
 subroutine setsym(indsym,irrzon,iscf,natom,nfft,ngfft,nspden,nsppol,nsym,phnons,&  
  &  symafm,symrec,symrel,tnons,typat,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: iscf
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ngfft(18)
  integer,intent(out) :: indsym(4,nsym,natom)
  integer,intent(inout) :: irrzon(nfft,2,(nspden/nsppol)-3*(nspden/4))
  real(dp),intent(out) :: phnons(2,nfft,(nspden/nsppol)-3*(nspden/4))
  integer,intent(in) :: symafm(nsym)
  integer,intent(out) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine setsym
end interface

interface
 subroutine smpbz(brav,iout,kptrlatt,mkpt,nkpt,nshiftk,option,shiftk,spkpt,downsampling)
  use defs_basis
  implicit none
  integer,intent(in) :: brav
  integer,intent(in) :: iout
  integer,intent(in) :: mkpt
  integer,intent(out) :: nkpt
  integer,intent(in) :: nshiftk
  integer,intent(in) :: option
  integer,optional,intent(in) :: downsampling(3)
  integer,intent(in) :: kptrlatt(3,3)
  real(dp),intent(in) :: shiftk(3,nshiftk)
  real(dp),intent(out) :: spkpt(3,mkpt)
 end subroutine smpbz
end interface

interface
 subroutine testkgrid(bravais,iout,kptrlatt,kptrlen,&  
  &  msym,nshiftk,nsym,prtkpt,rprimd,shiftk,symafm,symrel,vacuum)
  use defs_basis
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: msym
  integer,intent(out) :: nshiftk
  integer,intent(in) :: nsym
  integer,intent(in) :: prtkpt
  real(dp),intent(inout) :: kptrlen
  integer,intent(in) :: bravais(11)
  integer,intent(out) :: kptrlatt(3,3)
  integer,intent(in) :: vacuum(3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: shiftk(3,210)
  integer,intent(in) :: symafm(msym)
  integer,intent(in) :: symrel(3,3,msym)
 end subroutine testkgrid
end interface

end module interfaces_56_recipspace
!!***
