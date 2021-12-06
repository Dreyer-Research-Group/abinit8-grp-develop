!{\src2tex{textfont=tt}}
!!****f* ABINIT/initylmg
!! NAME
!! initylmg
!!
!! FUNCTION
!! Calculate the real spherical harmonics Ylm (and gradients)
!! over a set of (reciprocal space) (k+G) vectors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive
!!              translations (b^-1)
!!  kg(3,mpw)=integer coordinates of G vectors in basis sphere
!!  kptns(3,nkpt)=k points in terms of reciprocal translations
!!  mkmem =number of k points treated by this node
!!  mpi_enreg=information about MPI parallelization
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotential
!!  mpw   =maximum number of planewaves in basis sphere
!!         (large number)
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  nkpt  =number of k points
!!  npwarr(nkpt)=array holding npw for each k point
!!  nsppol=1 for unpolarized, 2 for polarized
!!  optder= 0=compute Ylm(K)
!!          1=compute Ylm(K) and dYlm/dKi
!!          2=compute Ylm(K), dYlm/dKi and d2Ylm/dKidKj
!!         -1=compute only dYlm/dKi
!!  rprimd(3,3)=dimensional primitive translations in real space
!!              (bohr)
!!  unkg=unit number for (k+G) sphere data
!!  unylm=unit number for storage of Ylm on disk
!!
!! OUTPUT
!!  if (optder>=0)
!!    ylm(mpw*mkmem,mpsang*mpsang) = real spherical harmonics
!!    for each G and k point
!!  if (optder>=1 or optder==-1)
!!    ylm_gr(mpw*mkmem,1:3,mpsang*mpsang)= gradients of real
!!    spherical harmonics wrt (G+k) in reduced coordinates
!!  if (optder>=2)
!!    ylm_gr(mpw*mkmem,4:9,mpsang*mpsang)= second gradients of
!!    real spherical harmonics wrt (G+k) in reduced coordinates
!!
!! NOTES
!! Remember the expression of complex spherical harmonics:
!! $Y_{lm}(%theta ,%phi)=sqrt{{(2l+1) over (4 %pi)}
!! {fact(l-m) over fact(l+m)} } P_l^m(cos(%theta))
!! func e^{i m %phi}$
!! Remember the expression of real spherical harmonics as
!!   linear combination of complex spherical harmonics:
!! $Yr_{lm}(%theta ,%phi)=(Re{Y_{l-m}}+(-1)^m Re{Y_{lm}})/sqrt{2}
!! $Yr_{l-m}(%theta ,%phi)=(Im{Y_{l-m}}-(-1)^m Im{Y_{lm}})/sqrt{2}
!!
!! CEDrev: This is altered version of initylmg to give COMPLEX Ylm 
!!        and CARTESIAN derivative up to 3rd
!!
!! PARENTS
!!      debug_tools,gstate,ks_ddiago,loop3dte,loper3,m_paw_pwij,m_shirley,m_wfs
!!      mover,nstpaw3,partial_dos_fractions,respfn,scfcv,wffile
!!
!! CHILDREN
!!      plm_coeff
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cinitylmgi3(gprimd,kg,kptns,mkmem,mpi_enreg,mpsang,mpw,&
&  nband,nkpt,npwarr,nsppol,optder,rprimd,ylm,ylm_gr)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_paw_sphharm, only : ass_leg_pol, plm_dtheta, plm_dphi, plm_coeff
 !use m_sphharm!, only : ass_leg_pol, plm_dtheta, plm_dphi, plm_coeff

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cinitylmgi3'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mkmem,mpsang,mpw,nkpt,nsppol,optder
! integer,intent(in) :: unkg,unylm
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),nband(nkpt*nsppol)
 integer,intent(in) :: npwarr(nkpt)
 real(dp),intent(in) :: gprimd(3,3),kptns(3,nkpt),rprimd(3,3)
 complex(dp),intent(out) :: ylm(mpw*mkmem,mpsang*mpsang)
 complex(dp),intent(out) :: ylm_gr(mpw*mkmem,3+6*(optder/2),mpsang*mpsang)
!Local variables ------------------------------
!scalars
 integer :: dimgr,ia,ib,ii,ikg,ikpt,ilang,ipw
 integer :: jj,kk,l0,ll
 integer :: me_distrb,mm,npw_k
 real(dp),parameter :: tol=1.d-10
 real(dp) :: cphi,ctheta,fact,onem,rr,sphi,stheta,work1,work2
 real(dp) :: xx,ylmcst,ylmcst2
 real(dp) :: yy,zz
 !character(len=500) :: message
!arrays

!CEDrev: 
! integer,parameter :: alpha(6)=(/1,2,3,3,3,2/)
! integer,parameter :: beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: alpha(6)=(/1,1,1,2,2,3/),beta(6)=(/1,2,3,2,3,3/)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: dphi(3),dtheta(3),iphase(mpsang-1),kpg(3)
 real(dp) :: rphase(mpsang-1)
 real(dp) :: blm(5,mpsang*mpsang)
 real(dp) :: ylmgr2_cart(3,3,2),ylmgr2_tmp(3,3)
 real(dp) :: ylmgr_cart(3,2)
 real(dp) :: ylmgr_red(10,2)

 !CEDrev: for thrid derivative
 integer,parameter :: alpha3(10)=(/1,1,1,1,1,1,2,2,2,3/)
 integer,parameter :: beta3(10)= (/1,1,1,2,2,3,2,2,3,3/)
 integer,parameter :: gamma3(10)=(/1,2,3,2,3,3,2,3,3,3/)
 integer :: pp,qq,ss,aa,igam
 real(dp) :: plm_d3t
 real(dp) :: plm_d2t(mpsang*mpsang)
 real(dp) :: dydth(2),d2ydth2(2),d3ydth3(2),dydph(2),d2ydph2(2)
 real(dp) :: d3ydph3(2),d2ydthdph(2),d3ydth2dph(2),d3ydthdph2(2)
 real(dp) :: ylmgr3_cart(3,3,3,2),ylmgr3_tmp(3,3,3)

 !TEST
 real(dp) :: ylmtstconst,ylmgr2_cart_tst(3,3,2)

!*****************************************************************

!Begin executable
 me_distrb=mpi_enreg%me_kpt
!Initialisation of spherical harmonics (and gradients)
 if (optder>=0) ylm(:,:)  =cmplx(zero,zero)
 if (optder/=0) ylm_gr(:,:,:)=cmplx(zero,zero)

!CEDrev: make room for third derivatives:
if (optder<3) then
   dimgr=3+6*(optder/2)
else if (optder==8) then
   dimgr=27
end if
! dimgr=3+6*(optder/2)

!Allocate some memory
! if (optder/=0) then
!   ABI_ALLOCATE(ylmgr_cart,(3,2))
! end if
! if (optder/=0.and.optder/=2) then
!   ABI_ALLOCATE(ylmgr_red,(3,2))
! end if
! if (optder==2) then
!   ABI_ALLOCATE(ylmgr2_cart,(3,3,2))
!   ABI_ALLOCATE(ylmgr2_tmp,(3,3))
!   ABI_ALLOCATE(ylmgr_red,(6,2))
!   ABI_ALLOCATE(blm,(5,mpsang*mpsang))
! end if

!CEDrev: For third derivative
! if (optder==8) then
!   ABI_ALLOCATE(ylmgr2_cart,(3,3,2))
!   ABI_ALLOCATE(ylmgr2_tmp,(3,3))
!   ABI_ALLOCATE(ylmgr_red,(10,2))
!   ABI_ALLOCATE(blm,(5,mpsang*mpsang))
! end if

!Loop over k-points:
 ikg=0
 do ikpt=1,nkpt

   ! Don't need this since its only for one kpoint
   !if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband(ikpt),-1,me_distrb)) cycle 

!  Get k+G-vectors, for this k-point:
   npw_k=npwarr(ikpt)
   ABI_ALLOCATE(kg_k,(3,npw_k))
   kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

!  Special case for l=0
   if (optder>=0) ylm(1+ikg:npw_k+ikg,1)=cmplx(1._dp/sqrt(four_pi),0)
   if (optder/=0) ylm_gr(1+ikg:npw_k+ikg,1:dimgr,1)=zero

   if (mpsang>1) then
!    Loop over all k+G
      do ipw=1,npw_k
        
!      Load k+G
       kpg(1)=kptns(1,ikpt)+real(kg_k(1,ipw),dp)
       kpg(2)=kptns(2,ikpt)+real(kg_k(2,ipw),dp)
       kpg(3)=kptns(3,ikpt)+real(kg_k(3,ipw),dp)

!      Calculate mod of k+G
       xx=two_pi*(gprimd(1,1)*kpg(1)+gprimd(1,2)*kpg(2)+gprimd(1,3)*kpg(3))
       yy=two_pi*(gprimd(2,1)*kpg(1)+gprimd(2,2)*kpg(2)+gprimd(2,3)*kpg(3))
       zz=two_pi*(gprimd(3,1)*kpg(1)+gprimd(3,2)*kpg(2)+gprimd(3,3)*kpg(3))
       rr=sqrt(xx**2+yy**2+zz**2)
       
       !TEST
       !if (ipw==57)then
       !   write(*,*) xx,yy,zz
       !end if

!      Continue only for k+G<>0
       if (rr>tol) then

!        Determine theta and phi
         cphi=one
         sphi=zero
         ctheta=zz/rr
         stheta=sqrt(abs((one-ctheta)*(one+ctheta)))
         if (stheta>tol) then
           cphi=xx/(rr*stheta)
           sphi=yy/(rr*stheta)
         end if
         do mm=1,mpsang-1
           rphase(mm)=dreal(dcmplx(cphi,sphi)**mm)
           iphase(mm)=aimag(dcmplx(cphi,sphi)**mm)
         end do

!        Determine gradients of theta and phi
         if (optder/=0) then
           dtheta(1)=ctheta*cphi
           dtheta(2)=ctheta*sphi
           dtheta(3)=-stheta
           dphi(1)=-sphi
           dphi(2)=cphi
           dphi(3)=zero
         end if

!        COMPUTE Ylm(K) 
!        ============================================
         if (optder>=0) then
!          Loop over angular momentum l
           do ilang=2,mpsang
             ll=ilang-1
             l0=ll**2+ll+1
             fact=1._dp/real(ll*(ll+1),dp)
             ylmcst=sqrt(real(2*ll+1,dp)/four_pi)
!            Special case m=0
             ylm(ikg+ipw,l0)=cmplx(ylmcst*ass_leg_pol(ll,0,ctheta),0)
!            Compute for m>0
             onem=one
             do mm=1,ll
               onem=-onem
               work1=ylmcst*sqrt(fact)*ass_leg_pol(ll,mm,ctheta)!*sqrt(2._dp)
               ylm(ikg+ipw,l0+mm)=cmplx(work1*rphase(mm),work1*iphase(mm))
               !TEST
               ylm(ikg+ipw,l0-mm)=-conjg(ylm(ikg+ipw,l0+mm))
!               ylm(ikg+ipw,l0-mm)=onem*conjg(cmplx(work1*rphase(mm),work1*iphase(mm)))
               if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
             end do ! End loop over m
           end do  ! End loop over l
         end if

!        COMPUTE dYlm/dKi
!        ============================================
         if (optder/=0) then
!          Loop over angular momentum l
           do ilang=2,mpsang
             ll=ilang-1
             l0=ll**2+ll+1
             fact=1._dp/real(ll*(ll+1),dp)
             ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/rr
!            === Special case m=0 ===
!            1-compute gradients in cartesian coordinates
             work1=ylmcst*plm_dtheta(ll,0,ctheta)
             ylmgr_cart(1:3,1)=work1*dtheta(1:3)
!            2-Transfer gradients into reduced coordinates
!             do ii=1,3
!               ylmgr_red(ii,1)=(rprimd(1,ii)*ylmgr_cart(1,1)+&
!&               rprimd(2,ii)*ylmgr_cart(2,1)+&
!&               rprimd(3,ii)*ylmgr_cart(3,1))
!             end do
!            3-Store gradients
             ylm_gr(ikg+ipw,1:3,l0) =cmplx(ylmgr_cart(1:3,1),0)
!             ylm_gr(ikg+ipw,1:3,l0) =cmplx(ylmgr_red(1:3,1),0)
!            === Compute for m>0 ===
             onem=one
             do mm=1,ll
               onem=-onem
!              1-compute gradients in cartesian coordinates
               work1=ylmcst*sqrt(fact)*plm_dtheta(ll,mm,ctheta)
               work2=ylmcst*sqrt(fact)*plm_dphi  (ll,mm,ctheta)
               ylmgr_cart(1:3,1)=rphase(mm)*work1*dtheta(1:3)-iphase(mm)*work2*dphi(1:3)
               ylmgr_cart(1:3,2)=iphase(mm)*work1*dtheta(1:3)+rphase(mm)*work2*dphi(1:3)
               
               !TEST
               !if (ll==1.and.mm==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   write(*,'(a15,4e16.4e2)') 'll=1,mm=1,ylm1 11', ylmgr_cart(1,1),ylmgr_cart(1,2)
               !else if (ll==2.and.mm==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   write(*,'(a15,4e16.4e2)') 'll=2,mm=1,ylm1 11', ylmgr_cart(1,1),ylmgr_cart(1,2)
               !else if (ll==2.and.mm==2.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   write(*,'(a15,4e16.4e2)') 'll=2,mm=2,ylm1 11', ylmgr_cart(1,1),ylmgr_cart(1,2)
               !end if

!              2-Transfer gradients into reduced coordinates
!               do kk=1,2
!                 do ii=1,3
!                   ylmgr_red(ii,kk)=(rprimd(1,ii)*ylmgr_cart(1,kk)+&
!&                   rprimd(2,ii)*ylmgr_cart(2,kk)+&
!&                   rprimd(3,ii)*ylmgr_cart(3,kk))
!                 end do
!               end do
!              3-Store gradients
               ylm_gr(ikg+ipw,1:3,l0+mm) =cmplx(ylmgr_cart(1:3,1),ylmgr_cart(1:3,2))
               ylm_gr(ikg+ipw,1:3,l0-mm) =-conjg(ylm_gr(ikg+ipw,1:3,l0+mm))
!               ylm_gr(ikg+ipw,1:3,l0-mm) =onem*conjg(cmplx(ylmgr_cart(1:3,1),ylmgr_cart(1:3,2)))
               if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
            end do ! End loop over m
          end do  ! End loop over l
        end if

!        COMPUTE d2Ylm/dKidKj
!        ============================================
         ! CEDrev:
        if (optder>=2) then 
!         if (optder==2) then
           call plm_coeff(blm,mpsang,ctheta)
!          Loop over angular momentum l
           do ilang=2,mpsang
             ll=ilang-1
             l0=ll**2+ll+1
             fact=1._dp/real(ll*(ll+1),dp)
             ylmcst=sqrt(real(2*ll+1,dp)/four_pi)/(rr**2)
!            === Special case m=0 ===
!            1-compute gradients in cartesian coordinates
             ylmgr2_cart(1,1,1)=ylmcst*(-blm(3,l0)*sphi*sphi+blm(4,l0)*cphi*cphi)
             ylmgr2_cart(2,2,1)=ylmcst*(-blm(3,l0)*cphi*cphi+blm(4,l0)*sphi*sphi)
             ylmgr2_cart(3,3,1)=ylmcst*blm(1,l0)
             ylmgr2_cart(3,1,1)=ylmcst*blm(2,l0)*cphi
             ylmgr2_cart(3,2,1)=ylmcst*blm(2,l0)*sphi
             ylmgr2_cart(2,1,1)=ylmcst*(blm(3,l0)+blm(4,l0))*sphi*cphi
             ylmgr2_cart(1,3,1)=ylmgr2_cart(3,1,1)
             ylmgr2_cart(1,2,1)=ylmgr2_cart(2,1,1)
             ylmgr2_cart(2,3,1)=ylmgr2_cart(3,2,1)

             !TEST: Check against my derivatives in cartesian coordinates
             !if (ll==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
             !   ylmtstconst=-sqrt(3./pi)/(2*rr**5)
             !   ylmgr2_cart_tst(1,1,1)=ylmtstconst*zz*(-2.*xx**2+yy**2+zz**2)
             !   ylmgr2_cart_tst(2,2,1)=ylmtstconst*zz*(xx**2-2.*yy**2+zz**2)
             !   ylmgr2_cart_tst(3,3,1)=ylmtstconst*3.*zz*(xx**2+yy**2)
             !   ylmgr2_cart_tst(3,1,1)=ylmtstconst*xx*(xx**2+yy**2-2.*zz**2)
             !   ylmgr2_cart_tst(3,2,1)=ylmtstconst*yy*(xx**2+yy**2-2.*zz**2)
             !   ylmgr2_cart_tst(2,1,1)=-ylmtstconst*(xx*yy*zz)
             !   ylmgr2_cart_tst(1,3,1)=ylmgr2_cart_tst(3,1,1)
             !   ylmgr2_cart_tst(1,2,1)=ylmgr2_cart_tst(2,1,1)
             !   ylmgr2_cart_tst(2,3,1)= ylmgr2_cart_tst(3,2,1)
                
             !   write(*,'(3i5)') kg_k(:,ipw)
                !write(*,'(a15,2e16.4e2)') 'ylm2 11',ylmgr2_cart(1,1,1),ylmgr2_cart_tst(1,1,1)
                !write(*,'(a15,2e16.4e2)') 'ylm2 22', ylmgr2_cart(2,2,1),ylmgr2_cart_tst(2,2,1)
                !write(*,'(a15,2e16.4e2)') 'ylm2 33', ylmgr2_cart(3,3,1),ylmgr2_cart_tst(3,3,1)
             !   write(*,'(a15,2e16.4e2)') 'ylm2 31', ylmgr2_cart(3,1,1),ylmgr2_cart_tst(3,1,1)
             !   write(*,*) 'l0',l0,'cphi',cphi,'blm(2,l0)',blm(2,l0)
             !   write(*,*) 'rr',rr,'xx',xx
                !write(*,'(a15,2e16.4e2)') 'ylm2 32', ylmgr2_cart(3,2,1),ylmgr2_cart_tst(3,2,1)
                !write(*,'(a15,2e16.4e2)') 'ylm2 21', ylmgr2_cart(2,1,1),ylmgr2_cart_tst(2,1,1)
                !stop
             !end if

!            2-Transfer gradients into reduced coordinates
!             do jj=1,3
!               do ii=1,3
!                 ylmgr2_tmp(ii,jj)=(rprimd(1,jj)*ylmgr2_cart(1,ii,1)+&
!&                 rprimd(2,jj)*ylmgr2_cart(2,ii,1)+&
!&                 rprimd(3,jj)*ylmgr2_cart(3,ii,1))
!               end do
          !end do
          do ii=1,6
               ia=alpha(ii);ib=beta(ii)
!               ylm_gr(ikg+ipw,4:9,l0) =cmplx(ylmgr2_cart(ia,ib,1),0)
               ylm_gr(ikg+ipw,3+ii,l0) =cmplx(ylmgr2_cart(ia,ib,1),0) 

!               ylmgr_red(ii,1)=(rprimd(1,ia)*ylmgr2_tmp(1,ib)+&
!&               rprimd(2,ia)*ylmgr2_tmp(2,ib)+&
!&               rprimd(3,ia)*ylmgr2_tmp(3,ib))
            end do
             !CEDrev: DO I NEED A SQRT 2 HERE????
!             ylm_gr(ikg+ipw,4:9,l0) =cmplx(ylmgr_cart(1:6,1),0)
!            === Compute for m>0 ===
             onem=one
             do mm=1,ll                
                onem=-onem
               ylmcst2=ylmcst*sqrt(fact)!*sqrt(two)

               ylmgr2_cart(1,1,1)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*rphase(mm)-&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
               ylmgr2_cart(1,1,2)=ylmcst2*((-blm(3,l0+mm)*sphi*sphi+blm(4,l0+mm)*cphi*cphi)*iphase(mm)+&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
               ylmgr2_cart(2,2,1)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*rphase(mm)+&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*iphase(mm))
               ylmgr2_cart(2,2,2)=ylmcst2*((-blm(3,l0+mm)*cphi*cphi+blm(4,l0+mm)*sphi*sphi)*iphase(mm)-&
&               blm(5,l0+mm)*2.d0*cphi*sphi*mm*rphase(mm))
               ylmgr2_cart(3,3,1)=ylmcst2*blm(1,l0+mm)*rphase(mm)
               ylmgr2_cart(3,3,2)=ylmcst2*blm(1,l0+mm)*iphase(mm)
               ylmgr2_cart(3,1,1)=ylmcst2*(blm(2,l0+mm)*cphi*rphase(mm)-&
&               mm*iphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(3,1,2)=ylmcst2*(blm(2,l0+mm)*cphi*iphase(mm)+&
&               mm*rphase(mm)*sphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(3,2,1)=ylmcst2*(blm(2,l0+mm)*sphi*rphase(mm)+&
&               mm*iphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
              ylmgr2_cart(3,2,2)=ylmcst2*(blm(2,l0+mm)*sphi*iphase(mm)-&
&               mm*rphase(mm)*cphi*onem*plm_dtheta(ll,mm,ctheta))
               ylmgr2_cart(2,1,1)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*rphase(mm)-&
&               blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*iphase(mm))
               ylmgr2_cart(2,1,2)=ylmcst2*((blm(3,l0+mm)+blm(4,l0+mm))*sphi*cphi*iphase(mm)+&
&               blm(5,l0+mm)*(sphi*sphi-cphi*cphi)*mm*rphase(mm))
               ylmgr2_cart(1,3,:)=ylmgr2_cart(3,1,:)
               ylmgr2_cart(1,2,:)=ylmgr2_cart(2,1,:)
               ylmgr2_cart(2,3,:)=ylmgr2_cart(3,2,:)
!              2-Transfer gradients into reduced coordinates
!               do kk=1,2
!                 do jj=1,3
!                   do ii=1,3
!                     ylmgr2_tmp(ii,jj)=(rprimd(1,jj)*ylmgr2_cart(1,ii,kk)+&
!&                     rprimd(2,jj)*ylmgr2_cart(2,ii,kk)+&
!&                     rprimd(3,jj)*ylmgr2_cart(3,ii,kk))
!                   end do
               !                 end do

               !TEST
               !if (ll==1.and.mm==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   ylmtstconst=sqrt(3./two_pi)/(2*rr**5)
               !   ylmgr2_cart_tst(3,1,1)=ylmtstconst*zz*(-2.*xx**2+yy**2+zz**2)
               !   ylmgr2_cart_tst(3,1,2)=ylmtstconst*zz*(-3.*xx*yy)
                  !write(*,'(a15,4e16.4e2)') 'ylm2 31', onem*ylmgr2_cart(3,1,1),ylmgr2_cart_tst(3,1,1),onem*ylmgr2_cart(3,1,2),ylmgr2_cart_tst(3,1,2)
               !   write(*,'(a15,4e16.4e2)') 'll=1,mm=1 ylm2 11', ylmgr2_cart(1,1,1),ylmgr2_cart(1,1,2)
                  !write(*,'(a15,4e16.4e2)') 'ylm2 12', ylmgr2_cart(2,1,1),ylmgr2_cart(2,1,2)
                  !write(*,*) 'l0+mm',l0+mm,'cphi',cphi,'blm(2,l0+mm)',blm(2,l0+mm)
                  !write(*,*) 'rr',rr,'xx',xx
                  !stop
               !end if
               !if (ll==2.and.mm==2.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
                  !write(*,'(a15,4e16.4e2)') 'ylm2 31', ylmgr2_cart(3,1,1),ylmgr2_cart(3,1,2)
               !   write(*,'(a15,4e16.4e2)') 'll=2,mm=2 ylm2 11', ylmgr2_cart(1,1,1),ylmgr2_cart(1,1,2)
                  !write(*,'(a15,4e16.4e2)') 'ylm2 12', ylmgr2_cart(2,1,1),ylmgr2_cart(2,1,2)
                  !stop
               !end if
               !if (ll==2.and.mm==1.and.kg_k(1,ipw)==1.and.kg_k(2,ipw)==1.and.kg_k(3,ipw)==1) then
               !   write(*,'(a15,4e16.4e2)') 'll=2,mm=1,ylm2 11', ylmgr2_cart(1,1,1),ylmgr2_cart(1,1,2)
               !end if


               do ii=1,6
                  ia=alpha(ii);ib=beta(ii)
                  
                  !TEST
                  ylm_gr(ikg+ipw,3+ii,l0+mm) =onem*cmplx(ylmgr2_cart(ia,ib,1),ylmgr2_cart(ia,ib,2))
                  ylm_gr(ikg+ipw,3+ii,l0-mm) =-conjg(ylm_gr(ikg+ipw,3+ii,l0+mm)) 

                  !ylm_gr(ikg+ipw,3+ii,l0+mm) =cmplx(ylmgr2_cart(ia,ib,1),ylmgr2_cart(ia,ib,2))
                  !ylm_gr(ikg+ipw,3+ii,l0-mm) =onem*conjg(cmplx(ylmgr2_cart(ia,ib,1),ylmgr2_cart(ia,ib,2)))

                  ! ylmgr_red(ii,kk)=(rprimd(1,ia)*ylmgr2_tmp(1,ib)+&
                  !& rprimd(2,ia)*ylmgr2_tmp(2,ib)+&
                  !& rprimd(3,ia)*ylmgr2_tmp(3,ib))
               end do
               !              end do
               !ylm_gr(ikg+ipw,4:9,l0+mm) =cmplx(ylmgr_cart(1:6,1),ylmgr_cart(1:6,2))
               !ylm_gr(ikg+ipw,4:9,l0-mm) =onem*conjg(cmplx(ylmgr_cart(1:6,1),ylmgr_cart(1:6,2)))
               if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
            end do ! End loop over m
         end do  ! End loop over l
      end if

!CEDrev: Since I'm not converting to reduced coordinates, I can clean this up significantly
!        COMPUTE d3Ylm/dKidKjKn
!        ============================================
        if (optder==8) then
           ! Loop over angular momentum l
           do ilang=2,mpsang
              ll=ilang-1
              l0=ll**2+ll+1
              ! === Special case m=0 ===
              
              ylmgr3_cart=zero

              if (ll==1) then
                 ylmcst=sqrt(3./(pi*4.))/(rr**7)
                 !xxx
                 ylmgr3_cart(1,1,1,1)=ylmcst*(-6.*zz*xx**3+9.*xx*zz*(yy**2+zz**2))
                 !xxy
                 ylmgr3_cart(1,1,2,1)=ylmcst*(3.*yy*zz*(-4.*xx**2+yy**2+zz**2))
                 ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                 ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                 !xxz
                 ylmgr3_cart(1,1,3,1)=ylmcst*(2.*xx**4-yy**4+(yy**2)*(zz**2)+2.*zz**4+(xx**2)*(yy**2-11.*zz**2))
                 ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                 ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                 !xyy
                 ylmgr3_cart(1,2,2,1)=ylmcst*(3.*xx*zz*(xx**2-4.*yy**2+zz**2))
                 ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                 ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                 !xyz
                 ylmgr3_cart(1,2,3,1)=ylmcst*(3.*xx*yy*(xx**2+yy**2-4.*zz**2))
                 ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                 !yyy
                 ylmgr3_cart(2,2,2,1)=ylmcst*(3.*yy*zz*(3.*xx**2-2.*yy**2+3.*zz**2))
                 !yyz
                 ylmgr3_cart(2,2,3,1)=ylmcst*(-xx**4+2.*yy**4-11.*(yy**2)*(zz**2)+2.*zz**4+(xx**2)*(yy**2+zz**2))
                 ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                 ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                 !yzz
                 ylmgr3_cart(2,3,3,1)=ylmcst*(3.*zz*yy*(3.*xx**2+3.*yy**2-2.*zz**2))
                 ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                 ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                 !zzx
                 ylmgr3_cart(3,3,1,1)=ylmcst*(3.*xx*zz*(3.*xx**2+3.*yy**2-2.*zz**2))
                 ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                 ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                 !zzz
                 ylmgr3_cart(3,3,3,1)=-ylmcst*3.*(xx**2+yy**2)*(xx**2+yy**2-4.*zz**2)

              else if (ll==2) then
                 ylmcst=sqrt(5./(pi*16.))/(rr**8)
                 !xxx
                 ylmgr3_cart(1,1,1,1)=ylmcst*72.*xx*(zz**2)*(-xx**2+yy**2+zz**2)
                 !xxy
                 ylmgr3_cart(1,1,2,1)=ylmcst*24.*(-5.*xx**2+yy**2+zz**2)
                 ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                 ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                 !xxz
                 ylmgr3_cart(1,1,3,1)=ylmcst*12.*zz*(3.*xx**4-yy**4+zz**4+2.*(xx**2)*(yy**2-4.*zz**2))
                 ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                 ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                 !xyy
                 ylmgr3_cart(1,2,2,1)=ylmcst*24.*xx*(zz**2)*(xx**2-5.*yy**2+zz**2)
                 ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                 ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                 !xyz
                 ylmgr3_cart(1,2,3,1)=ylmcst*48.*xx*yy*zz*(xx**2+yy**2-2.*zz**2)
                 ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                 !yyy
                 ylmgr3_cart(2,2,2,1)=ylmcst*72.*yy*(zz**2)*(xx**2-yy**2+zz**2)
                 !yyz
                 ylmgr3_cart(2,2,3,1)=ylmcst*12.*zz*(-xx**4+2.*(xx**2)*(yy**2)+3.*(yy**4)-8.*(yy**2)*(zz**2)+zz**4)
                 ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                 ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                 !yzz
                 ylmgr3_cart(2,3,3,1)=-ylmcst*12.*yy*(xx**4+yy**4-8.*(yy**2)*(zz**2)+3.*zz**4+2.*(xx**2)*(yy**2-4.*zz**2))
                 ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                 ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                 !zzx
                 ylmgr3_cart(3,3,1,1)=-ylmcst*12.*xx*(xx**4+yy**4-8.*(yy**2)*(zz**2)+3.*zz**4+2.*(xx**2)*(yy**2-4.*zz**2))
                 ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                 ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                 !zzz
                 ylmgr3_cart(3,3,3,1)=-ylmcst*72.*zz*(xx**2+yy**2)*(xx**2+yy**2-zz**2)

              else if (ll==3) then
                 !write(*,*) 'Third deriv of ylm not implemented for f electrons'
                 !stop

                 ylmcst=(3./4.)*sqrt(7./pi)/(rr**9)
                  
                 !xxx
                 ylmgr3_cart(1,1,1,1)=ylmcst*xx*zz*(6.*xx**4 - 9.*yy**4 + 57.*yy**2*zz**2 + 66.*zz**4 - xx**2*(3.*yy**2 + 103.*zz**2))
                 !xxy
                 ylmgr3_cart(1,1,2,1)=ylmcst*yy*zz*(12.*xx**4 - 3.*yy**4 + 19.*yy**2*zz**2 + 22.*zz**4 + 3.*xx**2*(3.*yy**2 - 47.*zz**2))
                 ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                 ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                 !xxz
                 ylmgr3_cart(1,1,3,1)=ylmcst*(-2.)*xx**6 + yy**6 - 15.*yy**4*zz**2 - 8.*yy**2*zz**4 + 8.*zz**6 - 3.*xx**4*(yy**2 - 23.*zz**2) + 6.*xx**2*(9.*yy**2*zz**2 - 16.*zz**4)
                 ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                 ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                 !xyy
                 ylmgr3_cart(1,2,2,1)=ylmcst*xx*zz*((-3.)*xx**4 + 12.*yy**4 - 141.*yy**2*zz**2 + 22.*zz**4 + xx**2*(9.*yy**2 + 19.*zz**2))
                 ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                 ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                 !xyz
                 ylmgr3_cart(1,2,3,1)=ylmcst*(-1)*(xx*yy*(3.*xx**4 + 3.*yy**4 - 84.*yy**2*zz**2 + 88.*zz**4 + 6.*xx**2*(yy**2 - 14.*zz**2)))
                 ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                 ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                 !yyy
                 ylmgr3_cart(2,2,2,1)=ylmcst*yy*zz*(-9.*xx**4 + 6.*yy**4 - 103.*yy**2*zz**2 + 66.*zz**4 - 3.*xx**2*(yy**2 - 19.*zz**2))
                 !yyz
                 ylmgr3_cart(2,2,3,1)=ylmcst*xx**6 - 2.*yy**6 - 15.*xx**4*zz**2 + 69.*yy**4*zz**2 - 96.*yy**2*zz**4 + 8.*zz**6 + xx**2*(-3.*yy**4 + 54.*yy**2*zz**2 - 8.*zz**4)
                 ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                 ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                 !yzz
                 ylmgr3_cart(2,3,3,1)=-ylmcst*yy*zz*(39.*xx**4 + 39.*yy**4 - 112.*yy**2*zz**2 + 24.*zz**4 + 2.*xx**2*(39.*yy**2 - 56.*zz**2))
                 ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                 ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                 !zzx
                 ylmgr3_cart(3,3,1,1)=-ylmcst*(xx*zz*(39.*xx**4 + 39.*yy**4 - 112.*yy**2*zz**2 + 24.*zz**4 + 2.*xx**2*(39.*yy**2 - 56.*zz**2)))
                 ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                 ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                 !zzz
                 ylmgr3_cart(3,3,3,1)= ylmcst*(xx**2 + yy**2)*(13.*xx**4 + 13.*yy**4 - 114.*yy**2*zz**2 + 48.*zz**4 + 2.*xx**2*(13.*yy**2 - 57.*zz**2))



              end if

              ! 2-Transfer gradients into reduced coordinates
              !ylmgr3_tmp=zero
              !do ii=1,3
              !   do jj=1,3
              !      do kk=1,3
              !         do pp=1,3
              !            do qq=1,3
              !               do ss=1,3
              !                  ylmgr3_tmp(ii,jj,kk)=ylmgr3_tmp(ii,jj,kk) &
              !                       & +rprimd(pp,ii)*rprimd(qq,jj)*rprimd(ss,kk)*ylmgr3_cart(pp,qq,ss,1)
              !               end do
              !            end do
              !         end do
              !      end do
              !   end do
              !end do
              do ii=1,10
                 ia=alpha3(ii);ib=beta3(ii);igam=gamma3(ii)
                 ylm_gr(ikg+ipw,9+ii,l0)=cmplx(ylmgr3_cart(ia,ib,igam,1),0)
                 !ylmgr_cart(ii,1)=ylmgr3_tmp(ia,ib,igam) 
                 !ylmgr_red(ii,1)=ylmgr3_tmp(ia,ib,igam)
              end do
              !ylm_gr(ikg+ipw,10:19,l0)=cmplx(ylmgr_red(1:10,1),0)
              ! === Compute for m>0 ===
              onem=one
              do mm=1,ll
                 onem=-onem

                 ! === elements ===
                 ylmgr3_cart=zero
                 
                 if (ll==1 .and. mm==1) then

                    ylmcst=sqrt(3./(pi*8.))/(rr**7)

                    !xxx
                    ylmgr3_cart(1,1,1,1)=ylmcst*3.*(-4.*(xx**2)*(yy**2+zz**2)+(yy**2+zz**2)**2)
                    ylmgr3_cart(1,1,1,2)=ylmcst*3.*(2.*(xx**3)*yy-3.*xx*yy*(yy**2+zz**2))
                    !xxy
                    ylmgr3_cart(1,1,2,1)=-ylmcst*(-6.*yy*xx**3+9.*xx*yy*(yy**2+zz**2))
                    ylmgr3_cart(1,1,2,2)=-ylmcst*(2.*xx**4+2.*yy**4+(yy**2)*(zz**2)-zz**4+(xx**2)*(-11.*yy**2+zz**2))
                    ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                    ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                    !xxz
                    ylmgr3_cart(1,1,3,1)=-ylmcst*3.*zz*(-2.*xx**3+3.*xx*(yy**2+zz**2))
                    ylmgr3_cart(1,1,3,2)=-ylmcst*3.*zz*(-4.*yy*xx**2+yy*(yy**2+zz**2))
                    ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                    ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                    !xyy
                    ylmgr3_cart(1,2,2,1)=ylmcst*(-2.*xx**4-2.*yy**4-(yy**2)*(zz**2)+zz**4+(xx**2)*(11.*yy**2-zz**2))
                    ylmgr3_cart(1,2,2,2)=ylmcst*(-9.*yy*xx**3+3.*xx*(2.*yy**3-3.*yy*zz**2))
                    ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                    ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                    !xyz
                    ylmgr3_cart(1,2,3,1)=-ylmcst*(-12.*zz*yy*xx**2+3.*zz*yy*(yy**2+zz**2))
                    ylmgr3_cart(1,2,3,2)=-ylmcst*3.*zz*(xx**3+xx*(-4.*yy**2+zz**2))
                    ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                    !yyy
                    ylmgr3_cart(2,2,2,1)=ylmcst*3.*(-3.*yy*xx**3+xx*(2.*yy**3-3.*yy*zz**2))
                    ylmgr3_cart(2,2,2,2)=ylmcst*3.*(xx**4-4.*(yy**2)*(zz**2)+zz**4+(xx**2)*(-4*yy**2+2.*zz**2))
                    !yyz
                    ylmgr3_cart(2,2,3,1)=-ylmcst*3.*zz*(xx**3+xx*(-4.*yy**2+zz**2))
                    ylmgr3_cart(2,2,3,2)=-ylmcst*3.*zz*(3.*yy*xx**2-2.*yy**3+3.*yy*zz**2)
                    ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                    ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                    !yzz
                    ylmgr3_cart(2,3,3,1)=ylmcst*(-3.*yy*xx**3-3.*xx*(yy**3-4.*yy*zz**2))
                    ylmgr3_cart(2,3,3,2)=ylmcst*(xx**4-2.*yy**4+11.*(yy**2)*(zz**2)-2.*zz**4-(xx**2)*(yy**2+zz**2))
                    ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                    ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                    !zzx
                    ylmgr3_cart(3,3,1,1)=ylmcst*(-2.*xx**4+yy**4-(yy**2)*(zz**2)-2.*zz**4-(xx**2)*(yy**2-11.*zz**2))
                    ylmgr3_cart(3,3,1,2)=ylmcst*(-3.*yy*xx**3-3.*xx*(yy**3-4.*yy*zz**2))
                    ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                    ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                    !zzz
                    ylmgr3_cart(3,3,3,1)=-ylmcst*3.*xx*zz*(3.*xx**2+3.*yy**2-2.*zz**2)
                    ylmgr3_cart(3,3,3,2)=-ylmcst*3.*yy*zz*(3.*xx**2+3.*yy**2-2.*zz**2)

                 else if (ll==2 .and. mm==1) then

                    ylmcst=sqrt(15./(pi*8.))/(rr**8)

                    !xxx
                    ylmgr3_cart(1,1,1,1)=ylmcst*6.*zz*(xx**4-6.*(xx**2)*(yy**2+zz**2)+(yy**2+zz**2)**2)
                    ylmgr3_cart(1,1,1,2)=ylmcst*6.*zz*(4.*yy*xx**3-4.*xx*yy*(yy**2+zz**2))
                    !xxy
                    ylmgr3_cart(1,1,2,1)=-ylmcst*2.*zz*(-12.*yy*xx**3+12.*xx*yy*(yy**2+zz**2))
                    ylmgr3_cart(1,1,2,2)=-ylmcst*2.*zz*(3.*xx**4+3.*yy**4+2.*(yy**2)*(zz**2)-zz**4+2.*(xx**2)*(-9*yy**2+zz**2))
                    ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                    ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                    !xxz
                    ylmgr3_cart(1,1,3,1)=-ylmcst*2.*(xx**5-2.*(xx**3)*(yy**2+7.*zz**2)+xx*(-3.*yy**4+6.*(yy**2)*(zz**2)+9.*zz**4))
                    ylmgr3_cart(1,1,3,2)=-ylmcst*2.*(3.*yy*xx**4+2.*(xx**2)*(yy**3-9.*yy*zz**2)-(yy**5-2.*(yy**3)*(zz**2)-3*yy*zz**4))
                    ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                    ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                    !xyy
                    ylmgr3_cart(1,2,2,1)=ylmcst*2.*zz*(-3.*xx**4-3.*yy**4-2.*(yy**2)*(zz**2)+zz**4+2.*(xx**2)*(9.*yy**2-zz**2))
                    ylmgr3_cart(1,2,2,2)=ylmcst*2.*zz*(-12.*yy*xx**3+12*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                    ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                    !xyz
                    ylmgr3_cart(1,2,3,1)=ylmcst*2.*(-3.*yy*xx**4-2.*(xx**2)*(yy**3-9*yy*zz**2)+yy**5-2.*(yy**3)*(zz**2)-3.*yy*zz**4)
                    ylmgr3_cart(1,2,3,2)=ylmcst*2.*(xx**5-2.*(xx**3)*(yy**2+zz**2)-3*xx*(yy**4-6.*(yy**2)*(zz**2)+zz**4))
                    ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                    !yyy
                    ylmgr3_cart(2,2,2,1)=ylmcst*6.*zz*(-4.*yy*xx**3+4.*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(2,2,2,2)=ylmcst*6.*zz*(xx**4+yy**4-6.*(yy**2)*(zz**2)+zz**4+(xx**2)*(-6.*yy**2+2.*zz**2))
                    !yyz
                    ylmgr3_cart(2,2,3,1)=ylmcst*2.*(xx**5-2.*(xx**3)*(yy**2+zz**2)-3.*xx*(yy**4-6.*(yy**2)*(zz**2)+zz**4))
                    ylmgr3_cart(2,2,3,2)=ylmcst*2.*(3.*yy*xx**4+2.*(xx**2)*(yy**3-3.*yy*zz**2)-(yy**5-14*(yy**3)*(zz**2)+9.*yy*zz**4))
                    ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                    ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                    !yzz  
                    ylmgr3_cart(2,3,3,1)=ylmcst*2.*zz*(-12.*yy*xx**3-12*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(2,3,3,2)=ylmcst*2.*zz*(3.*xx**4-9.*yy**4+14.*(yy**2)*(zz**2)-zz**4+(xx**2)*(-6.*yy**2+2.*zz**2))
                    ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                    ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                    !zzx
                    ylmgr3_cart(3,3,1,1)=-ylmcst*2.*zz*(9.*xx**4-3.*yy**4-2.*(yy**2)*(zz**2)+zz**4+2.*(xx**2)*(3.*yy**2-7.*zz**2))
                    ylmgr3_cart(3,3,1,2)=-ylmcst*2.*zz*(12.*yy*xx**3+12*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                    ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                    !zzz
                    ylmgr3_cart(3,3,3,1)=ylmcst*6.*xx*(xx**4+yy**4-6.*(yy**2)*(zz**2)+zz**4+2.*(xx**2)*(yy**2-3.*zz**2))
                    ylmgr3_cart(3,3,3,2)=ylmcst*6.*yy*(xx**4+yy**4-6.*(yy**2)*(zz**2)+zz**4+2.*(xx**2)*(yy**2-3.*zz**2))

                 else if (ll==2 .and. mm==2) then

                    ylmcst=sqrt(15./(pi*32.))/(rr**8)

                    !xxx
                    ylmgr3_cart(1,1,1,1)=-ylmcst*12.*(-2.*(xx**3)*(2.*yy**2+zz**2)+2.*xx*(2.*yy**4+3.*(yy**2)*(zz**2)+zz**4))
                    ylmgr3_cart(1,1,1,2)=-ylmcst*12.*(yy*xx**4-6.*yy*(xx**2)*(yy**2+zz**2)+yy*(yy**2+zz**2)**2)
                    !xxy
                    ylmgr3_cart(1,1,2,1)=ylmcst*4.*(-6.*yy*xx**4-2.*(yy**3)*(yy**2+zz**2)+2.*(xx**2)*(8.*yy**3+3.*yy*zz**2))
                    ylmgr3_cart(1,1,2,2)=ylmcst*4.*(xx**5-2.*(xx**3)*(7.*yy**2+zz**2)+xx*(9.*yy**4+6.*(yy**2)*(zz**2)-3.*zz**4))
                    ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                    ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                    !xxz
                    ylmgr3_cart(1,1,3,1)=-ylmcst*4.*zz*(3.*xx**4+3.*yy**4+4.*(yy**2)*(zz**2)+zz**4-2.*(xx**2)*(9.*yy**2+4.*zz**2))
                    ylmgr3_cart(1,1,3,2)=-ylmcst*4.*zz*(12.*yy*xx**3-12.*xx*yy*(yy**2+zz**2))
                    ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                    ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                    !xyy
                    ylmgr3_cart(1,2,2,1)=ylmcst*4.*(2.*xx**5+2.*(xx**3)*(-8.*yy**2+zz**2)+6.*xx*(yy**4-(yy**2)*(zz**2)))
                    ylmgr3_cart(1,2,2,2)=ylmcst*4.*(9.*yy*xx**4-2.*(xx**2)*(7.*yy**3-3.*yy*zz**2)+yy**5-2.*(yy**3)*(zz**2)-3.*yy*zz**4)
                    ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                    ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                    !xyz
                    ylmgr3_cart(1,2,3,1)=ylmcst*4.*zz*(-12.*yy*xx**3+12.*xx*yy**3)
                    ylmgr3_cart(1,2,3,2)=ylmcst*4.*zz*(3.*xx**4+3.*yy**4+2.*(yy**2)*(zz**2)-zz**4+2.*(xx**2)*(-9.*yy**2+zz**2))
                    ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                    !yyy
                    ylmgr3_cart(2,2,2,1)=-ylmcst*12.*(-4.*yy*xx**4+2.*yy*(zz**2)*(yy**2-zz**2)+2.*(xx**2)*(2.*yy**3-3.*yy*zz**2))
                    ylmgr3_cart(2,2,2,2)=-ylmcst*12.*(xx**5+(xx**3)*(-6.*yy**2+2.*zz**2)+xx*(yy**4-6.*(yy**2)*(zz**2)+zz**4))
                    !yyz
                    ylmgr3_cart(2,2,3,1)=ylmcst*4.*zz*(3.*xx**4+3.*yy**4-8.*(yy**2)*(zz**2)+zz**4+(xx**2)*(-18*yy**2+4.*zz**2))
                    ylmgr3_cart(2,2,3,2)=ylmcst*4.*zz*(12*yy*xx**3-12.*xx*(yy**3-yy*zz**2))
                    ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                    ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                    !yzz
                    ylmgr3_cart(2,3,3,1)=ylmcst*4.*(yy*(xx**4-yy**4-2.*(xx**2)*(zz**2)+8.*(yy**2)*(zz**2)-3.*zz**4)+xx*(2.*yy*xx**3+2.*xx*(yy**3-5.*yy*zz**2)))
                    ylmgr3_cart(2,3,3,2)=ylmcst*4.*(-xx*(xx**4-yy**4-2.*(xx**2)*(zz**2)+8.*(yy**2)*(zz**2)-3.*zz**4)+yy*(2.*yy*xx**3+2.*xx*(yy**3-5.*yy*zz**2)))
                    ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                    ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                    !zzx
                    ylmgr3_cart(3,3,1,1)=ylmcst*4.*(xx*(xx**4-yy**4-8.*(xx**2)*(zz**2)+2.*(yy**2)*(zz**2)+3.*zz**4)+yy*(-2.*yy*xx**3-2.*xx*(yy**3-5.*yy*zz**2)))
                    ylmgr3_cart(3,3,1,2)=ylmcst*4.*(yy*(xx**4-yy**4-8.*(xx**2)*(zz**2)+2.*(yy**2)*(zz**2)+3.*zz**4)+xx*(2.*yy*xx**3+2.*xx*(yy**3-5.*yy*zz**2)))
                    ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                    ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                    !zzz
                    ylmgr3_cart(3,3,3,1)=ylmcst*24.*zz*xx*(xx**2+yy**2-zz**2)
                    ylmgr3_cart(3,3,3,2)=ylmcst*24.*zz*yy*(xx**2+yy**2-zz**2)

                 else if (ll==3 .and. mm==1) then
                    ylmcst=sqrt(3./2.)*(3./4.)*sqrt(7./pi)/(rr**9)

                    !xxx
                    ylmgr3_cart(1,1,1,1)=ylmcst*(yy**2 - 14.*zz**2)*(yy**2 + zz**2)**2 - 4.*xx**4*(yy**2 + 11.*zz**2) - 3.*xx**2*(yy**4 - 38.*yy**2*zz**2 - 39.*zz**4)
                    !xxy
                    ylmgr3_cart(1,1,2,1)=ylmcst*xx*yy*(2.*xx**4 - 3.*yy**4 + 69.*yy**2*zz**2 + 72.*zz**4 - xx**2*(yy**2 + 101.*zz**2))
                    ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                    ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                    !xxz
                    ylmgr3_cart(1,1,3,1)=ylmcst*xx*zz*(22.*xx**4 - 33.*yy**4 + 9.*yy**2*zz**2 + 42.*zz**4 - xx**2*(11.*yy**2 + 111.*zz**2))
                    ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                    ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                    !xyy
                    ylmgr3_cart(1,2,2,1)=ylmcst*((-2.)*xx**6 - 2.*yy**6 + 57.*yy**4*zz**2 + 45.*yy**2*zz**4 - 14.*zz**6 + xx**4*(9.*yy**2 + 57.*zz**2) + 9.*xx**2*(yy**4 - 44.*yy**2*zz**2 + 5.*zz**4))/3.
                    ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                    ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                    !xyz
                    ylmgr3_cart(1,2,3,1)=ylmcst*yy*zz*(44.*xx**4 - 11.*yy**4 + 3.*yy**2*zz**2 + 14.*zz**4 + 3.*xx**2*(11.*yy**2 - 39.*zz**2))
                    ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                    !yyy
                    ylmgr3_cart(2,2,2,1)=-ylmcst*(xx*yy*(3.*xx**4 - 2.*yy**4 + 101.*yy**2*zz**2 - 72.*zz**4 + xx**2*(yy**2 - 69.*zz**2)))
                    !yyz
                    ylmgr3_cart(2,2,3,1)=ylmcst*xx*zz*((-11.)*xx**4 + 44.*yy**4 - 117.*yy**2*zz**2 + 14.*zz**4 + 3.*xx**2*(11.*yy**2 + zz**2))
                    ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                    ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                    !yzz
                    ylmgr3_cart(2,3,3,1)=-ylmcst*(xx*yy*(11.*xx**4 + 11.*yy**4 - 108.*yy**2*zz**2 + 56.*zz**4 + 2.*xx**2*(11.*yy**2 - 54.*zz**2)))
                    ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                    ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                    !zzx
                    ylmgr3_cart(3,3,1,1)=ylmcst*(-22.*xx**6 + 11.*yy**6 - 45.*yy**4*zz**2 - 48.*yy**2*zz**4 + 8.*zz**6 + 3.*xx**4*(-11.*yy**2 + 93.*zz**2) + 3.*xx**2*(78.*yy**2*zz**2 - 72.*zz**4))/3.
                    ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                    ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                    !zzz
                    ylmgr3_cart(3,3,3,1)= -ylmcst*(xx*zz*(63.*xx**4 + 63.*yy**4 - 104.*yy**2*zz**2 + 8.*zz**4 + 2.*xx**2*(63.*yy**2 - 52.*zz**2)))

                 else if (ll==3 .and. mm==2) then
                    ylmcst=sqrt(15.)*(3./4.)*sqrt(7./pi)/(rr**9)

                    !xxx
                    ylmgr3_cart(1,1,1,1)=ylmcst*xx*zz*(-2.*xx**4 + xx**2*(41.*yy**2 + 21.*zz**2) - 3.*(9.*yy**4 + 13.*yy**2*zz**2 + 4.*zz**4))
                    !xxy
                    ylmgr3_cart(1,1,2,1)=-ylmcst*5.*yy*zz*(4.*xx**4 + yy**2*(yy**2 + zz**2) - 3.*xx**2*(3.*yy**2 + zz**2))
                    ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                    ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                    !xxz
                    ylmgr3_cart(1,1,3,1)=ylmcst*(2.*xx**6)/3. + (5.*yy**6)/3. - 3.*yy**4*zz**2 - 6.*yy**2*zz**4 - (4.*zz**6)/3. - &
                         & xx**4*(7.*yy**2 + 15.*zz**2) - 6.*xx**2*(yy**4 - 7.*yy**2*zz**2 - 3.*zz**4)
                    ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                    ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                    !xyy
                    ylmgr3_cart(1,2,2,1)=ylmcst*5.*xx*zz*(xx**4 + 4.*yy**4 - 3.*yy**2*zz**2 + xx**2*(-9.*yy**2 + zz**2))
                    ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                    ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                    !xyz
                    ylmgr3_cart(1,2,3,1)=ylmcst*5.*xx*yy*(xx**2 - yy**2)*(xx**2 + yy**2 - 6.*zz**2)
                    ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                    !yyy
                    ylmgr3_cart(2,2,2,1)=ylmcst*yy*zz*(27.*xx**4 + 2.*yy**4 - 21.*yy**2*zz**2 + 12.*zz**4 + xx**2*(-41.*yy**2 + 39.*zz**2))
                    !yyz
                    ylmgr3_cart(2,2,3,1)=ylmcst*(-5.*xx**6 - 2.*yy**6 + 45.*yy**4*zz**2 - 54.*yy**2*zz**4 + 4.*zz**6 + 9.*xx**4*(2.*yy**2 + zz**2) + 3.*xx**2*(7.*yy**4 - 42.*yy**2*zz**2 + 6.*zz**4))/3.
                    ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                    ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                    !yzz
                    ylmgr3_cart(2,3,3,1)=ylmcst*yy*zz*(21.*xx**4 - 9.*yy**4 + 22.*yy**2*zz**2 - 4.*zz**4 + 6.*xx**2*(2.*yy**2 - 3.*zz**2))
                    ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                    ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                    !zzx
                    ylmgr3_cart(3,3,1,1)=ylmcst*xx*zz*(9.*xx**4 - 21.*yy**4 + 18.*yy**2*zz**2 + 4.*zz**4 - 2.*xx**2*(6.*yy**2 + 11.*zz**2))
                    ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                    ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                    !zzz
                    ylmgr3_cart(3,3,3,1)= -ylmcst*((xx**2 - yy**2)*(3.*xx**4 + 3.*yy**4 - 24.*yy**2*zz**2 + 8.*zz**4 + 6.*xx**2*(yy**2 - 4.*zz**2)))

                 else if (ll==3 .and. mm==3) then
                    ylmcst=sqrt(5./2.)*(3./4.)*sqrt(7./pi)/(rr**9)

                    !xxx
                    ylmgr3_cart(1,1,1,1)=ylmcst*12.*xx**4*(3.*yy**2 + zz**2) + (yy**2 + zz**2)**2*(11.*yy**2 + 2.*zz**2) - 3.*xx**2*(31.*yy**4 + 38.*yy**2*zz**2 + 7.*zz**4)
                    !xxy
                    ylmgr3_cart(1,1,2,1)= ylmcst*xx*yy*(-18.*xx**4 + xx**2*(89.*yy**2 + 29.*zz**2) - 3.*(11.*yy**4 + 7.*yy**2*zz**2 - 4.*zz**4))
                    ylmgr3_cart(1,2,1,:)= ylmgr3_cart(1,1,2,:)
                    ylmgr3_cart(2,1,1,:)= ylmgr3_cart(1,1,2,:)
                    !xxz
                    ylmgr3_cart(1,1,3,1)=-ylmcst*xx*zz*(6.*xx**4 + 51.*yy**4 + 57.*yy**2*zz**2 + 6.*zz**4 - xx**2*(83.*yy**2 + 23.*zz**2))
                    ylmgr3_cart(1,3,1,:)= ylmgr3_cart(1,1,3,:)
                    ylmgr3_cart(3,1,1,:)= ylmgr3_cart(1,1,3,:)
                    !xyy
                    ylmgr3_cart(1,2,2,1)=ylmcst*6.*xx**6 - 2.*yy**6 + 9.*yy**4*zz**2 + 9.*yy**2*zz**4 - 2.*zz**6 + xx**4*(-75.*yy**2 + 5.*zz**2) + xx**2*(57.*yy**4 - 36.*yy**2*zz**2 - 3.*zz**4)
                    ylmgr3_cart(2,1,2,:)= ylmgr3_cart(1,2,2,:)
                    ylmgr3_cart(1,2,2,:)= ylmgr3_cart(1,2,2,:)
                    !xyz
                    ylmgr3_cart(1,2,3,1)=-ylmcst*(yy*zz*(44.*xx**4 + 9.*yy**4 + 3.*yy**2*zz**2 - 6.*zz**4 + xx**2*(-87.*yy**2 + 3.*zz**2)))
                    ylmgr3_cart(3,2,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,1,3,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(2,3,1,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(3,1,2,:)= ylmgr3_cart(1,2,3,:)
                    ylmgr3_cart(1,3,2,:)= ylmgr3_cart(1,2,3,:)
                    !yyy
                    ylmgr3_cart(2,2,2,1)=ylmcst*xx*yy*(51.*xx**4 + 6.*yy**4 - 63.*yy**2*zz**2 + 36.*zz**4 + xx**2*(-83.*yy**2 + 87.*zz**2))
                    !yyz
                    ylmgr3_cart(2,2,3,1)=ylmcst*xx*zz*(11.*xx**4 + 36.*yy**4 - 63.*yy**2*zz**2 + 6.*zz**4 + xx**2*(-93.*yy**2 + 17.*zz**2))
                    ylmgr3_cart(2,3,2,:)= ylmgr3_cart(2,2,3,:)
                    ylmgr3_cart(3,2,2,:)= ylmgr3_cart(2,2,3,:)
                    !yzz
                    ylmgr3_cart(2,3,3,1)=ylmcst*xx*yy*(11.*xx**4 - 9.*yy**4 + 72.*yy**2*zz**2 - 24.*zz**4 + 2.*xx**2*(yy**2 - 24.*zz**2))
                    ylmgr3_cart(3,2,3,:)= ylmgr3_cart(2,3,3,:)
                    ylmgr3_cart(3,3,2,:)= ylmgr3_cart(2,3,3,:)
                    !zzx
                    ylmgr3_cart(3,3,1,1)=ylmcst*2.*xx**6 - xx**4*(13.*yy**2 + 21.*zz**2) + xx**2*(-12.*yy**4 + 90.*yy**2*zz**2 + 12.*zz**4) + 3.*(yy**6 - 3.*yy**4*zz**2 - 4.*yy**2*zz**4)
                    ylmgr3_cart(3,1,3,:)= ylmgr3_cart(3,3,1,:)
                    ylmgr3_cart(1,3,3,:)= ylmgr3_cart(3,3,1,:)
                    !zzz
                    ylmgr3_cart(3,3,3,1)= ylmcst*5.*xx*(xx**2 - 3.*yy**2)*zz*(3.*xx**2 + 3.*yy**2 - 4.*zz**2)

                 else
                    write(*,*) 'Illegal mm and ll in initylmg.'
                    write(*,*) 'll',ll,'mm',mm
                    stop
                 end if

                 ! 2-Transfer gradients into reduced coordinates
                 !ylmgr3_tmp=zero
                 !do aa=1,2
                 !   do ii=1,3
                 !      do jj=1,3
                 !         do kk=1,3
                 !            do pp=1,3
                 !               do qq=1,3
                 !                  do ss=1,3
                 !                     ylmgr3_tmp(ii,jj,kk)=ylmgr3_tmp(ii,jj,kk) &
                 !                          & +rprimd(pp,ii)*rprimd(qq,jj)*rprimd(ss,kk)*ylmgr3_cart(pp,qq,ss,aa)
                 !                  end do
                 !               end do
                 !            end do
                 !         end do
                 !      end do
                 !   end do
                 do ii=1,10
                    ia=alpha3(ii);ib=beta3(ii);igam=gamma3(ii)

                    ylm_gr(ikg+ipw,9+ii,l0+mm)=cmplx(ylmgr3_cart(ia,ib,igam,1),ylmgr3_cart(ia,ib,igam,2))
                    ylm_gr(ikg+ipw,9+ii,l0-mm)=-conjg(ylm_gr(ikg+ipw,9+ii,l0+mm))

                    !ylm_gr(ikg+ipw,9+ii,l0+mm)=cmplx(ylmgr3_cart(ia,ib,igam,1),ylmgr3_cart(ia,ib,igam,2))
                    !ylm_gr(ikg+ipw,9+ii,l0-mm)=conjg(cmplx(-1*ylmgr3_cart(ia,ib,igam,1),ylmgr3_cart(ia,ib,igam,2)))      
                    !ylmgr_red(ii,aa)=ylmgr3_tmp(ia,ib,igam)
                 end do
                 if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)
              end do !mm
                 !ylm_gr(ikg+ipw,10:19,l0+mm)=cmplx(ylmgr_red(1:10,1),ylmgr_red(1:10,2))
                 !ylm_gr(ikg+ipw,10:19,l0-mm)=onem*conjg(cmplx(ylmgr_red(1:10,1),ylmgr_red(1:10,2)))

                 !TEST: Check for NaN
!                 do ii=10,19
!                    if (real(ylm_gr(ikg+ipw,ii,l0+mm)) /= real(ylm_gr(ikg+ipw,ii,l0+mm))) then
                       !write(*,*) 'error: ylm_gr is NaN'
                       !write(*,*) 'ii',ii,'mm:',mm,'l0:',l0,'kg:',kg_k(:,ipw)
                      ! stop
!                    end if
                 !if (mm/=ll) fact=fact/real((ll+mm+1)*(ll-mm),dp)              
          ! end do ! loop over m
           end do ! loop over l

        end if !optder=8

!        End condition r<>0
     end if
      
!      End loop over k+G
    end do

!    End condition l<>0
   end if

   !************************************************************* 
   !TEST
!   open (unit=15, file='ylmgr.test', status='replace')
!   write(*,*) kptns(:,1)
!   do ipw=1,npw_k

      !      Load k+G 
!      kpg(1)=kptns(1,ikpt)+real(kg_k(1,ipw),dp)
!      kpg(2)=kptns(2,ikpt)+real(kg_k(2,ipw),dp)
!      kpg(3)=kptns(3,ikpt)+real(kg_k(3,ipw),dp)


      !      Calculate mod of k+G
!      xx=two_pi*(gprimd(1,1)*kpg(1)+gprimd(1,2)*kpg(2)+gprimd(1,3)*kpg(3))
!      yy=two_pi*(gprimd(2,1)*kpg(1)+gprimd(2,2)*kpg(2)+gprimd(2,3)*kpg(3))
!      zz=two_pi*(gprimd(3,1)*kpg(1)+gprimd(3,2)*kpg(2)+gprimd(3,3)*kpg(3))
!      rr=sqrt(xx**2+yy**2+zz**2)

      ! TEST: Use real spherical harmonics for p
      
      !Ylm:
!      ylm(ikg+ipw,2)=cmplx(sqrt(3./four_pi)*yy/rr,zero)
!      ylm(ikg+ipw,3)=cmplx(sqrt(3./four_pi)*zz/rr,zero)
!      ylm(ikg+ipw,4)=cmplx(sqrt(3./four_pi)*xx/rr,zero)
      
      !dYlm/dKz
!      ylm_gr(ikg+ipw,3,2)=cmplx(-sqrt(3./four_pi)*yy*zz/rr**3,zero)
!      ylm_gr(ikg+ipw,3,3)=cmplx(sqrt(3./four_pi)*(xx**2+yy**2)/rr**3,zero)
!      ylm_gr(ikg+ipw,3,4)=cmplx(-sqrt(3./four_pi)*xx*zz/rr**3,zero)

      !d2Ylm/d2Kz
!      ylm_gr(ikg+ipw,9,2)=cmplx(-sqrt(3./four_pi)*yy*(xx**2+yy**2-2*zz**2)/rr**5,zero)
!      ylm_gr(ikg+ipw,9,3)=cmplx(-3.*sqrt(3./four_pi)*zz*(xx**2+yy**2)/rr**5,zero)
!      ylm_gr(ikg+ipw,9,4)=cmplx(-sqrt(3./four_pi)*xx*(xx**2+yy**2-2*zz**2)/rr**5,zero)

      !d3Ylm/d3Kz
!      ylm_gr(ikg+ipw,19,2)=cmplx(3.*sqrt(3./four_pi)*yy*zz*(3.*(xx**2+yy**2)-2*zz**2)/rr**7,zero)
!      ylm_gr(ikg+ipw,19,3)=cmplx(-3.*sqrt(3./four_pi)*(xx**2+yy**2)*(xx**2+yy**2-4.*zz**2)/rr**7,zero)
!      ylm_gr(ikg+ipw,19,4)=cmplx(3.*sqrt(3./four_pi)*xx*zz*(3.*(xx**2+yy**2)-2*zz**2)/rr**7,zero)

      


!      write(15,'(9e16.4e2)') xx,yy,zz,ylm(ipw,2),ylm(ipw,3),ylm(ipw,4)!,ylm_gr(ipw,19,2),ylm_gr(ipw,19,3),ylm_gr(ipw,19,4)


!   end do
!   close(unit=15)
!   stop
   !*************************************************************

   ABI_DEALLOCATE(kg_k)

   ikg=ikg+npw_k
 end do !  End Loop over k-points

!Release the temporary memory
!Allocate some memory
!TEST
!write(*,*)'ylm(20,2):',ylm(57,4)
!write(*,*)'ylm_gr(20,1,2):',ylm_gr(57,1,4)
!write(*,*)'ylm_gr(20,4,2):',ylm_gr(57,4,4)
!write(*,*)'ylm_gr(20,10,2):',ylm_gr(57,10,4)
!stop

!write(*,*) 'here2'
! if (optder/=0) then
!   ABI_DEALLOCATE(ylmgr_cart)
! end if
!TEST
!write(*,*) 'here3'
! if (optder/=0.and.optder/=2) then
!   ABI_DEALLOCATE(ylmgr_red)
! end if
!TEST
!write(*,*) 'here4'
! if (optder>=2) then
!   ABI_DEALLOCATE(ylmgr2_cart)
!   ABI_DEALLOCATE(ylmgr2_tmp)
!   ABI_DEALLOCATE(ylmgr_red)
!   ABI_DEALLOCATE(blm)
! end if

end subroutine cinitylmgi3
!!***
