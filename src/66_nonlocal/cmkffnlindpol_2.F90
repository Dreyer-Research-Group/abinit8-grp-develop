!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkffnl
!! NAME
!! mkffnl
!!
!! FUNCTION
!! Make FFNL, nonlocal form factors, for each type of atom up to ntypat
!! and for each angular momentum.
!! When Legendre polynomials are used in the application of the
!!   nonlocal operator, FFNLs depend on (l,n) components; in this
!!   case, form factors are real and divided by |k+G|^l;
!! When spherical harmonics are used, FFNLs depend on (l,m,n)
!!   components; in this case, form factors are multiplied by Ylm(k+G).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (DCA, XG, GMR, MT, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=described below
!!
!!
!! NOTES
!!
!! CEDrev: FULL CARTESIAN VERSION
!!
!! TODO
!!  Some parts can be rewritten with BLAS1 calls.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cmkffnlindpol_2(calcpm,dtset,ffnl,iat,kg,kpt,mpi_enreg,npw,Psps,sign_dyad)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_splines
 use defs_datatypes
 use defs_abitypes
 use m_geometry

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cmkffnlindpol_2'
 use interfaces_18_timing
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw,calcpm,iat
! type(crystal_t),intent(in) :: Cryst
 type(pseudopotential_type),intent(in) :: Psps
 type(dataset_type), intent(in) :: dtset
 type(MPI_type), intent(in) :: mpi_enreg
! type(datafiles_type), intent(in) :: dtfil
!arrays
 integer,intent(in) :: kg(3,npw)
 integer,intent(out) :: sign_dyad(Psps%lnmax,Psps%ntypat)
 real(dp),intent(in) :: kpt(3) 
 complex(dp),intent(out) :: ffnl(npw,20,Psps%lnmax**2)
!type(kb_potential),intent(out) :: KBgrad_k

!Local variables-------------------------------
!scalars
 integer :: iffnl,ig,ig0,il,ilm,ilmn,iln,iln0,im,itypat,mu,mua,mub,muc,nlmn,imm
 real(dp),parameter :: renorm_factor=0.5d0/pi**2,tol_norm=tol10
 real(dp) :: ecut,ecutsm,effmass,fact,kpg1,kpg2,kpg3,kpgc1,kpgc2,kpgc3,rmetab,yp1,factor
 logical :: testnl=.false.
 character(len=500) :: message
!arrays

!CEDrev: FOR TESTING, NEED TO CHANGE BACK WHEN USING INDPOL
! integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer,parameter :: alpha(6)=(/1,1,1,2,2,3/),beta(6)=(/1,2,3,2,3,3/)
 integer,parameter :: alpha3(10)=(/1,1,1,1,1,1,2,2,2,3/)
 integer,parameter :: beta3(10)= (/1,1,1,2,2,3,2,2,3,3/)
 integer,parameter :: gamma3(10)=(/1,2,3,2,3,3,2,3,3,3/)
! integer,allocatable :: sign_dyad(:,:)

 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),ucvol,gam(3,3,npw),rprimd(3,3),rred(3,dtset%natom),rcart(3)
 real(dp) :: tsec(2)
 real(dp) :: xdotg(npw),gcart(3,npw)
 real(dp),allocatable :: kpgc(:,:),kpgn(:,:),kpgnorm(:),kpgnorm_inv(:),wk_ffnl1(:)
 real(dp),allocatable :: wk_ffnl2(:),wk_ffnl3(:),wk_ffspl(:,:),wk_ffnl11(:),wk_ffnl12(:),wk_ffnl13(:),wk_ffspl1(:,:)
 complex(dp) :: sfac(npw),ffnl2(3,3),ylm2(3,3)
 complex(dp) :: ffnlpm(npw,20,Psps%lnmax**2,dtset%natom)
 complex(dp),allocatable :: ylm(:,:),ylm_gr(:,:,:)

!CEDrev:
logical :: filexist

! *************************************************************************

! DBG_ENTER("COLL")

!Keep track of time spent in mkffnl
 call timab(16,1,tsec)

! Get metric, etc.
call mkrdim(dtset%acell_orig(1:3,1),dtset%rprim_orig(1:3,1:3,1),rprimd)
call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Get complex spherical harmonics:
ABI_ALLOCATE(ylm,(npw,Psps%mpsang*Psps%mpsang))
ABI_ALLOCATE(ylm_gr,(npw,27,Psps%mpsang*Psps%mpsang))
call cinitylmgi3(gprimd,kg,kpt,1,mpi_enreg,Psps%mpsang,npw,&
&  dtset%nband(1),1,(/npw/),dtset%nsppol,8,rprimd,ylm,ylm_gr)


! Get dyadic sign
do itypat=1,Psps%ntypat
   nlmn=count(Psps%indlmn(3,:,itypat)>0)
   iln0=0
   do ilmn=1,nlmn
      iln=Psps%indlmn(5,ilmn,itypat)
      if (iln>iln0) then
         iln0=iln
         sign_dyad(iln,itypat)=nint(DSIGN(one,Psps%ekb(ilmn,itypat)))
      end if
   end do
end do

!Get (k+G) and |k+G|:
 ABI_ALLOCATE(kpgnorm,(npw))
 ABI_ALLOCATE(kpgnorm_inv,(npw))
 ig0=-1 ! index of |k+g|=0 vector
 ABI_ALLOCATE(kpgc,(npw,3))
 ABI_ALLOCATE(kpgn,(npw,3))
 do ig=1,npw
    kpg1=kpt(1)+dble(kg(1,ig))
    kpg2=kpt(2)+dble(kg(2,ig))
    kpg3=kpt(3)+dble(kg(3,ig))

    ! leave out 2pi here for splfit, then add in to gam and kpgnorm_inv in the second order term
    kpgc1=kpg1*gprimd(1,1)+kpg2*gprimd(1,2)+kpg3*gprimd(1,3)
    kpgc2=kpg1*gprimd(2,1)+kpg2*gprimd(2,2)+kpg3*gprimd(2,3)
    kpgc3=kpg1*gprimd(3,1)+kpg2*gprimd(3,2)+kpg3*gprimd(3,3)
    kpgc(ig,1)=kpgc1
    kpgc(ig,2)=kpgc2
    kpgc(ig,3)=kpgc3
    kpgnorm(ig)=sqrt(kpgc1*kpgc1+kpgc2*kpgc2+kpgc3*kpgc3)
    if (kpgnorm(ig)<=tol_norm) ig0=ig
    kpgnorm_inv(ig)=1.d0/max(kpgnorm(ig),tol_norm)

    !All cartesian, all the time
    kpgn(ig,1)=kpgc1*kpgnorm_inv(ig)
    kpgn(ig,2)=kpgc2*kpgnorm_inv(ig)
    kpgn(ig,3)=kpgc3*kpgnorm_inv(ig)

    ! For structure factor
    gcart(:,ig)=real(kg(1,ig))*gprimd(:,1)+real(kg(2,ig))*gprimd(:,2)+real(kg(3,ig))*gprimd(:,3)

    ! For second derivative:
    do mua=1,3
       do mub=1,3
          !If EVERYTHING is cartesian, rmet is delta function
          if (mua==mub) then
             gam(mua,mub,ig)=(one-kpgn(ig,mua)*kpgn(ig,mub))*kpgnorm_inv(ig)/two_pi  
          else             
             gam(mua,mub,ig)=(-kpgn(ig,mua)*kpgn(ig,mub))*kpgnorm_inv(ig)/two_pi
          end if
       end do
    end do
  
 end do

! Arrays for radial part of ff and RADIAL derivatives 
 ABI_ALLOCATE(wk_ffnl1,(npw))
 ABI_ALLOCATE(wk_ffnl2,(npw))
 ABI_ALLOCATE(wk_ffnl3,(npw))
 ABI_ALLOCATE(wk_ffspl,(Psps%mqgrid_ff,2))
 ABI_ALLOCATE(wk_ffnl11,(npw))
 ABI_ALLOCATE(wk_ffnl12,(npw))
 ABI_ALLOCATE(wk_ffnl13,(npw))
 ABI_ALLOCATE(wk_ffspl1,(Psps%mqgrid_ff,2))

! NO LOOP OVER ATOMS!
! do ia=1,dtset%natom
    itypat=dtset%typat(iat)

    ! Structure factor
    rcart(:)=matmul(rprimd(:,:),dtset%xred_orig(:,iat,1))
    xdotg(:)=two_pi*(gcart(1,:)*rcart(1)+gcart(2,:)*rcart(2)+gcart(3,:)*rcart(3))
    sfac(:)=cmplx(cos(-1.*xdotg(:)),sin(-1.*xdotg(:)))

    !  Loop over (l,m,n) values
    iln0=0;nlmn=count(Psps%indlmn(3,:,itypat)>0)
    ! zero index for ffnl
    iffnl=0
    do ilmn=1,nlmn
       il=Psps%indlmn(1,ilmn,itypat)+1
       !     im=Psps%indlmn(2,ilmn,itypat)
       ilm =Psps%indlmn(4,ilmn,itypat)
       iln =Psps%indlmn(5,ilmn,itypat)

       factor=four_pi*sqrt(abs(Psps%ekb(iln,itypat))/ucvol)

       ! 884rev: For now I have removed the capability for third derivatives so I do not need 
       ! to modify the psp read in. Commented lines labeled (*)
       
       ! Store form factors (from ffspl)
       do ig=1,Psps%mqgrid_ff
          wk_ffspl(ig,1)=PsPs%ffspl(ig,1,iln,itypat)
          wk_ffspl(ig,2)=Psps%ffspl(ig,2,iln,itypat)
          ! For third derivative, 1st and 3rd of ff
          !wk_ffspl1(ig,1)=PsPs%ffspl1(ig,1,iln,itypat) (*)
          !wk_ffspl1(ig,2)=Psps%ffspl1(ig,2,iln,itypat) (*)
       end do
       call splfit(Psps%qgrid_ff,wk_ffnl2,wk_ffspl,1,kpgnorm,wk_ffnl1,Psps%mqgrid_ff,npw)
       call splfit(Psps%qgrid_ff,wk_ffnl3,wk_ffspl,2,kpgnorm,wk_ffnl1,Psps%mqgrid_ff,npw)
!       call splfit(Psps%qgrid_ff,wk_ffnl12,wk_ffspl1,1,kpgnorm,wk_ffnl11,Psps%mqgrid_ff,npw) (*)
!       call splfit(Psps%qgrid_ff,wk_ffnl13,wk_ffspl1,2,kpgnorm,wk_ffnl11,Psps%mqgrid_ff,npw) (*)
      

       ! To account for the 1/2pi in d/dK
       wk_ffnl2=wk_ffnl2/two_pi
       wk_ffnl3=wk_ffnl3/two_pi**2
       !wk_ffnl13=wk_ffnl13/two_pi**3 (*)

       ! imm is index for ylm
       imm=(il-1)*(il-1)
       ! Loop over m
       do im=1,2*(il-1)+1

          ! Indicies
          iffnl=iffnl+1
          imm=imm+1
          
          !TEST: Check the indicies
          !write(*,*) 'MKFFNL il',il,'ilm', ilm,'iln',iln,'iffnl',iffnl, 'im',im, 'imm',imm

          do ig=1,npw
             ! Direction of derivatives:
             do mu=1,10
                
                ! Potential (only do once):
                if (mu==1) ffnl(ig,1,iffnl)=factor*sfac(ig)*ylm(ig,imm)*wk_ffnl1(ig)*sign_dyad(iln,itypat)
                
                ! First derivative (three components):
                if (mu <= 3) ffnl(ig,mu+1,iffnl)=factor*sfac(ig)*(ylm(ig,imm)*wk_ffnl2(ig)*kpgn(ig,mu)+ylm_gr(ig,mu,imm)*wk_ffnl1(ig))
                
                ! Second derivative (6 components)
                if (mu <= 6) then
                   mua=alpha(mu); mub=beta(mu)
                   
                   ffnl(ig,mu+4,iffnl)=factor*sfac(ig)*(&
                        &   ylm_gr(ig,3+mu,imm)*wk_ffnl1(ig) &
                        & + gam(mua,mub,ig)*ylm(ig,imm)*wk_ffnl2(ig) &
                        & + ylm(ig,imm)*kpgn(ig,mua)*kpgn(ig,mub)*wk_ffnl3(ig) &
                        & + ylm_gr(ig,mua,imm)*kpgn(ig,mub)*wk_ffnl2(ig) &
                        & + ylm_gr(ig,mub,imm)*kpgn(ig,mua)*wk_ffnl2(ig))

                end if
                ylm2(1,1)=ylm_gr(ig,4,imm)
                ylm2(1,2)=ylm_gr(ig,5,imm)
                ylm2(2,1)=ylm_gr(ig,5,imm)
                ylm2(1,3)=ylm_gr(ig,6,imm)
                ylm2(3,1)=ylm_gr(ig,6,imm)
                ylm2(2,2)=ylm_gr(ig,7,imm)
                ylm2(2,3)=ylm_gr(ig,8,imm)
                ylm2(3,2)=ylm_gr(ig,8,imm)
                ylm2(3,3)=ylm_gr(ig,9,imm)

! (*) No third derivative for now, see above
                ffnl(ig,mu+10,iffnl)=zero
!                mua=alpha3(mu); mub=beta3(mu); muc=gamma3(mu)

!                ffnl(ig,mu+10,iffnl) = factor*sfac(ig)*( &
!                     &  kpgn(ig,mua)*kpgn(ig,mub)*kpgn(ig,muc)*ylm(ig,imm)*wk_ffnl13(ig) & !1
!                     & +ylm_gr(ig,mu+9,imm)*wk_ffnl1(ig) & !4
!                     !
!                     & +kpgn(ig,mua)*kpgn(ig,mub)*ylm_gr(ig,muc,imm)*wk_ffnl3(ig) & !1
!                     & +kpgn(ig,mua)*kpgn(ig,muc)*ylm_gr(ig,mub,imm)*wk_ffnl3(ig) & !2
!                     & +kpgn(ig,mub)*kpgn(ig,muc)*ylm_gr(ig,mua,imm)*wk_ffnl3(ig) & !3
!                     !
!                     & +kpgn(ig,mua)*ylm2(mub,muc)*wk_ffnl2(ig) & !2
!                     & +kpgn(ig,muc)*ylm2(mua,mub)*wk_ffnl2(ig) & !4
!                     & +kpgn(ig,mub)*ylm2(mua,muc)*wk_ffnl2(ig) & !3
!                     ! Remember, this is all in cartesian coordinates
!                     & +gam(mua,muc,ig)*kpgn(ig,mub)*ylm(ig,imm)*wk_ffnl3(ig) & !1
!                     & +gam(mub,muc,ig)*kpgn(ig,mua)*ylm(ig,imm)*wk_ffnl3(ig) & !1
!                     & +gam(mua,mub,ig)*kpgn(ig,muc)*ylm(ig,imm)*wk_ffnl3(ig) & !5
!                     !
!                     & +gam(mua,muc,ig)*ylm_gr(ig,mub,imm)*wk_ffnl2(ig) & !2
!                     & +gam(mub,muc,ig)*ylm_gr(ig,mua,imm)*wk_ffnl2(ig) & !3
!                     & +gam(mua,mub,ig)*ylm_gr(ig,muc,imm)*wk_ffnl2(ig) & !5
!                     ! WHY 1/2pi????
!                     & -gam(mua,muc,ig)*kpgn(ig,mub)*ylm(ig,imm)*wk_ffnl2(ig)*kpgnorm_inv(ig)/two_pi & !5.2
!                     & -gam(mub,muc,ig)*kpgn(ig,mua)*ylm(ig,imm)*wk_ffnl2(ig)*kpgnorm_inv(ig)/two_pi & !5.3
!                     & -gam(mua,mub,ig)*kpgn(ig,muc)*ylm(ig,imm)*wk_ffnl2(ig)*kpgnorm_inv(ig)/two_pi) !5


             end do !mu
             
          end do !ig          
       end do !im
    end do !ilmn
! end do !itypat

 call timab(16,2,tsec)

! DBG_EXIT("COLL")

end subroutine cmkffnlindpol_2
!!***
