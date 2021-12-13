!{\src2tex{textfont=tt}}
!!****f* ABINIT/joper
!! NAME
!! joper 
!!
!! FUNCTION
!! Applies the local and nonlocal expansion up to 2nd order  of the 
!! current density operator to a wavefunction
!!
!! Written by Cyrus Dreyer, Rutgers 2017
!! From sternab3, joper_sign_ffnl.F90
!!
!! INPUTS
!!  calcnl = 1 to calc NL part, 0 to just to local
!!  cg = GS wavefunction pw coefs
!!  cg1 = AD wavefunction pw coefs
!!  cryst <type(crystal)>=Unit cell and symmetries
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  kg = GS wavefunction G vectors
!!  kg = AD wavefunction G vectors 
!!  psps = Psuedopotential info  
!! 
!!
!!
!!
!!
!! OUTPUT
!!
!! PARENTS
!! 
!!
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine joper(calcnl,cwave0_npw1,dtset,gprimd,irfdir,kg1_k,kpt,mpi_enreg,npw1,psps,qpc,rprimd)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 !use m_crystal,    only : crystal_init, crystal_free, crystal_t,isalchemical
 use m_pawcprj, only : pawcprj_type
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'joper'
 use interfaces_18_timing
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(dataset_type), intent(in), target :: dtset
 type(pseudopotential_type), intent(in) :: psps
 type(MPI_type), intent(in) :: mpi_enreg
 integer,intent(in) :: npw1,calcnl,irfdir
 real(dp), intent(inout) :: cwave0_npw1(2,npw1*dtset%nspinor) ! NCrev
 real(dp), intent(in) :: gprimd(3,3),rprimd(3,3),kpt(3),qpc(3)
 integer, intent(in) :: kg1_k(3,npw1)


!Local variables-------------------------------
 real(dp) :: gsg(2,npw1*dtset%nspinor),gsgd(3,2,npw1*dtset%nspinor)
 real(dp) :: kplusg(3),qpc_(3)
 real(dp) :: cpuin,wallin,cpuout,wallout
 complex(dp) :: vnlg0(3,npw1*dtset%nspinor), vnlg1(3,npw1*dtset%nspinor),vnlg2(3,npw1*dtset%nspinor)
 complex(dp) :: fcta1,fcta2(3),fcta3(3),fcta4,fcta5(3,3),fcta6(3,3),fcta7(3,3,3),fcta8(3,3,3)
 complex(dp) :: ffnl1(3),ffnl2(3,3),ffnl3(3,3,3)
 complex(dp) :: ffnl(npw1*dtset%nspinor,20,psps%lnmax*psps%lnmax)
 integer :: sign_dyad(psps%lnmax,psps%ntypat)
 integer, allocatable :: nlmn(:),iln(:,:)
 integer :: ikg,ix,ix1,ix2,iat,ilm,il,ityp,nln,fid,tr,ispinor
 integer, parameter :: alpha(6)=(/1,1,1,2,2,3/),beta(6)=(/1,2,3,2,3,3/)
 integer, parameter :: alpha3(10)=(/1,1,1,1,1,1,2,2,2,3/)
 integer, parameter :: beta3(10)= (/1,1,1,2,2,3,2,2,3,3/)
 integer, parameter :: gamma3(10)=(/1,2,3,2,3,3,2,3,3,3/)
!------------------------------------------------

!************************************************
! Setup
!************************************************

! 884rev: 
!
! NOTE, for now third k derivative (order q^2) is not implemented 
! so psp read in did not have to be modified. Commented out lines labelled with (*)
!
! NCrev: edits for noncollinear

! if (mpi_enreg%me_kpt==0) then
!    if (calcnl==0) then
!       write(*,*)'LOCAL JOPER will be applied in direction',irfdir
!    else if (calcnl==2) then
!       write(*,*)'NONLOCAL JOPER will be applied in direction',irfdir
!    else if (calcnl==1) then
!       write(*,*)'DIAMAG JOPER will be applied in direction',irfdir
!    end if
! end if

!timing
call timein(cpuin,wallin)

!TESTs
tr=0 !Only TR systyems right now
!sign_dyadj(:,:)=sign_dyad(:,:) !Test removing dyadic sign
!write(*,*) 'joper WITHOUT SD, really'

!************************************************
! Local part
!************************************************
  
! Apply local part of current density operator (momentum) to wavefunction
do ispinor=1,dtset%nspinor ! NCrev 
   do ikg=1,npw1

      !NCrev
      gsg(:,ikg+npw1*(ispinor-1))=cwave0_npw1(:,ikg+npw1*(ispinor-1))
      !gsg(:,ikg)=cwave0_npw1(:,ikg)

      if (calcnl<2) then
         kplusg(:)=two_pi*matmul(gprimd(:,:),kpt(:)+real(kg1_k(:,ikg)))+0.5*qpc(:)

         ! NCrev and shortened
         gsgd(:,1,ikg+npw1*(ispinor-1))=kplusg(:)*gsg(1,ikg+npw1*(ispinor-1))
         gsgd(:,2,ikg+npw1*(ispinor-1))=kplusg(:)*gsg(2,ikg+npw1*(ispinor-1))

         !gsgd(1,1,ikg)=kplusg(1)*gsg(1,ikg)
         !gsgd(1,2,ikg)=kplusg(1)*gsg(2,ikg)
         !gsgd(2,1,ikg)=kplusg(2)*gsg(1,ikg)
         !gsgd(2,2,ikg)=kplusg(2)*gsg(2,ikg)
         !gsgd(3,1,ikg)=kplusg(3)*gsg(1,ikg)
         !gsgd(3,2,ikg)=kplusg(3)*gsg(2,ikg)

         !    write(*,'(3i5,7e16.5e2)') kg1_k(:,ikg),gsg(:,ikg),kplusg(:),gsgd(irfdir,:,ikg)
      end if
   end do
end do
 ! stop
 !************************************************
 ! Nonlocal potential
 !************************************************

 if (calcnl>0) then
    
    vnlg0=zero;vnlg1=zero;vnlg2=zero;

    ! determine number of nonlocal projectors to loop over  
    ABI_ALLOCATE(nlmn,(dtset%natom))
    ABI_ALLOCATE(iln,(dtset%natom,Psps%lnmax*Psps%lnmax))
    nlmn(:)=0
    do iat=1,dtset%natom
       ityp=dtset%typat(iat)
       nln=count(Psps%indlmn(3,:,ityp)>0)
       do ilm=1,nln
          il=Psps%indlmn(1,ilm,ityp)
          !Need this for sign_dyad
          iln(iat,nlmn(iat)+1:nlmn(iat)+2*il+1)=Psps%indlmn(5,ilm,ityp)
          nlmn(iat)= nlmn(iat)+2*il+1
       end do
    end do

    !TEST
    !write(*,*) 'joper nlmn',nlmn(1)
    !stop

    !loop over atoms and ang momentum channels
    do iat=1,dtset%natom
       ityp=dtset%typat(iat)

       ! Calculate ffnl here. This takes more time, but avoinds storing a 
       ! very large array which is prohibative for large cells
       ! TODO: ALSO CALC FFNL ON THE FLY FOR iln!!!!
       ffnl=zero
       call cmkffnlindpol_2(0,dtset,ffnl,iat,kg1_k,kpt,mpi_enreg,npw1,psps,sign_dyad)

       ! Loop over angular momentum channels
       do ilm=1,nlmn(iat)
 
          ! =====Calculate overlaps=====
          ! For now have them all. Since we have TR, we only need half.
          fcta1=0;fcta2=0;fcta3=0;fcta4=0;fcta5=0;fcta6=0;fcta7=0;fcta8=0

          if (tr==1) then
             do ikg=1,npw1

                ! For order in q >=0
                fcta1=fcta1+conjg(cmplx(gsg(1,ikg),gsg(2,ikg)))*ffnl(ikg,1,ilm)
                fcta3(:)=fcta3(:)+conjg(cmplx(gsg(1,ikg),gsg(2,ikg))) &
                     & *ffnl(ikg,2:4,ilm)

                ! For order in q >=1
                fid=4
                do ix=1,3 ! alpha
                   do ix1=ix,3 ! beta
                      fid=fid+1

                      fcta5(ix,ix1)=fcta5(ix,ix1)+conjg(cmplx(gsg(1,ikg),gsg(2,ikg))) &
                           & *ffnl(ikg,fid,ilm)
                      fcta5(ix1,ix)= fcta5(ix,ix1)

                   end do
                end do
                ! For order in q >=2
                fid=10
                do ix=1,3
                   do ix1=ix,3
                      do ix2=ix1,3
                         fid=fid+1

                         fcta7(ix,ix1,ix2)=fcta7(ix,ix1,ix2)+conjg(cmplx(gsg(1,ikg),gsg(2,ikg))) &
                              & *ffnl(ikg,fid,ilm)
                         fcta7(ix,ix2,ix1)=fcta7(ix,ix1,ix2)
                         fcta7(ix1,ix,ix2)=fcta7(ix,ix1,ix2)
                         fcta7(ix1,ix2,ix)=fcta7(ix,ix1,ix2)
                         fcta7(ix2,ix,ix1)=fcta7(ix,ix1,ix2)
                         fcta7(ix2,ix1,ix)=fcta7(ix,ix1,ix2)

                      end do
                   end do
                end do
             end do !ikg
          end if !tr

         do ikg=1,npw1

            ! For order in q >=0
            fcta2(:)=fcta2(:)+conjg(ffnl(ikg,2:4,ilm)) &
                 & *cmplx(gsg(1,ikg),gsg(2,ikg))
            fcta4=fcta4+conjg(ffnl(ikg,1,ilm))*cmplx(gsg(1,ikg),gsg(2,ikg))

            ! For order in q >=1
            fid=4
            do ix=1,3 ! alpha
               do ix1=ix,3 ! beta
                  fid=fid+1

                  fcta6(ix,ix1)=fcta6(ix,ix1)+conjg(ffnl(ikg,fid,ilm)) &
                       & *cmplx(gsg(1,ikg),gsg(2,ikg))
                  fcta6(ix1,ix)=fcta6(ix,ix1)

               end do
            end do
            ! For order in q >=2
            ! NO 3rd DERIV (*)
            !fid=10
            !do ix=1,3
            !   do ix1=ix,3
            !      do ix2=ix1,3
            !         fid=fid+1

             !        fcta8(ix,ix1,ix2)=fcta8(ix,ix1,ix2)+conjg(ffnl(ikg,fid,ilm)) &
             !             & *cmplx(gsg(1,ikg),gsg(2,ikg))
             !        fcta8(ix,ix2,ix1)=fcta8(ix,ix1,ix2)
             !        fcta8(ix1,ix,ix2)=fcta8(ix,ix1,ix2)
             !        fcta8(ix1,ix2,ix)=fcta8(ix,ix1,ix2)
             !        fcta8(ix2,ix,ix1)=fcta8(ix,ix1,ix2)
             !        fcta8(ix2,ix1,ix)=fcta8(ix,ix1,ix2)

             !     end do
             !  end do
             !end do
         end do !ikg

          
          !TEST: Check for NaN
          !if (fcta4 /= fcta4) then
          !   write(*,*) 'ilm',ilm,'iat',iat,'fta1',fcta1
          !   stop
          !end if
          !do ix=1,3
          !   if (fcta2(ix) /= fcta2(ix)) then
          !      write(*,*) 'ilm',ilm,'iat',iat,'fta3',fcta3
          !      stop
          !   end if
          !   do ix1=1,3
          !      if (fcta6(ix,ix1) /= fcta6(ix,ix1)) then
          !         write(*,*) 'ilm',ilm,'iat',iat,'fta5',fcta5
          !         stop
          !      end if
          !      do ix2=1,3
          !         if (fcta8(ix,ix1,ix2) /= fcta8(ix,ix1,ix2)) then
          !            write(*,*) 'ilm',ilm,'iat',iat,'fcta7',fcta7
          !            stop
          !         end if
          !      end do
          !   end do
          !end do
          
          !write(*,*) 'No NaN in ftas'
          !write (*,'(a15,6e12.2e2)') 'joper fcta',fcta5(3,:)
          !write (*,'(a10,6e16.4e2)') 'joper fcta3',fcta3
          !write (*,'(a10,2e16.4e2)') 'fcta5',fcta5
          !write (*,'(a10,2e16.4e2)') 'fcta7',fcta7
          
          
          ! Correction term expansion
          do ikg=1,npw1

             ! Convert ffnl to explict x,y,z indicies
             do ix=1,3
                ffnl1(ix)=ffnl(ikg,ix+1,ilm)
              end do
             do ix=1,6
                ffnl2(alpha(ix),beta(ix))=ffnl(ikg,ix+4,ilm)
                ffnl2(beta(ix),alpha(ix))=ffnl(ikg,ix+4,ilm)
             end do
           
             ! NO 3rd DERIV (*)
             !do ix=1,10
             !   ffnl3(alpha3(ix),beta3(ix),gamma3(ix))=ffnl(ikg,ix+10,ilm)
             !   ffnl3(alpha3(ix),gamma3(ix),beta3(ix))=ffnl(ikg,ix+10,ilm)
             !   ffnl3(beta3(ix),alpha3(ix),gamma3(ix))=ffnl(ikg,ix+10,ilm)
             !   ffnl3(beta3(ix),gamma3(ix),alpha3(ix))=ffnl(ikg,ix+10,ilm)
             !   ffnl3(gamma3(ix),alpha3(ix),beta3(ix))=ffnl(ikg,ix+10,ilm)
             !   ffnl3(gamma3(ix),beta3(ix),alpha3(ix))=ffnl(ikg,ix+10,ilm)
             !end do

             ! zeroth order nl induced polarization
             vnlg0(:,ikg)=vnlg0(:,ikg)-4*( &
                  & ffnl(ikg,1,ilm)*fcta2(:)+ffnl(ikg,2:4,ilm)*fcta4)

             do ix=1,3 ! alpha
                do ix1=1,3 !beta
                   
                   !ICL path
                   if (dtset%useric==0) then
                      vnlg1(ix,ikg)=vnlg1(ix,ikg)-2.*qpc(ix1)*( &
                           &  ffnl1(ix)*fcta2(ix1)*sign_dyad(iln(iat,ilm),ityp) &
                           & +ffnl1(ix1)*fcta2(ix)*sign_dyad(iln(iat,ilm),ityp) &
                           & +ffnl(ikg,1,ilm)*fcta6(ix,ix1) &
                           & +ffnl2(ix,ix1)*fcta4)
                     
                      !PM path
                   else if (dtset%useric==1) then
                      vnlg1(ix,ikg)=vnlg1(ix,ikg)-2.*qpc(ix1)*( &
                           &  2.*ffnl1(ix)*fcta2(ix1)*sign_dyad(iln(iat,ilm),ityp) &
                           & +ffnl(ikg,1,ilm)*fcta6(ix,ix1) &
                           & +ffnl2(ix,ix1)*fcta4)

                   end if

                   
                   ! NO 3rd DERIV (*)
                   !do ix2=1,3 ! gamma
                   !   if (dtset%useric==0) then
                   !      ! ICL path
                   !      vnlg2(ix,ikg)=vnlg2(ix,ikg)-(2./3.)*qpc(ix1)*qpc(ix2)*( &
                   !           &  ffnl(ikg,1,ilm)*fcta8(ix,ix1,ix2) &
                   !           & +ffnl1(ix)*fcta6(ix1,ix2)*sign_dyad(iln(iat,ilm),ityp) &
                   !           & +ffnl1(ix1)*fcta6(ix2,ix)*sign_dyad(iln(iat,ilm),ityp) &
                   !           & +ffnl1(ix2)*fcta6(ix,ix1)*sign_dyad(iln(iat,ilm),ityp) &
                   !           & +ffnl2(ix,ix1)*fcta2(ix2)*sign_dyad(iln(iat,ilm),ityp) &
                   !           & +ffnl2(ix,ix2)*fcta2(ix1)*sign_dyad(iln(iat,ilm),ityp) &
                   !           & +ffnl2(ix1,ix2)*fcta2(ix)*sign_dyad(iln(iat,ilm),ityp) &
                   !           & +ffnl3(ix,ix1,ix2)*fcta4)            
                   !      ! PM path
                   !   else if (dtset%useric==1) then
                   !      vnlg2(ix,ikg)=vnlg2(ix,ikg)-(2./3.)*qpc(ix1)*qpc(ix2)*( &
                   !           &  ffnl(ikg,1,ilm)*fcta8(ix,ix1,ix2) &
                   !           & +3.*ffnl1(ix)*fcta6(ix1,ix2)*sign_dyad(iln(iat,ilm),ityp) &
                   !           & +3.*ffnl2(ix,ix2)*fcta2(ix1)*sign_dyad(iln(iat,ilm),ityp) &
                   !           & +ffnl3(ix,ix1,ix2)*fcta4)
                   !   end if
                   !end do !ix2
                end do !ix1
             end do ! ix

          end do !ikg
       end do !ilm
    end do !iat

 end if !calcnl
 
!TEST: JUST CHECK VNLG0
!write(*,*) 'JUST vnlg0 and vnlg1!!!!!'
!vnlg1=zero
 !vnlg2=zero
 
 
 !Sum contributions. For now just in rfdir direction
 cwave0_npw1=zero
 do ispinor=1,dtset%nspinor !NCrev
    do ikg=1,npw1
       !TEST:
       !if (calcnl<2) cwave0_npw1(:,ikg)=-gsgd(irfdir,:,ikg) 
       !if (calcnl<2) cwave0_npw1(:,ikg)=-4*gsgd(irfdir,:,ikg)
       !NCrev
       if (calcnl<2) cwave0_npw1(:,ikg+npw1*(ispinor-1))=-4.*gsgd(irfdir,:,ikg+npw1*(ispinor-1))

       if (calcnl>0) then
          !cwave0_npw1(1,ikg)=cwave0_npw1(1,ikg) &
          !     & +real(real(vnlg0(irfdir,ikg)+vnlg1(irfdir,ikg)))!+vnlg2(irfdir,ikg))) (*)
          !cwave0_npw1(2,ikg)=cwave0_npw1(2,ikg) &
          !     & +real(aimag(vnlg0(irfdir,ikg)+vnlg1(irfdir,ikg)))!+vnlg2(irfdir,ikg))) (*)

          !NCrev
          cwave0_npw1(1,ikg+npw1*(ispinor-1))=cwave0_npw1(1,ikg+npw1*(ispinor-1)) &
               & +real(real(vnlg0(irfdir,ikg)+vnlg1(irfdir,ikg)))!+vnlg2(irfdir,ikg))) (*)
          cwave0_npw1(2,ikg+npw1*(ispinor-1))=cwave0_npw1(2,ikg+npw1*(ispinor-1)) &
               & +real(aimag(vnlg0(irfdir,ikg)+vnlg1(irfdir,ikg)))!+vnlg2(irfdir,ikg))) (*)
       end if
    end do
 end do
!TEST
!write(*,*) "4* TEST"


!timing

!call timab(101,4,tottim)
call timein(cpuout,wallout)
!write(*,*) 'time in joper:',wallout-wallin
!if (mpi_enreg%me_kpt==0) write(*,*) 'time in joper:',wallout-wallin

end subroutine joper
!!***

