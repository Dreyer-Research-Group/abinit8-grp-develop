!!{\src2tex{textfont=tt}}
!!****f* ABINIT/velfrc
!! NAME
!! velfrc
!!
!! FUNCTION
!! Applies the local and nonlocal expansion up to 2nd order  of the 
!! current density operator to a wavefunction
!!
!! Written by Cyrus Dreyer, Stony Brook and Flatiron CCQ, 2018
!!
!!
!! INPUTS
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

subroutine velfrc(cg,cg1_active,cg1_active_prev,cplex,doccde,docckqde,dtfil,dtset,eigen0,eigen1,eigenq,gh1c_set, &
     & gmet,gprimd,hdr,hdr0,idir_b,ipert_in,istwfk,kg,kg1,mcg,mcg1,mk1mem,mpi_enreg,mpw,mpw1,nband,nfftf,nkpt,npwarr,npwar1, &
     & occ,occkq,paw_ij,pawfgr,pawtab,pertcase,ph1d,psps,prt_eigen1_dk,rmet,rprimd,usecprj,useylmgr1,vtrial,wtk,xred)

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use m_pawcprj, only : pawcprj_type
  use m_xmpi
  use m_cgtools
  use m_occ,      only : occeig,newocc
  use m_paw_ij,     only : paw_ij_type
  use m_pawtab,   only : pawtab_type
  use m_hamiltonian
  use m_pawcprj,     only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_set_zero, pawcprj_axpby
  use m_pawfgr,     only : pawfgr_type
  use m_ioarr
  use m_hdr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'velfrc'
 use interfaces_32_util
 use interfaces_56_recipspace
 use interfaces_66_wfs
 use interfaces_72_response, except_this_one => velfrc
!End of the abilint section

  implicit none

  !Arguments ------------------------------------
  type(pseudopotential_type), intent(in) :: psps
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in), target :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg
  type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  type(pawfgr_type),intent(in) :: pawfgr
  type(hdr_type),intent(inout) :: hdr,hdr0
  integer,intent(in) :: mcg,mcg1,mk1mem,mpw,mpw1,nband(nkpt*dtset%nsppol),nkpt,nfftf
  integer,intent(in) :: istwfk(nkpt),npwarr(nkpt),npwar1(nkpt)
  integer,intent(in) :: usecprj,cplex,ipert_in,idir_b,pertcase,useylmgr1,prt_eigen1_dk
  integer,intent(in) :: kg(3,mpw*nkpt),kg1(3,mpw1*nkpt) ! To check: These were allocated wrong!!!!
  real(dp),intent(in) :: wtk(nkpt)
  real(dp),intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol),occkq(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) ::doccde(dtset%mband*nkpt*dtset%nsppol),docckqde(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: cg(2,mcg)
  real(dp),intent(in) :: cg1_active(2,mcg1)
  real(dp),intent(in) :: cg1_active_prev(2,mcg1)
  real(dp),intent(in) :: xred(3,dtset%natom),rprimd(3,3),rmet(3,3),gmet(3,3),gprimd(3,3)
  real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: eigen0(dtset%mband*nkpt*dtset%nsppol),eigenq(dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: eigen1(2*dtset%mband*dtset%mband*nkpt*dtset%nsppol)
  real(dp),intent(in) :: gh1c_set(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol)
  real(dp), intent(in) :: vtrial(nfftf,dtset%nspden)

  !Local variables-------------------------------
  type(gs_hamiltonian_type) :: gs_hamkq
  type(rf_hamiltonian_type) :: rf_hamkq_a,rf_hamkq_b
  integer :: calcnl_a,calcnl_b
  integer :: npw,npw1,nband_k
  integer :: index_k,index_k1,index_bnd,index_bnd2,index_bnd2_a,index_eij,index_eji,index_eij_a,index_eji_a
  integer :: index_kg,index_kg1
  integer :: ipert_a,ipert_b,ipert_j_a,ipert_j_b,rfdir_a(3),rfdir_b(3)
  integer :: ii 
  integer :: ierr,isppol,ikpt,iband,jband,idab,idir_a,ispden,finq!,idir_b
  integer :: bnd_start,bnd_fin,bnd_start1,bnd_fin1,bnd_start_j,bnd_fin_j,nbd_max,isp_max,nbd,isp
  integer,allocatable :: gbd(:,:),gbd1(:,:),kg_k(:,:),kg1_k(:,:)
  real(dp) :: dab_fac,dedw_fac,fermi_fac,fermi_fac_dab,fij,deleig,deleig_inv
  real(dp) :: dotr1,doti1,dotr2,doti2
  real(dp) :: eta
  real(dp) :: doccde2(dtset%mband*nkpt*dtset%nsppol),occ2(dtset%mband*nkpt*dtset%nsppol), eigen02(dtset%mband*nkpt*dtset%nsppol)
  real(dp) :: qpc(3),rcart_a(3),rcart_b(3),rmod_a,rmod_b,kcart_a(3),kcart_b(3),kmod_a,kmod_b, kpt(3)
  real(dp), allocatable :: denpot(:,:,:),dumr(:,:,:,:),eigen1_a(:)
  real(dp),allocatable :: rocceig(:,:)
  real(dp),allocatable :: cg1_bnd_dcov_a(:,:),cg1_bnd_dcov_b(:,:),cg_bnd(:,:),cg_bnd_j(:,:)
  real(dp),allocatable :: cg1_bndj_dcov_a(:,:),cg1_bndj_dcov_b(:,:)
  complex(dp) :: dedw1(2),dedw2(2),dab1(2),dab1_alt(2),dab2(2),dot_tst,dedw2_alt(2)
  character(len=fnlen) :: fiwf1o_tild

  ! For getgh1c
  type(pawcprj_type) :: cwaveprj(dtset%natom,dtset%nspinor)
  type(pawrhoij_type) :: pawrhoij(dtset%natom*psps%usepaw)
  real(dp),allocatable :: gberry(:,:),gs1c(:,:),gvnl(:,:),gsc(:,:)
  real(dp),allocatable :: gh1c_bnd_a(:,:),gh1c_bnd_b(:,:),ghc_bnd_dcov_b(:,:),gh1c_bnd_tst(:,:)
  real(dp),allocatable :: gh1c_bnd_a_dcov_b(:,:) !TEST
  real(dp),allocatable :: vlocal(:,:,:,:),vlocal1_a(:,:,:,:),vlocal1_b(:,:,:,:)
  real(dp) :: kpq(3),acell(3)
  real(dp) :: vtrial1_a(nfftf,dtset%nspden),vtrial1_b(nfftf,dtset%nspden)
  real(dp) :: eshift,etotal
  real(dp),allocatable :: ylm_k(:,:)
  real(dp),allocatable :: ylmgr1_k(:,:,:)
  real(dp),allocatable :: ylmgr_k(:,:,:)
  real(dp),allocatable :: ylm1_k(:,:)
  real(dp),allocatable :: dkinpw(:) 
  integer :: optlocal,optnl,opt_gvnl1,sij_opt,tim_getgh1c,usevnl,cpopt,optder
  integer :: fform,n4,n5,n6,nkpg,nkpg1
  integer :: simp_fij
  character(len=fnlen) :: fi1o

  real(dp), allocatable :: kpg_k(:,:),kpg1_k(:,:)
  real(dp), allocatable :: kinpw1(:),ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:)

  !TEST: Loop through occupancies
  integer :: iocc,nocc
  real(dp) :: nelect_occ,nelect
  real(dp) :: entropy,fermie_tune
  real(dp) :: occ_tune(dtset%mband*nkpt*dtset%nsppol), doccde_tune(dtset%mband*nkpt*dtset%nsppol)
  real(dp) :: spinmagntarget=-99.99_dp,stmbias=0.0

  ! Files to output info for Max's script
  integer :: max
  integer :: file_dcova_dcovb=30,file_dcova_Hb_ui=31,file_ui_Ha_dcovb=32
  integer :: file_uj_Hb_ui=33,file_ui_Ha_uj=34
  character(22) :: file_name
  real(dp) :: eigen_dcovab(2*dtset%mband*dtset%mband*nkpt*dtset%nsppol)
  real(dp) :: dcov_ab(2,dtset%nsppol,nkpt,dtset%mband,dtset%mband)
  integer :: bantot0,bantot1,dig
  real(dp) :: phasecg(2,dtset%mband*dtset%mband*nkpt*dtset%nsppol)
  character(len=fnlen) :: gkkfilnam


  !*************************************************************************************
  ! For reference, the useri's that I use:
  ! useria -> Triggers velfrc, 1 for first pert, 2 for second
  ! userib -> Triggers joper in cgwf3
  ! useric -> Turn on Drude weight, first digit is direction a, second is direction b
  ! userid -> Default (0) means joper uses local+nl; -1 mean joper uses only local
  ! userie -> 1 toggles simp_fij mod, 2 toggles output for max, -1 for prt_eigen1_dk
  ! userie > 2 toggle occupation tuning, number of steps
  !        userrb -> step size for nelect 
  !      
  !*************************************************************************************


  ! Toggle simple fij
  simp_fij=0
  if (dtset%userie==1)  simp_fij=1

  ! Toggle Max's outputs
  max=0
  if (dtset%userie==2) max=1

  ! idab==1: compute Drude weight; idab==0:no Drude weight
  if (dtset%useric > 0) then
     idab=1
  end if

  ! Set ipert_b and calcnl if we are doing joper
  if (dtset%userib == 2) then
     ipert_b = dtset%natom+99
     ipert_j_b = dtset%natom+1

     ! if userid = -1, we only do local part
     calcnl_b=dtset%userid+1

  else 
     ipert_b=ipert_in
  end if

  ! finq == 1: finite q, finq == 0: q=0
  ! Note: We assume the second perturbation is finite q!!! Also, we will assume TRS for now.
  finq=0
  if (dot_product(dtset%qptn,dtset%qptn) > 1.0d-10) then
     finq=1
     if (idab==1) then
        write(*,*) "WARNING: D_ab for finite q, be careful..."
        !stop
     end if

  end if

  n4=dtset%ngfft(4)
  n5=dtset%ngfft(5)
  n6=dtset%ngfft(6)

  acell(:)=dtset%acell_orig(:,1)

  ! Read in previous vtrial1's and send to all processes
  ! This is probably not done correctly with MPI....
  ! Should be able to use wfk_read_h1mat for this
  if (mpi_enreg%me_kpt==0) then
     open(19,file="vtrial1_1.dat",status="old")
     open(20,file="vtrial1_2.dat",status="old")
     read(19,*) idir_a
     read(20,*) !idir_b don't need this anymore

     do ii=1,nfftf
        read(19,*) vtrial1_a(ii,:)
        read(20,*) vtrial1_b(ii,:)
     end do
     close(19)
     close(20)

  end if

  call xmpi_barrier(mpi_enreg%comm_kpt)
  call xmpi_bcast(vtrial1_a,0,mpi_enreg%comm_cell,ierr)
  call xmpi_bcast(vtrial1_b,0,mpi_enreg%comm_cell,ierr)
  call xmpi_bcast(idir_a,0,mpi_enreg%comm_cell,ierr)
  if (ierr==1) write(*,*) 'ERROR: distribution of eigen1_a'

  ! Read in previous eigen1's and send to all processes.
  ! This is probably not done correctly with MPI....AND GKK!!!!!
  ! Should be able to use wfk_read_h1mat for this
  open(19,file="eigen1_dk.dat",status="old")
  read(19,*) ipert_a,rfdir_a(:)

  if ( prt_eigen1_dk == 1) then
     ABI_ALLOCATE(eigen1_a,(2*dtset%mband*dtset%mband*nkpt*dtset%nsppol))
     if (mpi_enreg%me_kpt==0) then
        do ii=1,2*dtset%mband*dtset%mband*nkpt*dtset%nsppol
           read(19,*) eigen1_a(ii)
        end do

     end if

  end if

  call xmpi_barrier(mpi_enreg%comm_kpt)

  if (prt_eigen1_dk == 1) call xmpi_bcast(eigen1_a,0,mpi_enreg%comm_cell,ierr)
  call xmpi_bcast(ipert_a,0,mpi_enreg%comm_cell,ierr)
  call xmpi_bcast(rfdir_a,0,mpi_enreg%comm_cell,ierr)
  if (ierr==1) write(*,*) 'ERROR: distribution of eigen1_a'

  close(19)


  ! Now that we have ipert_a and ipert_b, open files for optput for max's script
  if (max==1) then

     ! Zero this, will have to assemble over processors when we are done
     eigen_dcovab=zero
     dcov_ab=zero

  end if

  ! vtrial1 should be zero for d/dk, but is not if we use joper. Set by hand:
  if (ipert_a==dtset%natom+99) then
     vtrial1_a=zero
     ipert_j_a = dtset%natom+1
     calcnl_a = 1
  else if (ipert_a==dtset%natom+98) then ! joper, local only
     vtrial1_a=zero
     ipert_j_a = dtset%natom+1
     ipert_a = dtset%natom+99
     calcnl_a = 0
  end if

  ! Initialize hamiltonians for getgh1c and getghc 
  !if (idab==1) then

  ABI_ALLOCATE(gberry,(0,0))
  eshift=zero
  ABI_ALLOCATE(gs1c,(0,0))
  sij_opt=0
  tim_getgh1c=1
  usevnl=1; optlocal=1; optnl=2

  call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%natom,&
       & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
       & paw_ij=paw_ij,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
       & usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

  ! For both a and b perturbations, initialize rf to be dk if using joper
  if (ipert_a == dtset%natom+99) then
     call init_rf_hamiltonian(cplex,gs_hamkq,ipert_j_a,rf_hamkq_a)
  else
     call init_rf_hamiltonian(cplex,gs_hamkq,ipert_a,rf_hamkq_a)
  end if
  if (ipert_b == dtset%natom+99) then
     call init_rf_hamiltonian(cplex,gs_hamkq,ipert_j_b,rf_hamkq_b)
  else
     call init_rf_hamiltonian(cplex,gs_hamkq,ipert_b,rf_hamkq_b)
  end if

  ABI_ALLOCATE(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
  ABI_ALLOCATE(vlocal1_a,(cplex*n4,n5,n6,gs_hamkq%nvloc))
  ABI_ALLOCATE(vlocal1_b,(cplex*n4,n5,n6,gs_hamkq%nvloc))

  ! Turn idir_b into rfdir_b
  rfdir_b(:)=0  
  rfdir_b(idir_b)=1

  ! For nonorthogonal cells, we need to determine the cart magnitudes
  kcart_a=two_pi*matmul(gprimd(:,:),rfdir_a(:))
  kcart_b=two_pi*matmul(gprimd(:,:),rfdir_b(:))
  kmod_a=two_pi*sqrt(dot_product(matmul(gprimd(:,:),rfdir_a(:)),matmul(rfdir_a(:),gprimd(:,:))))
  kmod_b=two_pi*sqrt(dot_product(matmul(gprimd(:,:),rfdir_b(:)),matmul(rfdir_b(:),gprimd(:,:))))


  rcart_a=matmul(rprimd(:,:),rfdir_a(:))
  rcart_b=matmul(rprimd(:,:),rfdir_b(:))
  rmod_a=sqrt(abs(dot_product(matmul(rprimd(:,:),rfdir_a(:)),matmul(rfdir_a(:),rprimd(:,:)))))
  rmod_b=sqrt(abs(dot_product(matmul(rprimd(:,:),rfdir_b(:)),matmul(rfdir_b(:),rprimd(:,:)))))

  ! Also get cartesian q for joper
  qpc(:)=two_pi*matmul(gprimd(:,:),dtset%qptn(:))

  ! Setup factors depending on the perturbations we are calculating
  ! d/dtau contributes 2/a, d/dk contributes a/2\pi.
  ! need to test for nonorthogonal cells, also only one direction at a time

  if (dot_product(rfdir_a,rfdir_a)>1) then
     write(*,*) "ERROR: ONLY ONE DIRECTION AT A TIME FOR FIRST PERT!!"
     stop
  end if


  ! Take care of all of the factors
  dedw_fac=one
  if (ipert_a <= dtset%natom) then

     !dedw_fac=dedw_fac*two/(dot_product(rfdir_a,acell))
     dedw_fac=dedw_fac*two/rmod_a  

  else if (ipert_a==dtset%natom+1) then! d/dk

     !dedw_fac=dedw_fac*(dot_product(rfdir_a,acell))/(two_pi)
     dedw_fac=dedw_fac/kmod_a

     write(*,*) "WARNING, KMOD NOT TESTED!!!"

     ! to make consistant with d/dk:
  else if (ipert_a==dtset%natom+99) then! finite q d/dk  
     dedw_fac=dedw_fac/4.
  end if

  if (ipert_b <= dtset%natom.and.dtset%userib /= 2) then

     !dedw_fac=dedw_fac*two/(dot_product(rfdir_b,acell))
     dedw_fac=dedw_fac*two/rmod_b

  else if (ipert_b==dtset%natom+1) then! d/dk

     !dedw_fac=dedw_fac*(dot_product(rfdir_b,acell))/(two_pi)
     dedw_fac=dedw_fac/kmod_b

     write(*,*) "WARNING, KMOD NOT TESTED!!!"

     !TEST: DO I NEED THIS
  else if (ipert_b==dtset%natom+99) then! finite q d/dk  
     dedw_fac=dedw_fac/4.

  end if

  ! For cycle
  nbd_max=size(mpi_enreg%proc_distrb,2)
  isp_max=size(mpi_enreg%proc_distrb,3)


  ! For occ loop
  if (dtset%userie > 2) then
     nocc=dtset%userie
  else
     nocc=0
  end if

  ! Loop over occupancies
  do iocc=0,nocc


     ! Initialize indicies
     index_k=1; index_k1=1
     index_bnd=1; index_bnd2=1
     !index_bnd2_a=1 ! For a pert eigen1_dk
     dedw1=zero; dedw2=zero; dedw2_alt=zero
     dab1=zero; dab1_alt=zero; dab2=zero
     index_kg=0; index_kg1=0

     if (iocc>0) then
        nelect_occ=dtset%nelect+iocc*dtset%userrb
        call newocc(doccde_tune,eigen0,entropy,fermie_tune,spinmagntarget,dtset%mband,nband,&
             &  nelect_occ,nkpt,dtset%nspinor,dtset%nsppol,occ_tune,dtset%occopt,dtset%prtvol,&
             &  stmbias,dtset%tphysel,dtset%tsmear,wtk)

     else
        doccde_tune=doccde
        occ_tune=occ
     end if

     ! spin and kpoint loops
     do isppol=1,dtset%nsppol

        !  Continue to initialize the Hamiltonian, if not joper

        ! For a perturbation
        call rf_transgrid_and_pack(isppol,dtset%nspden,psps%usepaw,cplex,nfftf,dtset%nfft,&
             & dtset%ngfft,gs_hamkq%nvloc,pawfgr,mpi_enreg,vtrial,vtrial1_a,vlocal,vlocal1_a)        
        call load_spin_hamiltonian(gs_hamkq,isppol,vlocal=vlocal,with_nonlocal=.true.)
        call load_spin_rf_hamiltonian(rf_hamkq_a,isppol,vlocal1=vlocal1_a,with_nonlocal=.true.)

        ! For b perturbation     
        call rf_transgrid_and_pack(isppol,dtset%nspden,psps%usepaw,cplex,nfftf,dtset%nfft,&
             & dtset%ngfft,gs_hamkq%nvloc,pawfgr,mpi_enreg,vtrial,vtrial1_b,vlocal,vlocal1_b)
        call load_spin_rf_hamiltonian(rf_hamkq_b,isppol,vlocal1=vlocal1_b,with_nonlocal=.true.)

        do ikpt=1,nkpt

           kpt(:)=dtset%kptns(:,ikpt)

           ! get number of bands for this k and isspol
           nband_k=nband(ikpt+(isppol-1)*nkpt)

           ! I am now not sure whether we should be inncrementing by
           ! nband or mband. For now, stop if they are not the same
           if (nband_k /= dtset%mband) then
              write(*,*) "ERROR:nband_k does not equal mband. Stop to be safe!"
              stop
           end if

           ! Test if this k and band belong to me
           if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband(ikpt),isppol,mpi_enreg%me_kpt)) then
              index_bnd=index_bnd+nband(ikpt)
              index_bnd2=index_bnd2+2.*(nband(ikpt))**2
              cycle
           end if

           ! Allocate things that depend on npw, 
           npw=npwarr(ikpt)
           npw1=npwar1(ikpt)
           ABI_ALLOCATE(kg1_k,(3,npw1))
           ABI_ALLOCATE(kg_k,(3,npw))

           ! Assume that either BOTH are finite q or NEITHER
           ! is there a reason to generalize?? Then need to FFT
           ABI_ALLOCATE(cg1_bnd_dcov_a,(2,npw1*dtset%nspinor))
           ABI_ALLOCATE(cg1_bnd_dcov_b,(2,npw1*dtset%nspinor))

           ! Max wants off diagonal elements of covarient derivatives
           if (max==1) then
              ABI_ALLOCATE(cg1_bndj_dcov_a,(2,npw1*dtset%nspinor))
              ABI_ALLOCATE(cg1_bndj_dcov_b,(2,npw1*dtset%nspinor))
           end if

           ! Get kg's for this kpoint
           kg_k(:,:)=kg(:,1+index_kg:npw)
           kg1_k(:,:)=kg1(:,1+index_kg1:npw1)

           index_kg=index_kg+npw
           index_kg1=index_kg1+npw1

           ABI_ALLOCATE(cg_bnd,(2,npw*dtset%nspinor))
           ABI_ALLOCATE(gh1c_bnd_a,(2,npw1*dtset%nspinor))
           ABI_ALLOCATE(gh1c_bnd_b,(2,npw1*dtset%nspinor))
           ABI_ALLOCATE(ghc_bnd_dcov_b,(2,npw1*dtset%nspinor))          
           ABI_ALLOCATE(gh1c_bnd_tst,(2,npw1*dtset%nspinor))  !TEST        
           ABI_ALLOCATE(gh1c_bnd_a_dcov_b,(2,npw1*dtset%nspinor)) !TEST
           ABI_ALLOCATE(gsc,(2,npw*dtset%nspinor))
           ABI_ALLOCATE(gvnl,(2,npw*dtset%nspinor))

           !TEST: Alt Kubo
           ABI_ALLOCATE(cg_bnd_j,(2,npw*dtset%nspinor))


           ! YLM

           ABI_ALLOCATE(ylm_k,(npw,psps%mpsang*psps%mpsang*psps%useylm))
           ABI_ALLOCATE(ylmgr1_k,(npw1,3+6*((ipert_in-dtset%natom)/10),psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
           ABI_ALLOCATE(ylmgr_k,(npw,3+6*((ipert_in-dtset%natom)/10),psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
           ABI_ALLOCATE(ylm1_k,(npw1,psps%mpsang*psps%mpsang*psps%useylm))

           if (psps%useylm==1) then
              ! First, get spherical harmonics at this k point, if needed, this seems to give seg faults...
              optder=1
              call initylmg(gprimd,kg1_k,dtset%kptns(:,ikpt),1,mpi_enreg,psps%mpsang,npw1,&
                   &  nband(ikpt),1,npwar1(ikpt),dtset%nsppol,optder,rprimd,ylm1_k,ylmgr1_k)
              call initylmg(gprimd,kg_k,dtset%kptns(:,ikpt),1,mpi_enreg,psps%mpsang,npw,&
                   &  nband(ikpt),1,npwarr(ikpt),dtset%nsppol,optder,rprimd,ylm_k,ylmgr_k)
           else

              ylm_k=zero;ylmgr_k=zero;ylm1_k=zero;ylmgr1_k=zero

           end if

           ! Setup k-dependent hamiltonian stuff for getghc and getgh1c.
           kpq(:)=dtset%qptn(:)+dtset%kptns(:,ikpt)

           ! For a perturbation
           ! Take care of ipert for joper imp
           if (ipert_a==dtset%natom+99) then
              ipert_j_a = dtset%natom+1
           else
              ipert_j_a=ipert_a
           end if
           call getgh1c_setup(gs_hamkq,rf_hamkq_a,dtset,psps,dtset%kptns(:,ikpt),kpq,idir_a,ipert_j_a,& ! In
                & dtset%natom,rmet,gprimd,gmet,istwfk(ikpt),npw,npw1,&                          ! In
                & useylmgr1,kg_k,ylm_k,kg1_k,ylm1_k,ylmgr1_k,&                           ! In
                & dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)         ! Out

           ! For b perturbation
           ! This is set to npw to be consistent with getgh1c_setup...Why is it not npw1???
           call load_k_rf_hamiltonian(rf_hamkq_b,npw_k=npw,dkinpw_k=dkinpw)

           do iband=1,nband(ikpt)

              ! Check if this is a "buffer band"
              if (nband(ikpt)-iband < dtset%nbdbuf .and. occ(index_bnd+iband-1) > 1.d-4) then              
                 write(*,*) "WARNING: buffer band",iband,"occupied for kpt",ikpt
              end if

              ! First term in Eq. 37 of Max's Oct 4 notes
              ! _________________________________________

              ! Get one band out of cg1's
              bnd_start=index_k+npw*dtset%nspinor*(iband-1)
              bnd_fin=index_k+npw*dtset%nspinor*(iband)-1
              bnd_start1=index_k1+npw1*dtset%nspinor*(iband-1)
              bnd_fin1=index_k1+npw1*dtset%nspinor*(iband)-1

              ! Again, either NEITHER or BOTH are finite q
              cg1_bnd_dcov_a(:,:)=cg1_active_prev(:,bnd_start1:bnd_fin1)
              cg1_bnd_dcov_b(:,:)=cg1_active(:,bnd_start1:bnd_fin1)              


              ! Max's 11/19/18 equation, Berry term, Drude Weight
              ! This is much more efficient WRT kpt convergence
              ! _________________________________________________________
              if (idab==1) then

                 ! D_ab factor: each k derivative comes with a/2\pi, need to test for non-orthorhombic
                 !dab_fac=rprimd(idir_a,idir_a)*rprimd(idir_b,idir_b)/(4.*pi**2)
                 dab_fac=dedw_fac


                 ! STEP 1: Get one band out of cg, gh1c_set for testing
                 ! right now only q=0
                 cg_bnd(:,:)=cg(:,bnd_start:bnd_fin)
                 gh1c_bnd_tst(:,:)=gh1c_set(:,bnd_start:bnd_fin) ! This is for testing

                 ! STEP 2.1: H^(0)|\dcover_b u_i>:
                 cpopt=-1
                 call getghc(cpopt,cg1_bnd_dcov_b,cwaveprj,ghc_bnd_dcov_b,gsc,gs_hamkq,gvnl,eshift,mpi_enreg,&
                      &   1,dtset%prtvol,sij_opt,tim_getgh1c,0)


                 ! Subtract eigenvalues:
                 do ii=1,npw1*dtset%nspinor
                    ghc_bnd_dcov_b(:,ii)=ghc_bnd_dcov_b(:,ii)-eigen0(index_bnd+iband-1)*cg1_bnd_dcov_b(:,ii)
                 end do

                 ! STEP 3.1: make < dcover_a u_i | H-e_i | dcover_b u_i>
                 call dotprod_g(dotr1,doti1,istwfk(ikpt),npw*dtset%nspinor,2,cg1_bnd_dcov_a,ghc_bnd_dcov_b, &
                      & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

                 ! STEP 4.1: calc D_ab^berry
                 dab1_alt(isppol)=dab1_alt(isppol)-two*dab_fac*wtk(ikpt)* &
                      & occ_tune(index_bnd+iband-1)*cmplx(dotr1,doti1)

                 ! Dab using fermi-surface integral:
                 ! \overline{f}_{n\nk}\langle\unk\vert\hat{H}^{k_\alpha}_{\textbf{k}}\vert\unk\rangle
                 ! \times\langle\unk\hat{H}^{k_\beta}_{\textbf{k}}\vert\unk\rangle

                 if (prt_eigen1_dk == 1) then 

                    index_eij=2*(iband+(iband-1)*nband(ikpt))+index_bnd2-1

                    dab1(isppol)=dab1(isppol)+wtk(ikpt)*dab_fac* &
                         & doccde_tune(index_bnd+iband-1)* & 
                         & cmplx(eigen1_a(index_eij-1),-eigen1_a(index_eij))*cmplx(eigen1(index_eij-1),-eigen1(index_eij))

                 else

                    ! STEP 2.2: H^a|u_i> and H^b|u_i>. This is maybe not the best way since i need both H_a and H_b
                    ! H^a|u_i>
                    if (ipert_a==dtset%natom+99) then !Joper for finite q dk
                       ! TEST: with spinors
                       gh1c_bnd_a=cg_bnd
                       call joper(calcnl_a,gh1c_bnd_a,dtset,gprimd,idir_a,kg_k,kpt,mpi_enreg,npw,psps,qpc,rprimd)
                    else
                       call getgh1c(dtset%berryopt,cg_bnd,cwaveprj,gh1c_bnd_a,gberry,gsc,gs_hamkq,gvnl,idir_a,ipert_a,eshift,&
                            &   mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamkq_a,sij_opt,tim_getgh1c,usevnl)
                    end if


                    ! H^b|u_i>:
                    if (ipert_b==dtset%natom+99) then !Joper for finite q dk
                       gh1c_bnd_b=cg_bnd
                       call joper(calcnl_b,gh1c_bnd_b,dtset,gprimd,idir_b,kg_k,kpt,mpi_enreg,npw,psps,qpc,rprimd)
                    else
                       call getgh1c(dtset%berryopt,cg_bnd,cwaveprj,gh1c_bnd_b,gberry,gsc,gs_hamkq,gvnl,idir_b,ipert_b,eshift,&
                            &   mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamkq_b,sij_opt,tim_getgh1c,usevnl)
                    end if

                    ! STEP 3.2: < u_i | H_b | u_i > < u_i | H_a | u_i >
                    !gh1c_bnd_a(2,:)=-gh1c_bnd_a(2,:)               ! take CC's
                    call dotprod_g(dotr1,doti1,istwfk(ikpt),npw*dtset%nspinor,2,cg_bnd,gh1c_bnd_b, &
                         & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
                    call dotprod_g(dotr2,doti2,istwfk(ikpt),npw*dtset%nspinor,2,gh1c_bnd_a,cg_bnd, &
                         & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

                    dab1(isppol)=dab1(isppol)+wtk(ikpt)*dab_fac* &
                         & doccde_tune(index_bnd+iband-1)* & 
                         & cmplx(dotr1,doti1)*cmplx(dotr2,doti2)


                 end if ! prt_eigen1_dk

              end if !idab==1


              ! _____________________________________________

              ! Compute overlaps and sum
              call dotprod_g(dotr1,doti1,istwfk(ikpt),npw1*dtset%nspinor,2,cg1_bnd_dcov_a,cg1_bnd_dcov_b, &
                   & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
              call dotprod_g(dotr2,doti2,istwfk(ikpt),npw1*dtset%nspinor,2,cg1_bnd_dcov_b,cg1_bnd_dcov_a, &
                   & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

              dedw1(isppol)=dedw1(isppol)+wtk(ikpt)*half*dedw_fac* &
                   & occ_tune(index_bnd+iband-1)*cmplx(dotr1-dotr2,doti1-doti2)

              ! TEST: Check for NaN
              if (dotr1 /= dotr1) then
                 write(*,*) "dotr1 NaN at kpt",ikpt,"band",iband 
                 stop
              else if (doti1 /= doti1) then
                 write(*,*) "doti1 NaN at kpt",ikpt,"band",iband
                 stop
              else if (dotr2 /= dotr2) then
                 write(*,*) "dotr2 NaN at kpt",ikpt,"band",iband
                 stop
              else if (doti2 /= doti2) then
                 write(*,*) "doti2 NaN at kpt",ikpt,"band",iband
                 stop
              end if



              ! Second term (Kubo) in Eq. 37 of Max's Oct 4 notes
              ! _________________________________________________  
              do jband=1,nband(ikpt)

                 ! User defined eta:
                 eta=dtset%userra
                 if (abs(eta)<1.0e-20) eta=1.0d-5

                 ! Eq. 6 in Max's paper PRB 62 15283
                 ! _____________________________________________

                 ! Take care of finite q
                 if (finq==1) then
                    occ2=occkq
                    doccde2=docckqde
                    eigen02=eigenq
                 else
                    occ2=occ_tune
                    doccde2=doccde_tune
                    eigen02=eigen0
                 end if

                 deleig=(eigen02(index_bnd+jband-1)-eigen0(index_bnd+iband-1))

                 ! Old style fij
                 if (simp_fij==0) then
                    if (abs(deleig) < 1.d-6) then
                       fij=-half*(doccde_tune(index_bnd+jband-1)+doccde2(index_bnd+iband-1))

                       ! Step 2: Add small imaginary part to remaining energy denominator
                       deleig_inv=real(one/cmplx(real(eigen02(index_bnd+jband-1)-eigen0(index_bnd+iband-1)),real(eta)))

                    else
                       fij=(occ2(index_bnd+jband-1)-occ_tune(index_bnd+iband-1))/deleig
                       deleig_inv=one/deleig
                    end if

                    fermi_fac=wtk(ikpt)*half*fij*deleig_inv 

                 else ! Simplified fij

                    fermi_fac=wtk(ikpt)*half*(occ2(index_bnd+jband-1)-occ_tune(index_bnd+iband-1))/(deleig**2+eta**2)

                 end if

                 ! <i|daH|j><j|dbH|i> from eigenvalues
                 index_eij=2*(jband+(iband-1)*nband(ikpt))+index_bnd2-1
                 index_eji=2*(iband+(jband-1)*nband(ikpt))+index_bnd2-1
                 
                 !index_eij=(2*iband-1+(jband-1)*2*dtset%mband+band2tot_index
                 
                 ! TEST: For few kpoints, we can just use the eigenvalues 
                 if (prt_eigen1_dk == 1) then
                    dedw2(isppol)=dedw2(isppol)+fermi_fac*dedw_fac* &
                         & cmplx(eigen1_a(index_eij-1),-eigen1_a(index_eij))*cmplx(eigen1(index_eji-1),-eigen1(index_eji))

                 else
                    !************************************************************************
                    ! For many kpts, eigen1 is corrupted, so lets recalculate
                    ! In fact, for the Dab, we computed H^a|u_i> above and stored it in
                    ! gh1c_bnd_a (and gh1c_bnd_b).

                    cg_bnd(:,:)=cg(:,bnd_start:bnd_fin)
                    !gh1c_bnd_tst(:,:)=gh1c_set(:,bnd_start:bnd_fin) ! This is for testing

                    ! Step 1: get u_j
                    bnd_start_j=index_k+npw*dtset%nspinor*(jband-1)
                    bnd_fin_j=index_k+npw*dtset%nspinor*(jband)-1
                    cg_bnd_j(:,:)=cg(:,bnd_start_j:bnd_fin_j)

                    !TEST b perturbation
                    gh1c_bnd_tst(:,:)=gh1c_set(:,bnd_start_j:bnd_fin_j)

                    gh1c_bnd_a=zero; gh1c_bnd_b=zero

                    ! Step 2: Apply first order hamiltonian for pertubation a and b
                    ! H^a|u_i>
                    ! Use joper for "finite q dk"
                    if (ipert_a == dtset%natom+99) then
                       ! Does not work with spinors yet
                       gh1c_bnd_a=cg_bnd
                       call joper(calcnl_a,gh1c_bnd_a,dtset,gprimd,idir_a,kg_k,kpt,mpi_enreg,npw,psps,qpc,rprimd)
                    else
                       call getgh1c(dtset%berryopt,cg_bnd,cwaveprj,gh1c_bnd_a,gberry,gsc,gs_hamkq,gvnl,idir_a,ipert_a,eshift,&
                            &   mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamkq_a,sij_opt,tim_getgh1c,usevnl)

                    end if

                    ! H^b|u_i>
                    if (ipert_b == dtset%natom+99) then
                       gh1c_bnd_b=cg_bnd
                       call joper(calcnl_b,gh1c_bnd_b,dtset,gprimd,idir_b,kg_k,kpt,mpi_enreg,npw,psps,qpc,rprimd)
                    else
                       call getgh1c(dtset%berryopt,cg_bnd,cwaveprj,gh1c_bnd_b,gberry,gsc,gs_hamkq,gvnl,idir_b,ipert_b,eshift,&
                            &   mpi_enreg,optlocal,optnl,opt_gvnl1,rf_hamkq_b,sij_opt,tim_getgh1c,usevnl)       
                    end if

                    ! Step 3: dot with gh1c_bnd
                    ! <u_j|Hb|u_i>
                    call dotprod_g(dotr1,doti1,istwfk(ikpt),npw*dtset%nspinor,2,cg_bnd_j,gh1c_bnd_b, &
                         & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
                    ! <u_i|Ha|u_j>
                    call dotprod_g(dotr2,doti2,istwfk(ikpt),npw*dtset%nspinor,2,gh1c_bnd_a,cg_bnd_j, &
                         & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

                    ! NOTE: Joper seems to converge slower with, e.g., cutoff
                    ! Not sure why I need opposite sign here...
                    dedw2(isppol)=dedw2(isppol)+fermi_fac*dedw_fac* &
                         & cmplx(dotr1,doti1)*cmplx(dotr2,doti2)


                    ! TEST
                    !   write(*,'(3f6.2,2i5,8e20.10e2)') kpt(:),iband,jband,eigen1_a(index_eij-1),eigen1_a(index_eij),dotr2,doti2,eigen1(index_eji-1),eigen1(index_eji),dotr1,doti1
                    
                 end if ! eigenvals versus doing it myself

                 !**************************************************************************

                 ! useric toggles max's new equation 11/18/2018, so we do not include the other energy derivative
                 if (idab==1) then

                    if (simp_fij==0) then
                       ! TEST: Old style energy denominators
                       if (iband==jband) then
                          fij=zero
                       end if

                       fermi_fac_dab=dab_fac*wtk(ikpt)*fij

                    else ! Simplified fermi factor
                       fermi_fac_dab=dab_fac*wtk(ikpt)*(occ2(index_bnd+jband-1)-occ_tune(index_bnd+iband-1))/(deleig+eta)!fij

                    end if

                    if (prt_eigen1_dk == 1) then

                       dab2(isppol)=dab2(isppol)+fermi_fac_dab* &
                            & cmplx(eigen1_a(index_eij-1),-eigen1_a(index_eij))*cmplx(eigen1(index_eji-1),-eigen1(index_eji))

                    else
                       dab2(isppol)=dab2(isppol)+fermi_fac_dab* &
                            & cmplx(dotr1,doti1)*cmplx(dotr2,doti2)

                    end if

                 end if

                 ! Here we will get the whole matrix for Max's output
                 if (max==1) then

                    ! Step 1: Get covariant derivatives for band j
                    cg1_bndj_dcov_a(:,:)=cg1_active_prev(:,bnd_start_j:bnd_fin_j)
                    cg1_bndj_dcov_b(:,:)=cg1_active(:,bnd_start_j:bnd_fin_j)              

                    ! Compute overlaps: <d/dtau_ka u_m|d/dk_b u_n> - <dk_b u_m|d/dtau_ka u_n>
                    call dotprod_g(dotr1,doti1,istwfk(ikpt),npw1*dtset%nspinor,2,cg1_bndj_dcov_a,cg1_bnd_dcov_b, &
                         & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
                    call dotprod_g(dotr2,doti2,istwfk(ikpt),npw1*dtset%nspinor,2,cg1_bndj_dcov_b,cg1_bnd_dcov_a, &
                         & mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

                    !eigen_dcovab(index_eij-1)=dotr1-dotr2
                    !eigen_dcovab(index_eij)=doti1-doti2
                    dcov_ab(1,isppol,ikpt,iband,jband)=dotr1-dotr2
                    dcov_ab(2,isppol,ikpt,iband,jband)=doti1-doti2

                 end if

              end do ! jband                 
           end do ! iband         

           index_k=index_k+npw*dtset%nspinor*nband(ikpt)
           index_k1=index_k1+npw1*dtset%nspinor*nband(ikpt)
           index_bnd=index_bnd+nband(ikpt)
           index_bnd2=index_bnd2+2.*(nband(ikpt))**2

           index_bnd2_a=index_bnd2_a+2.*(nband(ikpt))**2

           ABI_DEALLOCATE(cg1_bnd_dcov_a)
           ABI_DEALLOCATE(cg1_bnd_dcov_b)

           if (max==1) then
              ABI_DEALLOCATE(cg1_bndj_dcov_a)
              ABI_DEALLOCATE(cg1_bndj_dcov_b)
           end if

           !ABI_DEALLOCATE(rocceig)
           ABI_DEALLOCATE(kg1_k)
           ABI_DEALLOCATE(kg_k)

           !if (idab==1) then
           ABI_DEALLOCATE(gsc)
           ABI_DEALLOCATE(gvnl)
           ABI_DEALLOCATE(cg_bnd)
           ABI_DEALLOCATE(cg_bnd_j)
           ABI_DEALLOCATE(gh1c_bnd_a)
           ABI_DEALLOCATE(gh1c_bnd_b)
           ABI_DEALLOCATE(ghc_bnd_dcov_b)
           ABI_DEALLOCATE(gh1c_bnd_tst) !TEST
           ABI_DEALLOCATE(gh1c_bnd_a_dcov_b)
           !ABI_DEALLOCATE(dkinpw)
           if (allocated(ylm_k)) then 
              ABI_DEALLOCATE(ylm_k)
              ABI_DEALLOCATE(ylmgr1_k)
              ABI_DEALLOCATE(ylmgr_k)
              ABI_DEALLOCATE(ylm1_k)
           end if
           !end if


           ! TEST: Lets us know where we are
           !write(*,*) "Done with KPT",ikpt, "of",nkpt

        end do ! ikpt

     end do ! isppol


#ifdef HAVE_MPI

     call xmpi_barrier(mpi_enreg%comm_kpt)
     call xmpi_sum_master(dedw1,0,mpi_enreg%comm_kpt,ierr)
     if (ierr==1) write(*,*) 'ERROR: Summation of dedw1'
     call xmpi_sum_master(dedw2,0,mpi_enreg%comm_kpt,ierr)
     if (ierr==1) write(*,*) 'ERROR: Summation of dedw2'
     call xmpi_sum_master(dab1,0,mpi_enreg%comm_kpt,ierr)
     if (ierr==1) write(*,*) 'ERROR: Summation of dab1'
     call xmpi_sum_master(dab1_alt,0,mpi_enreg%comm_kpt,ierr)
     if (ierr==1) write(*,*) 'ERROR: Summation of dab1_alt'
     call xmpi_sum_master(dab2,0,mpi_enreg%comm_kpt,ierr)
     if (ierr==1) write(*,*) 'ERROR: Summation of dab2'

     if (max==1) then
        !call xmpi_sum_master(eigen_dcovab,0,mpi_enreg%comm_kpt,ierr)
        !if (ierr==1) write(*,*) 'ERROR: Summation of eigen_dcovab'
        call xmpi_sum_master(dcov_ab,0,mpi_enreg%comm_kpt,ierr)
        if (ierr==1) write(*,*) 'ERROR: Summation of dcov_ab'

     end if

     if (mpi_enreg%me_kpt==0) then

        ! Warn about corrupted eigen1
        if (prt_eigen1_dk == 1) write(*,*) "WARNING: Using eigen1 ", &
             & "which may be corrupted for large number of kpoints"
        
        ! write out to *.out file
        do isppol=1,dtset%nsppol

           ! Set nelect for rigid band
           if (iocc==0) then
              nelect=dtset%nelect
           else
              nelect=nelect_occ
           end if

           ! Drude weight
           if (idab==1) then

              ! Leading term in drude weight
              if (idir_a==idir_b) then 

                 if (dtset%nsppol==2) then
                    dab_fac=-1.*nelect/2.
                 else
                    dab_fac=-1.*nelect
                 end if
              else 
                 dab_fac=0
              end if

              ! spin, re/img berry (alt way),re/img Kubo, re/img total
              write(*,'(a10,i5,f12.4,6e20.10e2)') 'D_ab',isppol,nelect,dab1_alt(isppol),dab2(isppol),dab_fac-(dab1_alt(isppol)+dab2(isppol))
              write(*,'(a20,i5,2e20.10e2)') 'Tot tst ',isppol,dab1(isppol)

           else
              
              write(*,'(a10,i5,f12.4,6e20.10e2)') 'dE/dw',isppol,nelect,dedw1(isppol),dedw2(isppol),dedw1(isppol)-dedw2(isppol)

           end if

        end do ! isppol

        ! Write out covarient FO eigenvalues.
        if (max==1) then
           ii=1
           do isppol=1,dtset%nsppol
              do ikpt=1,nkpt
                 do iband=1,nband(ikpt)
                    do jband=1,nband(ikpt)
                       eigen_dcovab(ii)=dcov_ab(1,isppol,ikpt,iband,jband)
                       eigen_dcovab(ii+1)=dcov_ab(2,isppol,ikpt,iband,jband)

                       ii=ii+2

                       ! Test for hermaticity
                       !write(*,'(a6,3i5,2e12.4)') "HERM",ikpt,iband,jband, &
                       !     & dcov_ab(1,isppol,ikpt,iband,jband)+dcov_ab(1,isppol,ikpt,jband,iband),&
                       !     & dcov_ab(2,isppol,ikpt,jband,iband)-dcov_ab(2,isppol,ikpt,iband,jband)
                    end do
                 end do
              end do
           end do

           ! Use outgkk to output in correct format
           bantot0=sum(nband(1:nkpt*dtset%nsppol))
           bantot1=dtset%mband*dtset%mband*nkpt*dtset%nsppol
           phasecg(1,:) = one ! This is not even really used
           phasecg(2,:) = zero
           dig=((ipert_a*100+idir_a)*1000+ipert_b)*100+idir_b

           write(*,*) "BEFORE OUTGKK",dig

           call appdig(dig,dtfil%fnameabo_gkk,gkkfilnam) 
           call outgkk(bantot0,bantot1,gkkfilnam,eigen0,eigen_dcovab,hdr0,hdr,mpi_enreg,phasecg)



        end if



     end if ! me_kpt==0
#endif

  end do ! iocc

end subroutine velfrc
!!***

