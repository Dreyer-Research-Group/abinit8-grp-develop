#
# Makefile for ABINIT                                      -*- Automake -*-
# Generated by make-makefiles-corelibs on 2021/06/01 15:12:57 +0000

#
# IMPORTANT NOTE
#
# Any manual change to this file will systematically be overwritten.
# Please modify the make-makefiles-corelibs script or its config file instead.
#

AM_CFLAGS = @ABI_CPPFLAGS@
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_67_common@

AM_CPPFLAGS = \
	@src_54_abiutil_incs@ \
	@src_59_ionetcdf_incs@ \
	@src_27_toolbox_oop_incs@ \
	@src_42_libpaw_incs@ \
	@src_65_paw_incs@ \
	@src_62_wvl_wfs_incs@ \
	@src_66_nonlocal_incs@ \
	@src_45_geomoptim_incs@ \
	@src_17_libtetra_ext_incs@ \
	@src_49_gw_toolbox_oop_incs@ \
	@src_62_iowfdenpot_incs@ \
	@src_18_timing_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_56_mixing_incs@ \
	@src_56_io_mpi_incs@ \
	@src_14_hidewrite_incs@ \
	@src_20_datashare_incs@ \
	@src_64_psp_incs@ \
	@src_32_util_incs@ \
	@src_28_numeric_noabirule_incs@ \
	@src_44_abitypes_defs_incs@ \
	@src_51_manage_mpi_incs@ \
	@src_43_wvl_wrappers_incs@ \
	@src_62_poisson_incs@ \
	@src_16_hideleave_incs@ \
	@src_56_recipspace_incs@ \
	@src_61_occeig_incs@ \
	@src_52_fft_mpi_noabirule_incs@ \
	@src_56_xc_incs@ \
	@src_66_wfs_incs@ \
	@src_41_xc_lowlevel_incs@ \
	@src_41_geometry_incs@ \
	@src_01_linalg_ext_incs@ \
	@src_53_spacepar_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_53_ffts_incs@ \
	@src_incs_incs@ \
	@lib_bigdft_incs@ \
	@lib_gpu_incs@ \
	@lib_netcdf_incs@ \
	@lib_mpi_incs@ \
	@lib_wannier90_incs@ \
	@fallbacks_incs@ \
	@abi_extra_incs@ \
	@fc_mod_incs@

# Regular source files
lib67_common_srcs = \
	berryphase.F90 \
	berryphase_new.F90 \
	calc_efg.F90 \
	calc_fc.F90 \
	cgq_builder.F90 \
	cgwf.F90 \
	chern_number.F90 \
	conducti_nc.F90 \
	conducti_paw.F90 \
	conducti_paw_core.F90 \
	constrf.F90 \
	dielmt.F90 \
	dieltcel.F90 \
	emispec_paw.F90 \
	energy.F90 \
	etotfor.F90 \
	evdw_wannier.F90 \
	extraprho.F90 \
	extrapwf.F90 \
	forces.F90 \
	forstr.F90 \
	forstrnps.F90 \
	fresid.F90 \
	fresidrsp.F90 \
	getcgqphase.F90 \
	init_e_field_vars.F90 \
	initberry.F90 \
	initorbmag.F90 \
	initmv.F90 \
	initro.F90 \
	ionion_realspace.F90 \
	ionion_surface.F90 \
	jellium.F90 \
	ks_ddiago.F90 \
	m_epjdos.F90 \
	mag_constr.F90 \
	mag_constr_e.F90 \
	make_efg_el.F90 \
	make_efg_ion.F90 \
	make_grad_berry.F90 \
	mkcore_inner.F90 \
	mkcore_paw.F90 \
	mkcore_wvl.F90 \
	mkgrid_fft.F90 \
	mklocl.F90 \
	mklocl_realspace.F90 \
	mklocl_recipspace.F90 \
	mklocl_wavelets.F90 \
	mkresi.F90 \
	mkrho.F90 \
	m_ioarr.F90 \
	m_iowf.F90 \
	mlwfovlp.F90 \
	mlwfovlp_proj.F90 \
	mlwfovlp_projpaw.F90 \
	mlwfovlp_pw.F90 \
	mlwfovlp_radial.F90 \
	mlwfovlp_seedname.F90 \
	mlwfovlp_setup.F90 \
	mlwfovlp_ylmfac.F90 \
	mlwfovlp_ylmfar.F90 \
	moddiel.F90 \
	msig.F90 \
	multipoles_out.F90 \
	m_eprenorms.F90 \
	m_plowannier.F90 \
	m_vcoul.F90 \
	newkpt.F90 \
	nres2vres.F90 \
	odamix.F90 \
	optics_vloc.F90 \
	partial_dos_fractions.F90 \
	partial_dos_fractions_paw.F90 \
	posdoppler.F90 \
	poslifetime.F90 \
	posratecore.F90 \
	prteigrs.F90 \
	prtene.F90 \
	prtefield.F90 \
	prtimg.F90 \
	prtrhomxmn.F90 \
	prtxf.F90 \
	prtxvf.F90 \
	rhotov.F90 \
	scprqt.F90 \
	setup_positron.F90 \
	setup1.F90 \
	setvtr.F90 \
	spatialchempot.F90 \
	spin_current.F90 \
	stress.F90 \
	strhar.F90 \
	sygrad.F90 \
	symrhg.F90 \
	uderiv.F90 \
	update_e_field_vars.F90 \
	vdw_dftd2.F90 \
	vdw_dftd3.F90 \
	vso_realspace_local.F90 \
	vtorhotf.F90 \
	waveformat.F90 \
	wf_mixing.F90 \
	wvl_initro.F90 \
	wvl_mkrho.F90 \
	interfaces_67_common.F90

# Source files depending on conditionals


# Library description
noinst_LIBRARIES = lib67_common.a

lib67_common_a_SOURCES= $(lib67_common_srcs)

CLEANFILES = \
	berryphase_cpp.f90 \
	berryphase_new_cpp.f90 \
	calc_efg_cpp.f90 \
	calc_fc_cpp.f90 \
	cgq_builder_cpp.f90 \
	cgwf_cpp.f90 \
	chern_number_cpp.f90 \
	conducti_nc_cpp.f90 \
	conducti_paw_cpp.f90 \
	conducti_paw_core_cpp.f90 \
	constrf_cpp.f90 \
	dielmt_cpp.f90 \
	dieltcel_cpp.f90 \
	emispec_paw_cpp.f90 \
	energy_cpp.f90 \
	etotfor_cpp.f90 \
	evdw_wannier_cpp.f90 \
	extraprho_cpp.f90 \
	extrapwf_cpp.f90 \
	forces_cpp.f90 \
	forstr_cpp.f90 \
	forstrnps_cpp.f90 \
	fresid_cpp.f90 \
	fresidrsp_cpp.f90 \
	getcgqphase_cpp.f90 \
	init_e_field_vars_cpp.f90 \
	initberry_cpp.f90 \
	initorbmag_cpp.f90 \
	initmv_cpp.f90 \
	initro_cpp.f90 \
	ionion_realspace_cpp.f90 \
	ionion_surface_cpp.f90 \
	jellium_cpp.f90 \
	ks_ddiago_cpp.f90 \
	m_epjdos_cpp.f90 \
	mag_constr_cpp.f90 \
	mag_constr_e_cpp.f90 \
	make_efg_el_cpp.f90 \
	make_efg_ion_cpp.f90 \
	make_grad_berry_cpp.f90 \
	mkcore_inner_cpp.f90 \
	mkcore_paw_cpp.f90 \
	mkcore_wvl_cpp.f90 \
	mkgrid_fft_cpp.f90 \
	mklocl_cpp.f90 \
	mklocl_realspace_cpp.f90 \
	mklocl_recipspace_cpp.f90 \
	mklocl_wavelets_cpp.f90 \
	mkresi_cpp.f90 \
	mkrho_cpp.f90 \
	m_ioarr_cpp.f90 \
	m_iowf_cpp.f90 \
	mlwfovlp_cpp.f90 \
	mlwfovlp_proj_cpp.f90 \
	mlwfovlp_projpaw_cpp.f90 \
	mlwfovlp_pw_cpp.f90 \
	mlwfovlp_radial_cpp.f90 \
	mlwfovlp_seedname_cpp.f90 \
	mlwfovlp_setup_cpp.f90 \
	mlwfovlp_ylmfac_cpp.f90 \
	mlwfovlp_ylmfar_cpp.f90 \
	moddiel_cpp.f90 \
	msig_cpp.f90 \
	multipoles_out_cpp.f90 \
	m_eprenorms_cpp.f90 \
	m_plowannier_cpp.f90 \
	m_vcoul_cpp.f90 \
	newkpt_cpp.f90 \
	nres2vres_cpp.f90 \
	odamix_cpp.f90 \
	optics_vloc_cpp.f90 \
	partial_dos_fractions_cpp.f90 \
	partial_dos_fractions_paw_cpp.f90 \
	posdoppler_cpp.f90 \
	poslifetime_cpp.f90 \
	posratecore_cpp.f90 \
	prteigrs_cpp.f90 \
	prtene_cpp.f90 \
	prtefield_cpp.f90 \
	prtimg_cpp.f90 \
	prtrhomxmn_cpp.f90 \
	prtxf_cpp.f90 \
	prtxvf_cpp.f90 \
	rhotov_cpp.f90 \
	scprqt_cpp.f90 \
	setup_positron_cpp.f90 \
	setup1_cpp.f90 \
	setvtr_cpp.f90 \
	spatialchempot_cpp.f90 \
	spin_current_cpp.f90 \
	stress_cpp.f90 \
	strhar_cpp.f90 \
	sygrad_cpp.f90 \
	symrhg_cpp.f90 \
	uderiv_cpp.f90 \
	update_e_field_vars_cpp.f90 \
	vdw_dftd2_cpp.f90 \
	vdw_dftd3_cpp.f90 \
	vso_realspace_local_cpp.f90 \
	vtorhotf_cpp.f90 \
	waveformat_cpp.f90 \
	wf_mixing_cpp.f90 \
	wvl_initro_cpp.f90 \
	wvl_mkrho_cpp.f90 \
	interfaces_67_common_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/67_common
#
#COPYRIGHT
#Copyright (C) 2010-2018 ABINIT group
#This file is distributed under the terms of the
#GNU General Public License, see ~abinit/COPYING
#or http://www.gnu.org/copyleft/gpl.txt .
#
#THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
#To do that: config/scripts/abilint --dependencies . .

CLEANFILES += \
	interfaces_67_common.$(MODEXT) \
	m_epjdos.$(MODEXT) \
	m_eprenorms.$(MODEXT) \
	m_ioarr.$(MODEXT) \
	m_iowf.$(MODEXT) \
	m_plowannier.$(MODEXT) \
	m_vcoul.$(MODEXT) 

berryphase_new.$(OBJEXT): m_iowf.$(OBJEXT) 

calc_efg.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

cgwf.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

conducti_nc.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

conducti_paw.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

constrf.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

energy.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

etotfor.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

evdw_wannier.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

extraprho.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

forces.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

forstr.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

init_e_field_vars.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

interfaces_67_common.$(OBJEXT): m_epjdos.$(OBJEXT) 

mkcore_paw.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

mkcore_wvl.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

mklocl.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

mklocl_realspace.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

mklocl_wavelets.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

mkrho.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

mlwfovlp.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

mlwfovlp_proj.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

mlwfovlp_projpaw.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

partial_dos_fractions.$(OBJEXT): m_epjdos.$(OBJEXT) 

partial_dos_fractions_paw.$(OBJEXT): m_epjdos.$(OBJEXT) 

posdoppler.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

prtimg.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

rhotov.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

scprqt.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

setup_positron.$(OBJEXT): interfaces_67_common.$(OBJEXT) m_ioarr.$(OBJEXT) 

setvtr.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

spin_current.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

stress.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

uderiv.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

update_e_field_vars.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

vtorhotf.$(OBJEXT): interfaces_67_common.$(OBJEXT) 

EXTRA_DIST += _67_common_
