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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_44_abitypes_defs@

AM_CPPFLAGS = \
	@src_32_util_incs@ \
	@src_27_toolbox_oop_incs@ \
	@src_42_libpaw_incs@ \
	@src_41_geometry_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_16_hideleave_incs@ \
	@src_14_hidewrite_incs@ \
	@src_incs_incs@ \
	@lib_gpu_incs@ \
	@lib_netcdf_incs@ \
	@fallbacks_incs@ \
	@abi_extra_incs@ \
	@fc_mod_incs@

# Regular source files
lib44_abitypes_defs_srcs = \
	defs_abitypes.F90 \
	m_args_gs.F90 \
	m_distribfft.F90 \
	m_dtset.F90 \
	m_efield.F90 \
	m_efmas_defs.F90 \
	m_energies.F90 \
	m_orbmag.F90 \
	m_results_gs.F90 \
	m_results_img.F90 \
	m_results_respfn.F90 \
	interfaces_44_abitypes_defs.F90

# Source files depending on conditionals


# Library description
noinst_LIBRARIES = lib44_abitypes_defs.a

lib44_abitypes_defs_a_SOURCES= $(lib44_abitypes_defs_srcs)

CLEANFILES = \
	defs_abitypes_cpp.f90 \
	m_args_gs_cpp.f90 \
	m_distribfft_cpp.f90 \
	m_dtset_cpp.f90 \
	m_efield_cpp.f90 \
	m_efmas_defs_cpp.f90 \
	m_energies_cpp.f90 \
	m_orbmag_cpp.f90 \
	m_results_gs_cpp.f90 \
	m_results_img_cpp.f90 \
	m_results_respfn_cpp.f90 \
	interfaces_44_abitypes_defs_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/44_abitypes_defs
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
	defs_abitypes.$(MODEXT) \
	interfaces_44_abitypes_defs.$(MODEXT) \
	m_args_gs.$(MODEXT) \
	m_distribfft.$(MODEXT) \
	m_dtset.$(MODEXT) \
	m_efield.$(MODEXT) \
	m_efmas_defs.$(MODEXT) \
	m_energies.$(MODEXT) \
	m_orbmag.$(MODEXT) \
	m_results_gs.$(MODEXT) \
	m_results_img.$(MODEXT) \
	m_results_respfn.$(MODEXT) 

defs_abitypes.$(OBJEXT): m_distribfft.$(OBJEXT) 

m_dtset.$(OBJEXT): defs_abitypes.$(OBJEXT) 

m_energies.$(OBJEXT): defs_abitypes.$(OBJEXT) 

m_results_gs.$(OBJEXT): m_energies.$(OBJEXT) 

m_results_img.$(OBJEXT): defs_abitypes.$(OBJEXT) m_energies.$(OBJEXT) m_results_gs.$(OBJEXT) 

m_results_respfn.$(OBJEXT): defs_abitypes.$(OBJEXT) 

EXTRA_DIST += _44_abitypes_defs_
