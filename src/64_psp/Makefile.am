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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_64_psp@

AM_CPPFLAGS = \
	@src_32_util_incs@ \
	@src_28_numeric_noabirule_incs@ \
	@src_44_abitypes_defs_incs@ \
	@src_42_libpaw_incs@ \
	@src_51_manage_mpi_incs@ \
	@src_43_wvl_wrappers_incs@ \
	@src_16_hideleave_incs@ \
	@src_56_recipspace_incs@ \
	@src_57_iopsp_parser_incs@ \
	@src_27_toolbox_oop_incs@ \
	@src_42_parser_incs@ \
	@src_11_qespresso_ext_incs@ \
	@src_18_timing_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_14_hidewrite_incs@ \
	@src_57_iovars_incs@ \
	@src_53_ffts_incs@ \
	@src_incs_incs@ \
	@lib_bigdft_incs@ \
	@lib_psml_incs@ \
	@lib_netcdf_incs@ \
	@lib_libxc_incs@ \
	@lib_math_incs@ \
	@fallbacks_incs@ \
	@abi_extra_incs@ \
	@fc_mod_incs@

# Regular source files
lib64_psp_srcs = \
	atm2fft.F90 \
	cc_derivatives.F90 \
	dfpt_atm2fft.F90 \
	m_psps.F90 \
	pspatm.F90 \
	pspini.F90 \
	psp1cc.F90 \
	psp1in.F90 \
	psp2in.F90 \
	psp2lo.F90 \
	psp3in.F90 \
	psp4cc.F90 \
	psp5in.F90 \
	psp5lo.F90 \
	psp5nl.F90 \
	psp6in.F90 \
	psp8cc.F90 \
	psp8in.F90 \
	psp8lo.F90 \
	psp8nl.F90 \
	psp9in.F90 \
	psp10in.F90 \
	psp10nl.F90 \
	psp11lo.F90 \
	psp11nl.F90 \
	psp_dump_outputs.F90 \
	upf2abinit.F90 \
	interfaces_64_psp.F90

# Source files depending on conditionals
if DO_BUILD_PSML
  lib64_psp_srcs += \
    psp9cc.F90
endif

# Library description
noinst_LIBRARIES = lib64_psp.a

lib64_psp_a_SOURCES= $(lib64_psp_srcs)

CLEANFILES = \
	atm2fft_cpp.f90 \
	cc_derivatives_cpp.f90 \
	dfpt_atm2fft_cpp.f90 \
	m_psps_cpp.f90 \
	pspatm_cpp.f90 \
	pspini_cpp.f90 \
	psp1cc_cpp.f90 \
	psp1in_cpp.f90 \
	psp2in_cpp.f90 \
	psp2lo_cpp.f90 \
	psp3in_cpp.f90 \
	psp4cc_cpp.f90 \
	psp5in_cpp.f90 \
	psp5lo_cpp.f90 \
	psp5nl_cpp.f90 \
	psp6in_cpp.f90 \
	psp8cc_cpp.f90 \
	psp8in_cpp.f90 \
	psp8lo_cpp.f90 \
	psp8nl_cpp.f90 \
	psp9in_cpp.f90 \
	psp10in_cpp.f90 \
	psp10nl_cpp.f90 \
	psp11lo_cpp.f90 \
	psp11nl_cpp.f90 \
	psp_dump_outputs_cpp.f90 \
	upf2abinit_cpp.f90 \
	interfaces_64_psp_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/64_psp
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
	interfaces_64_psp.$(MODEXT) \
	m_psps.$(MODEXT) 

psp10in.$(OBJEXT): interfaces_64_psp.$(OBJEXT) 

psp1cc.$(OBJEXT): interfaces_64_psp.$(OBJEXT) 

psp1in.$(OBJEXT): interfaces_64_psp.$(OBJEXT) 

psp2in.$(OBJEXT): interfaces_64_psp.$(OBJEXT) 

psp3in.$(OBJEXT): interfaces_64_psp.$(OBJEXT) 

psp5in.$(OBJEXT): interfaces_64_psp.$(OBJEXT) 

psp6in.$(OBJEXT): interfaces_64_psp.$(OBJEXT) 

psp8in.$(OBJEXT): interfaces_64_psp.$(OBJEXT) m_psps.$(OBJEXT) 

psp9in.$(OBJEXT): interfaces_64_psp.$(OBJEXT) m_psps.$(OBJEXT) 

pspatm.$(OBJEXT): interfaces_64_psp.$(OBJEXT) m_psps.$(OBJEXT) 

pspini.$(OBJEXT): interfaces_64_psp.$(OBJEXT) m_psps.$(OBJEXT) 

upf2abinit.$(OBJEXT): interfaces_64_psp.$(OBJEXT) 

EXTRA_DIST += _64_psp_
