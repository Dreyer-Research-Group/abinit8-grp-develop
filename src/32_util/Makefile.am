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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_32_util@

AM_CPPFLAGS = \
	@src_17_libtetra_ext_incs@ \
	@src_16_hideleave_incs@ \
	@src_28_numeric_noabirule_incs@ \
	@src_01_linalg_ext_incs@ \
	@src_18_timing_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_14_hidewrite_incs@ \
	@src_02_clib_incs@ \
	@src_27_toolbox_oop_incs@ \
	@src_incs_incs@ \
	@lib_fft_incs@ \
	@lib_linalg_incs@ \
	@lib_etsf_io_incs@ \
	@lib_netcdf_incs@ \
	@fallbacks_incs@ \
	@abi_extra_incs@ \
	@fc_mod_incs@

# Regular source files
lib32_util_srcs = \
	appdig.F90 \
	fappnd.F90 \
	isfile.F90 \
	littlegroup_q.F90 \
	m_angles.F90 \
	m_atomdata.F90 \
	m_dummy.F90 \
	m_exp_mat.F90 \
	m_gwls_utility.F90 \
	m_matrix.F90 \
	m_optic_tools.F90 \
	m_pptools.F90 \
	m_special_funcs.F90 \
	m_wfutils.F90 \
	mati3det.F90 \
	mati3inv.F90 \
	matr3inv.F90 \
	mknormpath.F90 \
	proc_distrb_cycle.F90 \
	radsintr.F90 \
	radsmear.F90 \
	smatrix.F90 \
	status.F90 \
	interfaces_32_util.F90

# Source files depending on conditionals


# Library description
noinst_LIBRARIES = lib32_util.a

lib32_util_a_SOURCES= $(lib32_util_srcs)

CLEANFILES = \
	appdig_cpp.f90 \
	fappnd_cpp.f90 \
	isfile_cpp.f90 \
	littlegroup_q_cpp.f90 \
	m_angles_cpp.f90 \
	m_atomdata_cpp.f90 \
	m_dummy_cpp.f90 \
	m_exp_mat_cpp.f90 \
	m_gwls_utility_cpp.f90 \
	m_matrix_cpp.f90 \
	m_optic_tools_cpp.f90 \
	m_pptools_cpp.f90 \
	m_special_funcs_cpp.f90 \
	m_wfutils_cpp.f90 \
	mati3det_cpp.f90 \
	mati3inv_cpp.f90 \
	matr3inv_cpp.f90 \
	mknormpath_cpp.f90 \
	proc_distrb_cycle_cpp.f90 \
	radsintr_cpp.f90 \
	radsmear_cpp.f90 \
	smatrix_cpp.f90 \
	status_cpp.f90 \
	interfaces_32_util_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/32_util
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
	interfaces_32_util.$(MODEXT) \
	m_angles.$(MODEXT) \
	m_atomdata.$(MODEXT) \
	m_dummy_tests.$(MODEXT) \
	m_exp_mat.$(MODEXT) \
	m_gwls_utility.$(MODEXT) \
	m_matrix.$(MODEXT) \
	m_optic_tools.$(MODEXT) \
	m_pptools.$(MODEXT) \
	m_special_funcs.$(MODEXT) \
	m_wfutils.$(MODEXT) 

m_angles.$(OBJEXT): m_special_funcs.$(OBJEXT) 

test_radsintr.$(OBJEXT): interfaces_32_util.$(OBJEXT) 

check_PROGRAMS = \
	test_radsintr

test_radsintr_SOURCES = test_radsintr.F90
test_radsintr_CPPFLAGS = @lib_netcdf_incs@ @lib_etsf_io_incs@ @lib_linalg_incs@ -I$(top_srcdir)/src/incs
test_radsintr_LDADD = \
	lib32_util.a \
	../27_toolbox_oop/lib27_toolbox_oop.a \
	../28_numeric_noabirule/lib28_numeric_noabirule.a \
	../16_hideleave/lib16_hideleave.a \
	../14_hidewrite/lib14_hidewrite.a \
	../12_hide_mpi/lib12_hide_mpi.a \
	../11_memory_mpi/lib11_memory_mpi.a \
	../10_dumpinfo/lib10_dumpinfo.a \
	../10_defs/lib10_defs.a \
	../01_linalg_ext/lib01_linalg_ext.a \
	@lib_etsf_io_libs@ \
	@lib_netcdf_libs@ \
	@lib_linalg_libs@

check-local:
	./test_radsintr >test_radsintr.log 2>&1; grep '^TEST FAILED' test_radsintr.log && echo 'Unit test test_radsintr FAILED' || echo 'Unit test test_radsintr OK'

CLEANFILES += \
	test_radsintr.log \
	test_radsintr.tmp

EXTRA_DIST += _32_util_
