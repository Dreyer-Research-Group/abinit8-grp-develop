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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_62_wvl_wfs@

AM_CPPFLAGS = \
	@src_16_hideleave_incs@ \
	@src_44_abitypes_defs_incs@ \
	@src_42_libpaw_incs@ \
	@src_51_manage_mpi_incs@ \
	@src_43_wvl_wrappers_incs@ \
	@src_12_hide_mpi_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_14_hidewrite_incs@ \
	@src_41_geometry_incs@ \
	@src_incs_incs@ \
	@lib_bigdft_incs@ \
	@lib_netcdf_incs@ \
	@fallbacks_incs@ \
	@abi_extra_incs@ \
	@fc_mod_incs@

# Regular source files
lib62_wvl_wfs_srcs = \
	wvl_hpsitopsi.F90 \
	wvl_nl_gradient.F90 \
	wvl_psitohpsi.F90 \
	wvl_rwwf.F90 \
	wvl_tail_corrections.F90 \
	interfaces_62_wvl_wfs.F90

# Source files depending on conditionals


# Library description
noinst_LIBRARIES = lib62_wvl_wfs.a

lib62_wvl_wfs_a_SOURCES= $(lib62_wvl_wfs_srcs)

CLEANFILES = \
	wvl_hpsitopsi_cpp.f90 \
	wvl_nl_gradient_cpp.f90 \
	wvl_psitohpsi_cpp.f90 \
	wvl_rwwf_cpp.f90 \
	wvl_tail_corrections_cpp.f90 \
	interfaces_62_wvl_wfs_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/62_wvl_wfs
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
	interfaces_62_wvl_wfs.$(MODEXT) 

EXTRA_DIST += _62_wvl_wfs_
