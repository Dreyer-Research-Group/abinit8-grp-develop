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
AM_FCFLAGS = @FPPFLAGS@ @FCFLAGS_FREEFORM@ @FCFLAGS_MODDIR@ @fcflags_opt_27_toolbox_oop@

AM_CPPFLAGS = \
	@src_12_hide_mpi_incs@ \
	@src_16_hideleave_incs@ \
	@src_11_memory_mpi_incs@ \
	@src_10_defs_incs@ \
	@src_14_hidewrite_incs@ \
	@src_incs_incs@ \
	@lib_netcdf_incs@ \
	@fallbacks_incs@ \
	@abi_extra_incs@ \
	@fc_mod_incs@

# Regular source files
lib27_toolbox_oop_srcs = \
	m_copy.F90 \
	m_iterators.F90 \
	m_mpiotk.F90 \
	m_nctk.F90 \
	interfaces_27_toolbox_oop.F90

# Source files depending on conditionals


# Library description
noinst_LIBRARIES = lib27_toolbox_oop.a

lib27_toolbox_oop_a_SOURCES= $(lib27_toolbox_oop_srcs)

CLEANFILES = \
	m_copy_cpp.f90 \
	m_iterators_cpp.f90 \
	m_mpiotk_cpp.f90 \
	m_nctk_cpp.f90 \
	interfaces_27_toolbox_oop_cpp.f90

EXTRA_DIST = abinit.src

EXTRA_DIST += abinit.dep

#Dependencies (inside the directory) of the directory ./src/27_toolbox_oop
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
	interfaces_27_toolbox_oop.$(MODEXT) \
	m_copy.$(MODEXT) \
	m_iterators.$(MODEXT) \
	m_mpiotk.$(MODEXT) \
	m_nctk.$(MODEXT) 

EXTRA_DIST += _27_toolbox_oop_
