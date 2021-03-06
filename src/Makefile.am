#
# Makefile for ABINIT                                      -*- Automake -*-
# Generated by make-makefiles-inter on 2021/06/01 15:12:57 +0000

#
# IMPORTANT NOTE
#
# Any manual change to this file will systematically be overwritten.
# Please modify the make-makefiles-inter script or its config file instead.
#

if DO_BUILD_01_GSL_EXT
 core_01_gsl_ext = 01_gsl_ext
else
 core_01_gsl_ext =
endif
if DO_BUILD_01_TRIQS_EXT
 core_01_triqs_ext = 01_triqs_ext
else
 core_01_triqs_ext =
endif
if DO_BUILD_15_GPU_TOOLBOX
 core_15_gpu_toolbox = 15_gpu_toolbox
else
 core_15_gpu_toolbox =
endif
if DO_BUILD_52_MANAGE_CUDA
 core_52_manage_cuda = 52_manage_cuda
else
 core_52_manage_cuda =
endif
if DO_BUILD_65_LOTF_BASE
 core_65_lotf_base = 65_lotf_base
else
 core_65_lotf_base =
endif
if DO_BUILD_68_LOTF
 core_68_lotf = 68_lotf
else
 core_68_lotf =
endif

SUBDIRS = \
	incs \
	mods \
	$(core_01_gsl_ext) \
	01_linalg_ext \
	01_macroavnew_ext \
	$(core_01_triqs_ext) \
	02_clib \
	10_defs \
	10_dumpinfo \
	11_memory_mpi \
	11_qespresso_ext \
	12_hide_mpi \
	14_hidewrite \
	$(core_15_gpu_toolbox) \
	16_hideleave \
	17_libtetra_ext \
	18_timing \
	20_datashare \
	21_hashfuncs \
	27_toolbox_oop \
	28_numeric_noabirule \
	29_xgTools \
	30_diago \
	32_util \
	41_geometry \
	41_xc_lowlevel \
	42_libpaw \
	42_nlstrain \
	42_parser \
	43_ptgroups \
	43_wvl_wrappers \
	44_abitypes_defs \
	45_geomoptim \
	49_gw_toolbox_oop \
	51_manage_mpi \
	52_fft_mpi_noabirule \
	$(core_52_manage_cuda) \
	53_ffts \
	53_spacepar \
	54_abiutil \
	56_io_mpi \
	56_mixing \
	56_recipspace \
	56_xc \
	57_iopsp_parser \
	57_iovars \
	59_ionetcdf \
	61_occeig \
	62_cg_noabirule \
	62_ctqmc \
	62_iowfdenpot \
	62_poisson \
	62_wvl_wfs \
	63_bader \
	64_atompaw \
	64_psp \
	$(core_65_lotf_base) \
	65_paw \
	66_nonlocal \
	66_wfs \
	67_common \
	68_dmft \
	$(core_68_lotf) \
	68_recursion \
	68_rsprc \
	69_wfdesc \
	70_gw \
	71_bse \
	72_response \
	77_ddb \
	77_suscep \
	78_effpot \
	78_eph \
	79_seqpar_mpi \
	80_tdep \
	83_cut3d \
	84_fold2Bloch \
	94_scfcv \
	95_drive \
	98_main
if DO_BUILD_EXPORTS
 SUBDIRS += libs
endif

# ---------------------------------------------------------------------------- #

# Enforce ABINIT Coding Style (the so-called ABIRULES)
abirules:
	@echo 'Reporting possible errors in the abirules procedure' > tmp-abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 01_gsl_ext >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 01_triqs_ext >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 41_geometry >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 42_libpaw >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 42_nlstrain >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 42_parser >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 43_wvl_wrappers >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 45_geomoptim >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 53_spacepar >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 56_mixing >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 56_recipspace >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 57_iopsp_parser >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 59_ionetcdf >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 61_occeig >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 62_wvl_wfs >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 63_bader >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 66_wfs >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 68_rsprc >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 72_response >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 77_suscep >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 78_effpot >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 78_eph >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 83_cut3d >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 84_fold2Bloch >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 94_scfcv >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 95_drive >> abirules.log
	@$(PERL) ../developers/maintainers/abirules.pl -d 98_main >> abirules.log
	cd .. ; ./config/scripts/abilint . . > src/abilint.log ; cd src
	@grep 'Error' abirules.log 

# ---------------------------------------------------------------------------- #

# Help
help:
	@cat $(top_srcdir)/doc/help_make/help_make_src.txt

# ---------------------------------------------------------------------------- #

# Creation of parent subroutines
parents:
	$(PERL) $(top_srcdir)/developers/maintainers/parents.pl -s
	rm -f *.old */*.old

# ---------------------------------------------------------------------------- #

# Beautify ABINIT sources (indentation)
# Note : the _ext, _noabirule should not treated because they are not expected to follow the abirules
# Note : 02_clib only includes C files
# Note : 43_ptgroup should not be treated because the routines have been generated automatically
# Note : 69_wfdesc, 70_gw and 71_bse  should not be treated due a special request of Matteo
# Note : the "Warning" occurences have been eliminated on 10 May 2011.
# Note : the "ERROR" occurences have been eliminated on 25 December 2011, also documenting the solution for each type of problem !!
# Note : 26 Feb 2012 : error messages have been improved, recursive has been allowed.
# Should write the documentation, so that to imitate the "make parents" treatment on shiva ...
# The treatment of $omp is still incorrect !!
#
# 1) ERROR(41_geometry/prt_cif.F90): found end statement at line 273 for '' subroutine symrel2string '' before executable section
#      The structuration between the declarations and the executable section was missing
# 2) ERROR(67_common/mklocl_realspace.F90): found end statement at line 849 for '' subroutine createIonicPotential_new '' instead of '' if construct ''
#      There was a (forbidden) on-line if ... call :
#       if (condition)  call subroutine
#      Solution : replace by
#        if (condition) then
#            call subroutine
#        endif
# 3) ERROR(67_common/rhotov.F90): found end statement at line 412 for '' subroutine rhotov '' instead of '' if construct ''
#      There was a (forbidden) on-line secondary if construct
#        command ; if( condition) command2
#      Solution : replace by
#        command
#        if(condition) command2
# 4) Similar message of error might be caused by sections of the code that start from the first column !! 
#      Solution : by hand, insert a blank. Abiauty will do the rest ...
# 
# When happening in a module, ignore the problems (modules are not yet treated by abiauty).
# 
# 1) ERROR, executable section not found
# 2) ERROR(10_defs/defs_datatypes.F90): found type at line 99 while expecting program/subroutine
# 3) ERROR(10_defs/defs_datatypes.F90): end of file 10_defs/defs_datatypes.F90 hit prematurely while expecting:
#  program/subroutine/function statement
# 4) ERROR(10_defs/defs_wannier90.F90): found end statement at line 62 for '' subroutine wannier_setup '' before executable section
abiauty:
	@echo 'Reporting possible errors in the abiauty procedure' > tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/10_defs >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/11_memory_mpi >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/12_hide_mpi >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/14_hidewrite >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/15_gpu_toolbox >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/16_hideleave >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/18_timing >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/20_datashare >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/27_toolbox_oop >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/32_util >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/41_geometry >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/41_xc_lowlevel >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/42_libpaw >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/42_nlstrain >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/42_parser >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/43_wvl_wrappers >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/44_abitypes_defs >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/45_geomoptim >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/49_gw_toolbox_oop >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/51_manage_mpi >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/52_manage_cuda >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/53_ffts >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/53_spacepar >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/54_abiutil >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/56_io_mpi >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/56_mixing >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/56_recipspace >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/56_xc >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/57_iopsp_parser >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/57_iovars >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/59_ionetcdf >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/61_occeig >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/62_ctqmc >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/62_iowfdenpot >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/62_poisson >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/62_wvl_wfs >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/63_bader >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/64_atompaw >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/64_psp >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/65_lotf_base >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/65_paw >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/66_wfs >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/66_nonlocal >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/67_common >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/68_dmft >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/68_lotf >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/68_recursion >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/68_rsprc >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/72_response >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/77_ddb >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/77_suscep >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/79_seqpar_mpi >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/83_cut3d >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/84_fold2Bloch >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/94_scfcv >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/95_drive >> tmp-abiauty.log
	@$(PERL) $(top_srcdir)/developers/maintainers/abiauty.pl -d $(top_srcdir)/src/98_main >> tmp-abiauty.log
	grep ERROR tmp-abiauty.log > error_make_abiauty  || true

# ---------------------------------------------------------------------------- #
# Note concerning abirules (as of 5.8.0)
# The order of config.h vs abi_common.h is not preserved.
# This should be corrected by hand. It is to be done, at present,
# in get_irredg, merge_kgirr, setshells, sigma

# Thorough clean-up
maintainer-clean-local:
	rm -f tmp-abiauty.log tmp-abirules.log
	rm -f abirules.log
	find . -name '*.old' -exec rm {} \;
	find . -name '*.orig' -exec rm {} \;
	find . -name '*.rej' -exec rm {} \;
