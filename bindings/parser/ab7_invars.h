#ifndef AB7_INVARS_H
#define AB7_INVARS_H

#include <stdlib.h>

#include "ab7_base.h"

/**
 * Ab7InvarsTypes:
 * @_INT_SCALAR: a 32 bits integer.
 * @_DOUBLE_SCALAR: a 64 bits float.
 * @_INT_ARRAY: an array of 32 bits integers.
 * @_DOUBLE_ARRAY: an array of 64 bits floats.
 *
 * The possible types of the attributes of datasets.
 */
typedef enum
  {
    _INT_SCALAR,
    _INT_ARRAY,
    _DOUBLE_SCALAR,
    _DOUBLE_ARRAY,
    _OTHER
  } Ab7InvarsTypes;

/* This file has been automatically generated, do not modify. */
typedef enum
{
  AB7_INVARS_SYMCHI         ,  /* _INT_SCALAR     */
  AB7_INVARS_PVELMAX        ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_BS_INTERP_RL_NB,  /* _INT_SCALAR     */
  AB7_INVARS_SLABWSRAD      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_RFSTRS         ,  /* _INT_SCALAR     */
  AB7_INVARS_EXCHN2N3D      ,  /* _INT_SCALAR     */
  AB7_INVARS_D3E_PERT3_ELFD ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTQMC_N      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_WVL_NPRCCG     ,  /* _INT_SCALAR     */
  AB7_INVARS_GETDDB         ,  /* _INT_SCALAR     */
  AB7_INVARS_GWGAMMA        ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_NGPTS   ,  /* _INT_SCALAR     */
  AB7_INVARS_DDB_NGQPT      ,  /* _INT_ARRAY      */
  AB7_INVARS_GW_TOLDFEIG    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_USEEXEXCH      ,  /* _INT_SCALAR     */
  AB7_INVARS_TOLMXF         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_LOTF_VERSION   ,  /* _INT_SCALAR     */
  AB7_INVARS_EFMAS_N_DIRS   ,  /* _INT_SCALAR     */
  AB7_INVARS_CD_FRQIM_METHOD,  /* _INT_SCALAR     */
  AB7_INVARS_PRT1DM         ,  /* _INT_SCALAR     */
  AB7_INVARS_USEPOTZERO     ,  /* _INT_SCALAR     */
  AB7_INVARS_PREPANL        ,  /* _INT_SCALAR     */
  AB7_INVARS_GWLS_KMAX_ANALYTIC,  /* _INT_SCALAR     */
  AB7_INVARS_NSPINOR        ,  /* _INT_SCALAR     */
  AB7_INVARS_RFATPOL        ,  /* _INT_ARRAY      */
  AB7_INVARS_DMFTCTQMC_GMOVE,  /* _INT_SCALAR     */
  AB7_INVARS_NDTSET         ,  /* _INT_SCALAR     */
  AB7_INVARS_USEPAWU        ,  /* _INT_SCALAR     */
  AB7_INVARS_MPW            ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWCROSS       ,  /* _INT_SCALAR     */
  AB7_INVARS_OCCOPT         ,  /* _INT_SCALAR     */
  AB7_INVARS_GETDDK         ,  /* _INT_SCALAR     */
  AB7_INVARS_IMGMOV         ,  /* _INT_SCALAR     */
  AB7_INVARS_POSTOLDFF      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_VDW_DF_TOLERANCE,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_POSTOLDFE      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_BOXCENTER      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_BUILTINTEST    ,  /* _INT_SCALAR     */
  AB7_INVARS_TOLMXDE        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_BXCTMINDG      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_QPTDM          ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_TD_MEXCIT      ,  /* _INT_SCALAR     */
  AB7_INVARS_GA_FITNESS     ,  /* _INT_SCALAR     */
  AB7_INVARS_NCTIME         ,  /* _INT_SCALAR     */
  AB7_INVARS_FRZFERMI       ,  /* _INT_SCALAR     */
  AB7_INVARS_HYB_MIXING     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_RFPHON         ,  /* _INT_SCALAR     */
  AB7_INVARS_PLOWAN_BANDI   ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTDEN         ,  /* _INT_SCALAR     */
  AB7_INVARS_GWPARA         ,  /* _INT_SCALAR     */
  AB7_INVARS_IRDBSRESO      ,  /* _INT_SCALAR     */
  AB7_INVARS_RECPTROTT      ,  /* _INT_SCALAR     */
  AB7_INVARS_RED_EFIELDBAR  ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_JDTSET         ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_HAYDOCK_NITER,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTQMC_L      ,  /* _INT_SCALAR     */
  AB7_INVARS_MEP_SOLVER     ,  /* _INT_SCALAR     */
  AB7_INVARS_NSTEP          ,  /* _INT_SCALAR     */
  AB7_INVARS_IOMODE         ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTGSR         ,  /* _INT_SCALAR     */
  AB7_INVARS_USE_GPU_CUDA   ,  /* _INT_SCALAR     */
  AB7_INVARS_ALGALCH        ,  /* _INT_ARRAY      */
  AB7_INVARS_USE_GEMM_NONLOP,  /* _INT_SCALAR     */
  AB7_INVARS_NC_XCCC_GSPACE ,  /* _INT_SCALAR     */
  AB7_INVARS_KPTRLATT       ,  /* _INT_ARRAY      */
  AB7_INVARS_BS_ALGORITHM   ,  /* _INT_SCALAR     */
  AB7_INVARS_GWMEM          ,  /* _INT_SCALAR     */
  AB7_INVARS_NLOALG         ,  /* _INT_ARRAY      */
  AB7_INVARS_NPKPT          ,  /* _INT_SCALAR     */
  AB7_INVARS_GPU_DEVICES    ,  /* _INT_ARRAY      */
  AB7_INVARS_PRTXML         ,  /* _INT_SCALAR     */
  AB7_INVARS_USEKDEN        ,  /* _INT_SCALAR     */
  AB7_INVARS_NPPERT         ,  /* _INT_SCALAR     */
  AB7_INVARS_XRED_ORIG      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_GWLS_BAND_INDEX,  /* _INT_SCALAR     */
  AB7_INVARS_PH_NQPATH      ,  /* _INT_SCALAR     */
  AB7_INVARS_KPTNS          ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NBERRY         ,  /* _INT_SCALAR     */
  AB7_INVARS_CHKSYMBREAK    ,  /* _INT_SCALAR     */
  AB7_INVARS_FREQREMIN      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_MGFFTDG        ,  /* _INT_SCALAR     */
  AB7_INVARS_EXCHMIX        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_D3E_PERT3_ATPOL,  /* _INT_ARRAY      */
  AB7_INVARS_TFKINFUNC      ,  /* _INT_SCALAR     */
  AB7_INVARS_PTCHARGE       ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_KPTNRM         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_IRDBSEIG       ,  /* _INT_SCALAR     */
  AB7_INVARS_SLABZBEG       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GW_CUSTOMNFREQSP,  /* _INT_SCALAR     */
  AB7_INVARS_RED_EFIELD     ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_GW_SCTYPE      ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_INTERP_KMULT,  /* _INT_ARRAY      */
  AB7_INVARS_GWLS_DIEL_MODEL,  /* _INT_SCALAR     */
  AB7_INVARS_VIS            ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_CD_IMFRQS      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_USERRE         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_DMFTCTQMC_CORREL,  /* _INT_SCALAR     */
  AB7_INVARS_PRTWF_FULL     ,  /* _INT_SCALAR     */
  AB7_INVARS_NPWWFN         ,  /* _INT_SCALAR     */
  AB7_INVARS_EFMAS_NTHETA   ,  /* _INT_SCALAR     */
  AB7_INVARS_TIM1REV        ,  /* _INT_SCALAR     */
  AB7_INVARS_MFFMEM         ,  /* _INT_SCALAR     */
  AB7_INVARS_DFIELD         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_DMFTCTQMC_GRNNS,  /* _INT_SCALAR     */
  AB7_INVARS_BS_EXCHANGE_TERM,  /* _INT_SCALAR     */
  AB7_INVARS_CHARGE         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GETWFKFINE     ,  /* _INT_SCALAR     */
  AB7_INVARS_EPH_TASK       ,  /* _INT_SCALAR     */
  AB7_INVARS_MK1MEM         ,  /* _INT_SCALAR     */
  AB7_INVARS_HYB_RANGE_FOCK ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_VDW_DF_DAMIN   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_D3E_PERT2_DIR  ,  /* _INT_ARRAY      */
  AB7_INVARS_PRTKPT         ,  /* _INT_SCALAR     */
  AB7_INVARS_FREQIM_ALPHA   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_VDW_DF_DCUT    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NPFFT          ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTVCLMB       ,  /* _INT_SCALAR     */
  AB7_INVARS_QMASS          ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_RF2_PERT2_DIR  ,  /* _INT_ARRAY      */
  AB7_INVARS_ECUTSIGX       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GA_OPT_PERCENT ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_POLCEN         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_RECRCUT        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_WVL_HGRID      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_USEREC         ,  /* _INT_SCALAR     */
  AB7_INVARS_GETXCART       ,  /* _INT_SCALAR     */
  AB7_INVARS_SPMETH         ,  /* _INT_SCALAR     */
  AB7_INVARS_NPVEL          ,  /* _INT_SCALAR     */
  AB7_INVARS_MKMEM          ,  /* _INT_SCALAR     */
  AB7_INVARS_QPRTRB         ,  /* _INT_ARRAY      */
  AB7_INVARS_GWLS_NPT_GAUSS_QUAD,  /* _INT_SCALAR     */
  AB7_INVARS_CINEB_START    ,  /* _INT_SCALAR     */
  AB7_INVARS_GETDEN         ,  /* _INT_SCALAR     */
  AB7_INVARS_AUXC_SCAL      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_RECEFERMI      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_IRDPAWDEN      ,  /* _INT_SCALAR     */
  AB7_INVARS_PLOWAN_IATOM   ,  /* _INT_ARRAY      */
  AB7_INVARS_VACWIDTH       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_EPH_NGQPT_FINE ,  /* _INT_ARRAY      */
  AB7_INVARS_USEDMATPU      ,  /* _INT_SCALAR     */
  AB7_INVARS_TD_MAXENE      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_FOCKDOWNSAMPLING,  /* _INT_ARRAY      */
  AB7_INVARS_GW_INVALID_FREQ,  /* _INT_SCALAR     */
  AB7_INVARS_PARAL_KGB      ,  /* _INT_SCALAR     */
  AB7_INVARS_STRING_ALGO    ,  /* _INT_SCALAR     */
  AB7_INVARS_ADPIMD         ,  /* _INT_SCALAR     */
  AB7_INVARS_GET1WF         ,  /* _INT_SCALAR     */
  AB7_INVARS_IRDDEN         ,  /* _INT_SCALAR     */
  AB7_INVARS_ZCUT           ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GETBSEIG       ,  /* _INT_SCALAR     */
  AB7_INVARS_GETXRED        ,  /* _INT_SCALAR     */
  AB7_INVARS_BMASS          ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GWLS_KMAX_POLES,  /* _INT_SCALAR     */
  AB7_INVARS_W90PRTUNK      ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTCTQMC_MRKA ,  /* _INT_SCALAR     */
  AB7_INVARS_PLOWAN_COMPUTE ,  /* _INT_SCALAR     */
  AB7_INVARS_TYPAT          ,  /* _INT_ARRAY      */
  AB7_INVARS_RATSPH         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_ICUTCOUL       ,  /* _INT_SCALAR     */
  AB7_INVARS_MEP_MXSTEP     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PRTSPCUR       ,  /* _INT_SCALAR     */
  AB7_INVARS_USEYLM         ,  /* _INT_SCALAR     */
  AB7_INVARS_MACRO_UJ       ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFT_RSLF      ,  /* _INT_SCALAR     */
  AB7_INVARS_NPIMAGE        ,  /* _INT_SCALAR     */
  AB7_INVARS_GETCELL        ,  /* _INT_SCALAR     */
  AB7_INVARS_MAX_NCPUS      ,  /* _INT_SCALAR     */
  AB7_INVARS_MAXNSYM        ,  /* _INT_SCALAR     */
  AB7_INVARS_TOLSYM         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_RECTOLDEN      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_KSSFORM        ,  /* _INT_SCALAR     */
  AB7_INVARS_USERRD         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_USERRA         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_USERRC         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_DDAMP          ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_EFFMASS_FREE   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PRTPHBANDS     ,  /* _INT_SCALAR     */
  AB7_INVARS_D3E_PERT3_PHON ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFT_MXSF      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_WFOPTALG       ,  /* _INT_SCALAR     */
  AB7_INVARS_TFW_TOLDFE     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_VCUTGEO        ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_TNONS          ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_EFMAS_DIRS     ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_IRDWFK         ,  /* _INT_SCALAR     */
  AB7_INVARS_IRD1DEN        ,  /* _INT_SCALAR     */
  AB7_INVARS_NTIMIMAGE      ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWLCUTD       ,  /* _INT_SCALAR     */
  AB7_INVARS_NBAND          ,  /* _INT_ARRAY      */
  AB7_INVARS_EXTRAPWF       ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWUJV         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_DMFTQMC_SEED   ,  /* _INT_SCALAR     */
  AB7_INVARS_ENUNIT         ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_NFRAG      ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTNEST        ,  /* _INT_SCALAR     */
  AB7_INVARS_D3E_PERT2_PHON ,  /* _INT_SCALAR     */
  AB7_INVARS_IRDWFQ         ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWSTGYLM      ,  /* _INT_SCALAR     */
  AB7_INVARS_RF2_PERT1_DIR  ,  /* _INT_ARRAY      */
  AB7_INVARS_IPRCEL         ,  /* _INT_SCALAR     */
  AB7_INVARS_IATFIX         ,  /* _INT_ARRAY      */
  AB7_INVARS_DMFTCTQMC_MEAS ,  /* _INT_SCALAR     */
  AB7_INVARS_RESTARTXF      ,  /* _INT_SCALAR     */
  AB7_INVARS_PH_NGQPT       ,  /* _INT_ARRAY      */
  AB7_INVARS_IRDDDK         ,  /* _INT_SCALAR     */
  AB7_INVARS_SPINMAGNTARGET ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_CHEMPOT        ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PRTPSPS        ,  /* _INT_SCALAR     */
  AB7_INVARS_AUTOPARAL      ,  /* _INT_SCALAR     */
  AB7_INVARS_DYNIMAGE       ,  /* _INT_ARRAY      */
  AB7_INVARS_SPGROUP        ,  /* _INT_SCALAR     */
  AB7_INVARS_NPWSIGX        ,  /* _INT_SCALAR     */
  AB7_INVARS_VEL_ORIG       ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NOMEGASF       ,  /* _INT_SCALAR     */
  AB7_INVARS_IATSPH         ,  /* _INT_ARRAY      */
  AB7_INVARS_PW_UNBAL_THRESH,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NLINE          ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTQMC_THERM  ,  /* _INT_SCALAR     */
  AB7_INVARS_KPT            ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_BDBERRY        ,  /* _INT_ARRAY      */
  AB7_INVARS_RFMETH         ,  /* _INT_SCALAR     */
  AB7_INVARS_USERIA         ,  /* _INT_SCALAR     */
  AB7_INVARS_NFREQSP        ,  /* _INT_SCALAR     */
  AB7_INVARS_USERIC         ,  /* _INT_SCALAR     */
  AB7_INVARS_USERID         ,  /* _INT_SCALAR     */
  AB7_INVARS_USERIE         ,  /* _INT_SCALAR     */
  AB7_INVARS_EFMAS_DEG      ,  /* _INT_SCALAR     */
  AB7_INVARS_GW_FREQSP      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_RANDOM_ATPOS   ,  /* _INT_SCALAR     */
  AB7_INVARS_NTYPALCH       ,  /* _INT_SCALAR     */
  AB7_INVARS_LOCALRDWF      ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTDOSM        ,  /* _INT_SCALAR     */
  AB7_INVARS_NIMAGE         ,  /* _INT_SCALAR     */
  AB7_INVARS_MDTEMP         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NP_SLK         ,  /* _INT_SCALAR     */
  AB7_INVARS_CD_HALFWAY_FREQ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_DIEMIX         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_LOTF_NNEIGX    ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_NSMOOTH ,  /* _INT_SCALAR     */
  AB7_INVARS_DENSTY         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PRTVPSP        ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTEBANDS      ,  /* _INT_SCALAR     */
  AB7_INVARS_MAXESTEP       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_EFMAS_CALC_DIRS,  /* _INT_SCALAR     */
  AB7_INVARS_GWENCOMP       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_RECTESTEG      ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTNABLA       ,  /* _INT_SCALAR     */
  AB7_INVARS_NTYPAT         ,  /* _INT_SCALAR     */
  AB7_INVARS_ICOULOMB       ,  /* _INT_SCALAR     */
  AB7_INVARS_IPRCFC         ,  /* _INT_SCALAR     */
  AB7_INVARS_NPHF           ,  /* _INT_SCALAR     */
  AB7_INVARS_CORECS         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_DMFT_TOLLC     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GETHAYDOCK     ,  /* _INT_SCALAR     */
  AB7_INVARS_ASR            ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFT_ENTROPY   ,  /* _INT_SCALAR     */
  AB7_INVARS_GW_NSTEP       ,  /* _INT_SCALAR     */
  AB7_INVARS_EPH_EXTRAEL    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_ECUTSM         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_EPH_TRANSPORT  ,  /* _INT_SCALAR     */
  AB7_INVARS_FREQREMAX      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NBDBUF         ,  /* _INT_SCALAR     */
  AB7_INVARS_EINTERP        ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NATPAWU        ,  /* _INT_SCALAR     */
  AB7_INVARS_DIELAM         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_EFIELD         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_GW_FRQIM_INZGRID,  /* _INT_SCALAR     */
  AB7_INVARS_ISCF           ,  /* _INT_SCALAR     */
  AB7_INVARS_BDEIGRF        ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWUJAT        ,  /* _INT_SCALAR     */
  AB7_INVARS_IXC            ,  /* _INT_SCALAR     */
  AB7_INVARS_DELAYPERM      ,  /* _INT_SCALAR     */
  AB7_INVARS_XC_TB09_C      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NSCFORDER      ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWOVLP        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PLOWAN_IT      ,  /* _INT_ARRAY      */
  AB7_INVARS_ATVSHIFT       ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PRTPMP         ,  /* _INT_SCALAR     */
  AB7_INVARS_CPUS           ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_VDW_DF_QCUT    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_ECUTEPS        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NZCHEMPOT      ,  /* _INT_SCALAR     */
  AB7_INVARS_GW_QLWL        ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PH_NDIVSM      ,  /* _INT_SCALAR     */
  AB7_INVARS_PLOWAN_LCALC   ,  /* _INT_ARRAY      */
  AB7_INVARS_IRDBSCOUP      ,  /* _INT_SCALAR     */
  AB7_INVARS_GWLS_MODEL_PARAMETER,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_EPH_MUSTAR     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GENAFM         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NPSPINOR       ,  /* _INT_SCALAR     */
  AB7_INVARS_NGFFTDG        ,  /* _INT_ARRAY      */
  AB7_INVARS_EPH_FSMEAR     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_DMFT_SOLV      ,  /* _INT_SCALAR     */
  AB7_INVARS_PREPGKK        ,  /* _INT_SCALAR     */
  AB7_INVARS_UCRPA          ,  /* _INT_SCALAR     */
  AB7_INVARS_CD_MAX_FREQ    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_RPRIMD_ORIG    ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PRTVXC         ,  /* _INT_SCALAR     */
  AB7_INVARS_NFFTDG         ,  /* _INT_SCALAR     */
  AB7_INVARS_PARAL_RF       ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_TOL_3BT    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_ISTATIMG       ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_TYPFRAG    ,  /* _INT_ARRAY      */
  AB7_INVARS_USEFOCK        ,  /* _INT_SCALAR     */
  AB7_INVARS_NBDBLOCK       ,  /* _INT_SCALAR     */
  AB7_INVARS_IRDWFKFINE     ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_TWEAKS  ,  /* _INT_SCALAR     */
  AB7_INVARS_IRDSCR         ,  /* _INT_SCALAR     */
  AB7_INVARS_ORBMAG         ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_FREQ_MESH   ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_GETGAM_EIG2NKQ ,  /* _INT_SCALAR     */
  AB7_INVARS_GW_FRQRE_INZGRID,  /* _INT_SCALAR     */
  AB7_INVARS_GWLS_FIRST_SEED,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_TOL        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_USE_SLK        ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_RSOFT   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_LEXEXCH        ,  /* _INT_ARRAY      */
  AB7_INVARS_RCUT           ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PRTSUSCEP      ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWSUSHAT      ,  /* _INT_SCALAR     */
  AB7_INVARS_DIEMIXMAG      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NSYM           ,  /* _INT_SCALAR     */
  AB7_INVARS_BERRYSTEP      ,  /* _INT_SCALAR     */
  AB7_INVARS_IRANDOM        ,  /* _INT_SCALAR     */
  AB7_INVARS_RFUSER         ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTPHDOS       ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTCIF         ,  /* _INT_SCALAR     */
  AB7_INVARS_DMATPAWU       ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PLOWAN_BANDF   ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTFSURF       ,  /* _INT_SCALAR     */
  AB7_INVARS_ACCURACY       ,  /* _INT_SCALAR     */
  AB7_INVARS_NPBAND         ,  /* _INT_SCALAR     */
  AB7_INVARS_D3E_PERT2_ELFD ,  /* _INT_SCALAR     */
  AB7_INVARS_GETBSRESO      ,  /* _INT_SCALAR     */
  AB7_INVARS_GWLS_KMAX_NUMERIC,  /* _INT_SCALAR     */
  AB7_INVARS_IXC_SIGMA      ,  /* _INT_SCALAR     */
  AB7_INVARS_NFFT           ,  /* _INT_SCALAR     */
  AB7_INVARS_RFASR          ,  /* _INT_SCALAR     */
  AB7_INVARS_QPTRLATT       ,  /* _INT_ARRAY      */
  AB7_INVARS_FFTGW          ,  /* _INT_SCALAR     */
  AB7_INVARS_NEB_SPRING     ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PLOWAN_NT      ,  /* _INT_SCALAR     */
  AB7_INVARS_RFELFD         ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTVHXC        ,  /* _INT_SCALAR     */
  AB7_INVARS_MBPT_SCISS     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_MGFFT          ,  /* _INT_SCALAR     */
  AB7_INVARS_DMATUDIAG      ,  /* _INT_SCALAR     */
  AB7_INVARS_RECGRATIO      ,  /* _INT_SCALAR     */
  AB7_INVARS_MDWALL         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GETPAWDEN      ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_QRATIO  ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PLOWAN_PROJCALC,  /* _INT_ARRAY      */
  AB7_INVARS_NPWKSS         ,  /* _INT_SCALAR     */
  AB7_INVARS_FRICTION       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GWLS_EXCHANGE  ,  /* _INT_SCALAR     */
  AB7_INVARS_GW_NQLWL       ,  /* _INT_SCALAR     */
  AB7_INVARS_GWCALCTYP      ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_HAYD_TERM   ,  /* _INT_SCALAR     */
  AB7_INVARS_IRDSUSCEP      ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_INTERP_METHOD,  /* _INT_SCALAR     */
  AB7_INVARS_ZEEMANFIELD    ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_FOCKTOLDFE     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_CHKPRIM        ,  /* _INT_SCALAR     */
  AB7_INVARS_NQPT           ,  /* _INT_SCALAR     */
  AB7_INVARS_DOSDELTAE      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PRTEFG         ,  /* _INT_SCALAR     */
  AB7_INVARS_IRDVDW         ,  /* _INT_SCALAR     */
  AB7_INVARS_KPTRLEN        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_KPTGW          ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_BS_COUPLING    ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTCTQMC_CHECK,  /* _INT_SCALAR     */
  AB7_INVARS_TOLRDE         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GETVEL         ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTBANDF      ,  /* _INT_SCALAR     */
  AB7_INVARS_BDGW           ,  /* _INT_ARRAY      */
  AB7_INVARS_DMFT_READ_OCCND,  /* _INT_SCALAR     */
  AB7_INVARS_PAWXCDEV       ,  /* _INT_SCALAR     */
  AB7_INVARS_ESMEAR         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PH_NQSHIFT     ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTCHECK      ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTWF          ,  /* _INT_SCALAR     */
  AB7_INVARS_KPTBOUNDS      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_MEM_TEST       ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_GCUT    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_IRDDDB         ,  /* _INT_SCALAR     */
  AB7_INVARS_KPTNS_HF       ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_GWCOMP         ,  /* _INT_SCALAR     */
  AB7_INVARS_WVL_CRMULT     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_MQGRIDDG       ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_NSTATES     ,  /* _INT_SCALAR     */
  AB7_INVARS_NDYNIMAGE      ,  /* _INT_SCALAR     */
  AB7_INVARS_AMU_ORIG       ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_GWLS_NSEEDS    ,  /* _INT_SCALAR     */
  AB7_INVARS_QPTN           ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_TL_NPRCCG      ,  /* _INT_SCALAR     */
  AB7_INVARS_BFIELD         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PAWPRTWF       ,  /* _INT_SCALAR     */
  AB7_INVARS_STRFACT        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_USERRB         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_CHNEUT         ,  /* _INT_SCALAR     */
  AB7_INVARS_SYMDYNMAT      ,  /* _INT_SCALAR     */
  AB7_INVARS_IRD1WF         ,  /* _INT_SCALAR     */
  AB7_INVARS_GWLS_PRINT_DEBUG,  /* _INT_SCALAR     */
  AB7_INVARS_NELECT         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_VDW_DF_DSOFT   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_POSITRON       ,  /* _INT_SCALAR     */
  AB7_INVARS_DIELNG         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NGKPT          ,  /* _INT_ARRAY      */
  AB7_INVARS_CHKDILATMX     ,  /* _INT_SCALAR     */
  AB7_INVARS_JPAWU          ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_OMEGASIMAX     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PAWLMIX        ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTLDEN        ,  /* _INT_SCALAR     */
  AB7_INVARS_OPTSTRESS      ,  /* _INT_SCALAR     */
  AB7_INVARS_WTQ            ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_RFDDK          ,  /* _INT_SCALAR     */
  AB7_INVARS_GWLS_STERN_KMAX,  /* _INT_SCALAR     */
  AB7_INVARS_PRTATLIST      ,  /* _INT_ARRAY      */
  AB7_INVARS_PRTPHSURF      ,  /* _INT_SCALAR     */
  AB7_INVARS_GA_RULES       ,  /* _INT_ARRAY      */
  AB7_INVARS_WTATCON        ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_USEPAW         ,  /* _INT_SCALAR     */
  AB7_INVARS_SMDELTA        ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWUJRAD       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_BANDPP         ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTCTQMC_BASIS,  /* _INT_SCALAR     */
  AB7_INVARS_LDAMINUSHALF   ,  /* _INT_ARRAY      */
  AB7_INVARS_GETWFK         ,  /* _INT_SCALAR     */
  AB7_INVARS_DIISMEMORY     ,  /* _INT_SCALAR     */
  AB7_INVARS_GETWFQ         ,  /* _INT_SCALAR     */
  AB7_INVARS_FBAND          ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_MAGCON_LAMBDA  ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NATSPH_EXTRA   ,  /* _INT_SCALAR     */
  AB7_INVARS_NBANDKSS       ,  /* _INT_SCALAR     */
  AB7_INVARS_FREQSPMAX      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_TL_RADIUS      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GW_SIGXCORE    ,  /* _INT_SCALAR     */
  AB7_INVARS_PITRANSFORM    ,  /* _INT_SCALAR     */
  AB7_INVARS_DIEGAP         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_ZIONTYPAT      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PAWMIXDG       ,  /* _INT_SCALAR     */
  AB7_INVARS_STRPRECON      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PRTDIPOLE      ,  /* _INT_SCALAR     */
  AB7_INVARS_POL            ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_DMFTCTQMC_MOV  ,  /* _INT_SCALAR     */
  AB7_INVARS_ISTATSHFT      ,  /* _INT_SCALAR     */
  AB7_INVARS_NPWEPS         ,  /* _INT_SCALAR     */
  AB7_INVARS_SPNORBSCL      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GA_ALGOR       ,  /* _INT_SCALAR     */
  AB7_INVARS_WTK            ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_GA_N_RULES     ,  /* _INT_SCALAR     */
  AB7_INVARS_ECUTWFN        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GWLS_RECYCLE   ,  /* _INT_SCALAR     */
  AB7_INVARS_AWTR           ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_ACUTMIN ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NKPTGW         ,  /* _INT_SCALAR     */
  AB7_INVARS_OCC_ORIG       ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NPULAYIT       ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFT_NWLI      ,  /* _INT_SCALAR     */
  AB7_INVARS_RATSPH_EXTRA   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_INCLVKB        ,  /* _INT_SCALAR     */
  AB7_INVARS_NNOS           ,  /* _INT_SCALAR     */
  AB7_INVARS_EFMAS_DIM      ,  /* _INT_SCALAR     */
  AB7_INVARS_GWRPACORR      ,  /* _INT_SCALAR     */
  AB7_INVARS_DILATMX        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_DMATPUOPT      ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWECUTDG      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GWLS_KMAX_COMPLEMENT,  /* _INT_SCALAR     */
  AB7_INVARS_TOLIMG         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_WVL_BIGDFT_COMP,  /* _INT_SCALAR     */
  AB7_INVARS_PRTSTM         ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_SUPERCELL  ,  /* _INT_ARRAY      */
  AB7_INVARS_ACELL_ORIG     ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_DMFTCTQMC_ORDER,  /* _INT_SCALAR     */
  AB7_INVARS_DIEMAC         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_EFMAS_DEG_TOL  ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_SIGNPERM       ,  /* _INT_SCALAR     */
  AB7_INVARS_BRAVAIS        ,  /* _INT_ARRAY      */
  AB7_INVARS_EFMAS          ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_PHISOFT ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PRTGDEN        ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTFC          ,  /* _INT_SCALAR     */
  AB7_INVARS_BERRYOPT       ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTVDW         ,  /* _INT_SCALAR     */
  AB7_INVARS_IXCPOSITRON    ,  /* _INT_SCALAR     */
  AB7_INVARS_GETQPS         ,  /* _INT_SCALAR     */
  AB7_INVARS_NDIVSM         ,  /* _INT_SCALAR     */
  AB7_INVARS_TIMOPT         ,  /* _INT_SCALAR     */
  AB7_INVARS_GWLS_LIST_PROJ_FREQ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_RFDIR          ,  /* _INT_ARRAY      */
  AB7_INVARS_SHIFTK         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_VDW_DF_THRESHOLD,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_KBERRY         ,  /* _INT_ARRAY      */
  AB7_INVARS_KPTOPT         ,  /* _INT_SCALAR     */
  AB7_INVARS_NSPDEN         ,  /* _INT_SCALAR     */
  AB7_INVARS_D3E_PERT1_ATPOL,  /* _INT_ARRAY      */
  AB7_INVARS_BS_LOBAND      ,  /* _INT_ARRAY      */
  AB7_INVARS_USEWVL         ,  /* _INT_SCALAR     */
  AB7_INVARS_GETOCC         ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTPOSCAR      ,  /* _INT_SCALAR     */
  AB7_INVARS_PH_INTMETH     ,  /* _INT_SCALAR     */
  AB7_INVARS_GWLS_CORRELATION,  /* _INT_SCALAR     */
  AB7_INVARS_NOMEGASRD      ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_NQPTS   ,  /* _INT_SCALAR     */
  AB7_INVARS_IBOXCUT        ,  /* _INT_SCALAR     */
  AB7_INVARS_PARAL_ATOM     ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWNTHETA      ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFT_T2G       ,  /* _INT_SCALAR     */
  AB7_INVARS_USEDMFT        ,  /* _INT_SCALAR     */
  AB7_INVARS_UCRPA_WINDOW   ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NSHIFTK        ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTELF         ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_RCUT    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GPU_LINALG_LIMIT,  /* _INT_SCALAR     */
  AB7_INVARS_PAWFATBND      ,  /* _INT_SCALAR     */
  AB7_INVARS_OMEGASRDMAX    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NATRD          ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFT_TOLFREQ   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NOSEINERT      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_WVL_NGAUSS     ,  /* _INT_ARRAY      */
  AB7_INVARS_INTXC          ,  /* _INT_SCALAR     */
  AB7_INVARS_USEXCNHAT_ORIG ,  /* _INT_SCALAR     */
  AB7_INVARS_IXCROT         ,  /* _INT_SCALAR     */
  AB7_INVARS_GETDELFD       ,  /* _INT_SCALAR     */
  AB7_INVARS_SPGAXOR        ,  /* _INT_SCALAR     */
  AB7_INVARS_DFPT_SCISS     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GET1DEN        ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWSPNORB      ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTBANDI      ,  /* _INT_SCALAR     */
  AB7_INVARS_NNSCLO         ,  /* _INT_SCALAR     */
  AB7_INVARS_RFMAGN         ,  /* _INT_SCALAR     */
  AB7_INVARS_IRDQPS         ,  /* _INT_SCALAR     */
  AB7_INVARS_NOMEGASI       ,  /* _INT_SCALAR     */
  AB7_INVARS_OPTNLXCCC      ,  /* _INT_SCALAR     */
  AB7_INVARS_DIECUT         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_EPH_FSEWIN     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_JELLSLAB       ,  /* _INT_SCALAR     */
  AB7_INVARS_PPMODEL        ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWOPTMIX      ,  /* _INT_SCALAR     */
  AB7_INVARS_OPTCELL        ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWPRTDOS      ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTGEO         ,  /* _INT_SCALAR     */
  AB7_INVARS_MDF_EPSINF     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_CD_SUBSET_FREQ ,  /* _INT_ARRAY      */
  AB7_INVARS_BERRYSAV       ,  /* _INT_SCALAR     */
  AB7_INVARS_TOLWFR         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PH_SMEAR       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PRTVOLIMG      ,  /* _INT_SCALAR     */
  AB7_INVARS_LOTF_CLASSIC   ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWNHATXC      ,  /* _INT_SCALAR     */
  AB7_INVARS_NATVSHIFT      ,  /* _INT_SCALAR     */
  AB7_INVARS_IONMOV         ,  /* _INT_SCALAR     */
  AB7_INVARS_STRTARGET      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_VDW_DF_DRATIO  ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_LOTF_NITEX     ,  /* _INT_SCALAR     */
  AB7_INVARS_W90INIPRJ      ,  /* _INT_SCALAR     */
  AB7_INVARS_NSPPOL         ,  /* _INT_SCALAR     */
  AB7_INVARS_NFREQRE        ,  /* _INT_SCALAR     */
  AB7_INVARS_SUPERCELL      ,  /* _INT_ARRAY      */
  AB7_INVARS_BS_COULOMB_TERM,  /* _INT_SCALAR     */
  AB7_INVARS_NBANDHF        ,  /* _INT_SCALAR     */
  AB7_INVARS_PTGROUPMA      ,  /* _INT_SCALAR     */
  AB7_INVARS_HYB_RANGE_DFT  ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_TPHYSEL        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_SYMREL         ,  /* _INT_ARRAY      */
  AB7_INVARS_PAWCPXOCC      ,  /* _INT_SCALAR     */
  AB7_INVARS_WVL_FRMULT     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_BS_INTERP_PREP ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWUSECP       ,  /* _INT_SCALAR     */
  AB7_INVARS_USERIB         ,  /* _INT_SCALAR     */
  AB7_INVARS_RF2_DKDK       ,  /* _INT_SCALAR     */
  AB7_INVARS_SO_PSP         ,  /* _INT_ARRAY      */
  AB7_INVARS_PRTPOT         ,  /* _INT_SCALAR     */
  AB7_INVARS_PLOWAN_REALSPACE,  /* _INT_SCALAR     */
  AB7_INVARS_RF2_DKDE       ,  /* _INT_SCALAR     */
  AB7_INVARS_DIPDIP         ,  /* _INT_SCALAR     */
  AB7_INVARS_OPTFORCES      ,  /* _INT_SCALAR     */
  AB7_INVARS_BRVLTT         ,  /* _INT_SCALAR     */
  AB7_INVARS_D3E_PERT1_PHON ,  /* _INT_SCALAR     */
  AB7_INVARS_BOXCUTMIN      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NTIME          ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTVHA         ,  /* _INT_SCALAR     */
  AB7_INVARS_NPSP           ,  /* _INT_SCALAR     */
  AB7_INVARS_SPBROAD        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GOPRECPRM      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PRTVOL         ,  /* _INT_SCALAR     */
  AB7_INVARS_NQPTDM         ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTBLTZTRP     ,  /* _INT_SCALAR     */
  AB7_INVARS_DENSFOR_PRED   ,  /* _INT_SCALAR     */
  AB7_INVARS_ISECUR         ,  /* _INT_SCALAR     */
  AB7_INVARS_RHOQPMIX       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NATSPH         ,  /* _INT_SCALAR     */
  AB7_INVARS_GETBSCOUP      ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFT_NLAMBDA   ,  /* _INT_SCALAR     */
  AB7_INVARS_PIMD_CONSTRAINT,  /* _INT_SCALAR     */
  AB7_INVARS_NWFSHIST       ,  /* _INT_SCALAR     */
  AB7_INVARS_FREQSPMIN      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_ECUT           ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PAWPRT_K       ,  /* _INT_SCALAR     */
  AB7_INVARS_PAWOPTOSC      ,  /* _INT_SCALAR     */
  AB7_INVARS_NKPT           ,  /* _INT_SCALAR     */
  AB7_INVARS_SYMAFM         ,  /* _INT_ARRAY      */
  AB7_INVARS_PAWPRT_B       ,  /* _INT_SCALAR     */
  AB7_INVARS_NCONEQ         ,  /* _INT_SCALAR     */
  AB7_INVARS_RPRIM_ORIG     ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_EFMAS_BANDS    ,  /* _INT_ARRAY      */
  AB7_INVARS_VDW_XC         ,  /* _INT_SCALAR     */
  AB7_INVARS_IEIG2RF        ,  /* _INT_SCALAR     */
  AB7_INVARS_USE_NONSCF_GKK ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_INTERP_MODE ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_EH_CUTOFF   ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_GETSUSCEP      ,  /* _INT_SCALAR     */
  AB7_INVARS_WFK_TASK       ,  /* _INT_SCALAR     */
  AB7_INVARS_GWLS_N_PROJ_FREQ,  /* _INT_SCALAR     */
  AB7_INVARS_NATOM          ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_NDPTS   ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTDOS         ,  /* _INT_SCALAR     */
  AB7_INVARS_SPGORIG        ,  /* _INT_SCALAR     */
  AB7_INVARS_NUCDIPMOM      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_D3E_PERT1_DIR  ,  /* _INT_ARRAY      */
  AB7_INVARS_EPH_INTMETH    ,  /* _INT_SCALAR     */
  AB7_INVARS_SPINAT         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_XREDSPH_EXTRA  ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NNSCLOHF       ,  /* _INT_SCALAR     */
  AB7_INVARS_SYMMORPHI      ,  /* _INT_SCALAR     */
  AB7_INVARS_CD_CUSTOMNIMFRQS,  /* _INT_SCALAR     */
  AB7_INVARS_FERMIE_NEST    ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_BS_HAYDOCK_TOL ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NPSPALCH       ,  /* _INT_SCALAR     */
  AB7_INVARS_AUXC_IXC       ,  /* _INT_SCALAR     */
  AB7_INVARS_D3E_PERT2_ATPOL,  /* _INT_ARRAY      */
  AB7_INVARS_ZNUCL          ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PAWNZLM        ,  /* _INT_SCALAR     */
  AB7_INVARS_NFREQIM        ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFT_NWLO      ,  /* _INT_SCALAR     */
  AB7_INVARS_UPAWU          ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PRTWANT        ,  /* _INT_SCALAR     */
  AB7_INVARS_MBAND          ,  /* _INT_SCALAR     */
  AB7_INVARS_GETSCR         ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTKDEN        ,  /* _INT_SCALAR     */
  AB7_INVARS_OPTDRIVER      ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_DAMAX   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NSHIFTK_ORIG   ,  /* _INT_SCALAR     */
  AB7_INVARS_GW_QPRANGE     ,  /* _INT_SCALAR     */
  AB7_INVARS_QPTOPT         ,  /* _INT_SCALAR     */
  AB7_INVARS_EPH_FERMIE     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PRTBBB         ,  /* _INT_SCALAR     */
  AB7_INVARS_PH_QSHIFT      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_RED_DFIELD     ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_NKPTHF         ,  /* _INT_SCALAR     */
  AB7_INVARS_DTION          ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PH_QPATH       ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_VDW_DF_ARATIO  ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PRTEIG         ,  /* _INT_SCALAR     */
  AB7_INVARS_VACNUM         ,  /* _INT_SCALAR     */
  AB7_INVARS_D3E_PERT3_DIR  ,  /* _INT_ARRAY      */
  AB7_INVARS_GOPRECON       ,  /* _INT_SCALAR     */
  AB7_INVARS_ISTWFK         ,  /* _INT_ARRAY      */
  AB7_INVARS_MIXALCH_ORIG   ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PRTGKK         ,  /* _INT_SCALAR     */
  AB7_INVARS_QUADMOM        ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_SHIFTK_ORIG    ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_PAWPRTVOL      ,  /* _INT_SCALAR     */
  AB7_INVARS_TOLVRS         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_MAGCONON       ,  /* _INT_SCALAR     */
  AB7_INVARS_F6OF2_SLA      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_WFMIX          ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NEB_ALGO       ,  /* _INT_SCALAR     */
  AB7_INVARS_POSDOPPLER     ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_INTERP_M3_WIDTH,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_RECNREC        ,  /* _INT_SCALAR     */
  AB7_INVARS_TSMEAR         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_KPTRLATT_ORIG  ,  /* _INT_ARRAY      */
  AB7_INVARS_CD_FULL_GRID   ,  /* _INT_SCALAR     */
  AB7_INVARS_BS_CALCTYPE    ,  /* _INT_SCALAR     */
  AB7_INVARS_POSNSTEP       ,  /* _INT_SCALAR     */
  AB7_INVARS_IRDHAYDOCK     ,  /* _INT_SCALAR     */
  AB7_INVARS_FXCARTFACTOR   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_RECNPATH       ,  /* _INT_SCALAR     */
  AB7_INVARS_ELPH2_IMAGDEN  ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_NKPATH         ,  /* _INT_SCALAR     */
  AB7_INVARS_GW_FRQRE_TANGRID,  /* _INT_SCALAR     */
  AB7_INVARS_PLOWAN_NBL     ,  /* _INT_ARRAY      */
  AB7_INVARS_ESHIFT         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_ADPIMD_GAMMA   ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_ISTATR         ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_ZAB     ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PAWNPHI        ,  /* _INT_SCALAR     */
  AB7_INVARS_D3E_PERT1_ELFD ,  /* _INT_SCALAR     */
  AB7_INVARS_LPAWU          ,  /* _INT_ARRAY      */
  AB7_INVARS_XCLEVEL        ,  /* _INT_SCALAR     */
  AB7_INVARS_FOCKOPTMIX     ,  /* _INT_SCALAR     */
  AB7_INVARS_NTYPPURE       ,  /* _INT_SCALAR     */
  AB7_INVARS_MQGRID         ,  /* _INT_SCALAR     */
  AB7_INVARS_DDB_SHIFTQ     ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_UCRPA_BANDS    ,  /* _INT_ARRAY      */
  AB7_INVARS_PLOWAN_NATOM   ,  /* _INT_SCALAR     */
  AB7_INVARS_JFIELDDIR      ,  /* _INT_ARRAY      */
  AB7_INVARS_TOLRFF         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PH_WSTEP       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_SYMSIGMA       ,  /* _INT_SCALAR     */
  AB7_INVARS_POSOCC         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_SLABZEND       ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PIMASS         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_DMFT_ITER      ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFT_DC        ,  /* _INT_SCALAR     */
  AB7_INVARS_F4OF2_SLA      ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_VEL_CELL_ORIG  ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_HYB_MIXING_SR  ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_PPMFRQ         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_TMESH          ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_CHKEXIT        ,  /* _INT_SCALAR     */
  AB7_INVARS_PRTDENSPH      ,  /* _INT_SCALAR     */
  AB7_INVARS_VDW_DF_NRPTS   ,  /* _INT_SCALAR     */
  AB7_INVARS_TOLDFF         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_TOLDFE         ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_ORTALG         ,  /* _INT_SCALAR     */
  AB7_INVARS_VPRTRB         ,  /* _DOUBLE_ARRAY   */
  AB7_INVARS_MKQMEM         ,  /* _INT_SCALAR     */
  AB7_INVARS_GETDKDE        ,  /* _INT_SCALAR     */
  AB7_INVARS_STMBIAS        ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_XC_DENPOS      ,  /* _DOUBLE_SCALAR  */
  AB7_INVARS_GETDKDK        ,  /* _INT_SCALAR     */
  AB7_INVARS_DMFTCTQMC_TRIQS_NLEG,  /* _INT_SCALAR     */
  AB7_INVARS_NGFFT          ,  /* _INT_ARRAY      */
  AB7_INVARS_N_IDS
} Ab7InvarsIds;



Ab7InvarsTypes ab7_invars_get_type_from_id(Ab7InvarsIds id);
/**
 * AB7_INVARS_TYPE:
 * @A: an #Ab7InvarsIds id.
 *
 * Get the type of a given attribute of Dtset structure.
 *
 * Returns: a #Ab7InvarsTypes id.
 */
#define AB7_INVARS_TYPE(A) ab7_invars_get_type_from_id(A)
/**
 * AB7_INVARS_STR:
 * @A: an #Ab7InvarsIds id.
 *
 * Get a string corresponding to the attribute name.
 *
 * Returns: a string owned by ABINIT, do not free or modify it.
 */
#define AB7_INVARS_STR(A) #A

/**
 * Ab7Invars:
 *
 * An object to handle an array of ABINIT datasets, read from a file.
 */
typedef int Ab7Invars;

/**
 * ab7_invars_new_from_file:
 * @filename: a string, NULL terminated.
 *
 * Parse the given file using ABINIT routines and allocate a
 * dtsets array. This array must be deallocated after use with
 * ab7_invars_free().
 *
 * Returns: an #Ab7Invars object or NULL on failure.
 */
Ab7Invars* ab7_invars_new_from_file(const char *filename);
/**
 * ab7_invars_new_from_file_with_pseudo:
 * @filename: a string, NULL terminated.
 * @pspfiles: an array of strings, NULL terminated. Can be NULL.
 *
 * Parse the given file using ABINIT routines and allocate a
 * dtsets array. This array must be deallocated after use with
 * ab7_invars_free(). If pseudo files are provided with @pspfiles,
 * some further initialisations of dtset are permitted.
 *
 * Returns: an #Ab7Invars object or NULL on failure.
 */
Ab7Invars* ab7_invars_new_from_file_with_pseudo(const char *filename, const char **pspfiles);
/**
 * ab7_invars_new_from_string:
 * @string: a string, NULL terminated.
 *
 * Parse the given string using ABINIT routines and allocate a
 * dtsets array. This array must be deallocated after use with
 * ab7_invars_free().
 *
 * Returns: an #Ab7Invars object or NULL on failure.
 */
Ab7Invars* ab7_invars_new_from_string(const char *string);
/**
 * ab7_invars_free:
 * @ptr: the dataset array to handle.
 *
 * Clean all allocated memory from the data set allocation.
 */
void ab7_invars_free(Ab7Invars *ptr);

/**
 * ab7_invars_get_ndtset:
 * @ptr: the dataset array to handle.
 * @ndtset: a location to store the returned value.
 *
 * An array of datasets may contain more than one. Test it with this
 * routine. @ndtset will contains the number of allocated datasets (in
 * addition to the default one).
 *
 * Returns: #AB7_NO_ERROR if @ptr is valid and correctly parsed.
 */
Ab7Error ab7_invars_get_ndtset(Ab7Invars *ptr, int *ndtset);
/**
 * ab7_invars_get_integer:
 * @ptr: the dataset array to handle.
 * @id: an attribute id, see dtset_c.h.
 * @idtset: the number of the dtset to read, 0 is default value.
 * @value: a location to store the returned value.
 *
 * Use this method to get the value of an integer attribute. @idtset
 * must be in [0;n] where n is the returned value of
 * ab7_invars_get_ndtset(). If @id is unknown, return value is
 * 0. For real attributes, see ab7_invars_get_real().
 *
 * Returns: #AB7_NO_ERROR if values are correctly read.
 */
Ab7Error ab7_invars_get_integer(Ab7Invars *ptr, Ab7InvarsIds id,
                                int idtset, int *value);
/**
 * ab7_invars_get_real:
 * @ptr: the dataset array to handle.
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 * @value: a location to store the returned value.
 *
 * Use this method to get the value of a double attribute. @idtset
 * must be in [0;n] where n is the return value of
 * ab7_invars_get_ndtset(). If @id is unknown, return value is
 * undefined. For integer attributes, see ab7_invars_get_integer().
 *
 * Returns: #AB7_NO_ERROR if values are correctly read.
 */
Ab7Error ab7_invars_get_real(Ab7Invars *ptr, Ab7InvarsIds id,
                             int idtset, double *value);

/**
 * ab7_invars_get_shape:
 * @ptr: the dataset array to handle.
 * @n: a location to store the number of dimensions.
 * @dims: an array with 7 integers ;
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 *
 * This method is used to poll the size of an array attribute. The
 * shape of the attribute is stored in @dims. Only the @n first values
 * of @dims are relevant.
 *
 * Returns: #AB7_NO_ERROR if values are correctly read.
 */
Ab7Error ab7_invars_get_shape(Ab7Invars *ptr, int *n, int dims[7],
			      Ab7InvarsIds id, int idtset);
/**
 * ab7_invars_get_integer_array:
 * @ptr: the dataset array to handle.
 * @values: an allocated array of @n values ;
 * @n: the size of the given array ;
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 *
 * This method is used to read the values of an array. The array must
 * already be allocated. To know its size, use ab7_invars_get_shape().
 *
 * Returns: #AB7_NO_ERROR if values are correctly read.
 */
Ab7Error ab7_invars_get_integer_array(Ab7Invars *ptr, int *values, size_t n,
				      Ab7InvarsIds id, int idtset);
/**
 * ab7_invars_get_real_array:
 * @ptr: the dataset array to handle.
 * @values: an allocated array of @n values ;
 * @n: the size of the given array ;
 * @id: an attribute id, see dtset_c.h ;
 * @idtset: the number of the dtset to read, 0 is default.
 *
 * This method is used to read the values of an array. The array must
 * already be allocated. To know its size, use ab7_invars_get_shape().
 *
 * Returns: #AB7_NO_ERROR if values are correctly read.
 */
Ab7Error ab7_invars_get_real_array(Ab7Invars *ptr, double *values, size_t n,
				   Ab7InvarsIds id, int idtset);

#endif
