################################################################################
#
# For Odyssey
#
################################################################################
#PJM --name "S1380_CPL"
#PJM -L rscgrp=regular-o
#PJM -L node=2
#PJM -g ga41
#PJM --rsc-list "elapse=01:00:00"
#PJM --mpi "proc=5"
#PJM --mpi "rank-map-bynode"
#PJM -S

##################################################

SolarConst=1380
exp_name=S${SolarConst}_CPL
TOPDIR=/work/gz06/n49000/Dennou-CCM2/Dennou-CCM/
EXPDIR=${TOPDIR}/exp/APEI07Couple/couple_AlbMod
RUNDIR=${EXPDIR}/S${SolarConst}
MPIRUN=mpiexec

###################################################
# coupled model configuration

DCCMConfPath=${EXPDIR}/../common/DCCM_ATM_T21-OCN_Pl42.conf

###################################################
# Atmosphere model configuration
atm_pe=${TOPDIR}/bin/atm_driver
atm_wdir=${RUNDIR}/atm

DCPAM_BIN_DIR=/work/gz06/n49000/Dennou-CCM2/dcpam5_v2015-08-04_ykawai_ext/src/main
atm_init_data_pe=${DCPAM_BIN_DIR}/dcpam_init_data
atm_init_data_sfc_pe=${DCPAM_BIN_DIR}/dcpam_init_data_surface

atm_nml_template=${EXPDIR}/../common/dcpam/dcpam_APEI07Couple_T21L16.conf
atm_init_data_nml=${EXPDIR}/../common/dcpam/dcpam_init_data_APEI07_T21L16.conf
atm_init_data_sfc_nml=${EXPDIR}/../common/dcpam/dcpam_surface_data_E_T21_280K.conf

atm_PE_NUM=4
atm_THREADS_NUM=12

####################################################
# Ocean model configuration
ocn_pe=${TOPDIR}/bin/ocn_driver
ocn_wdir=${RUNDIR}/ocn
ocn_nml_template=${EXPDIR}/../common/dogcm/dogcm_APEI07Couple_Pl64L60_I07SfcAlbMod.conf
ocn_PE_NUM=1
ocn_THREADS_NUM=12
#ocn_nodefile="ocn_nodefile"

DCPOM_BIN_DIR=/work/gz06/n49000/Dennou-CCM2/lib/DCPOM/bin/
ocn_standalone_pe=${DCPOM_BIN_DIR}/dogcm_axisym
ocn_standalone_libdir=${DCPOM_BIN_DIR}/../lib/

ocn_standalone_PE_NUM=1
ocn_standalone_THREADS_NUM=12

##################################

nCycle=1
StartCycleNum=2
coupledTimeIntrvPerCycle=$((2*365))
standaloneTimeIntrvPerCycle=$((50*365)) 
coupleODelTimeHour=4
standaloneODelTimeHour=12 
HistIntValueDayCPLRun=146
#HistIntValueDaySTDAloneRun=365
EXP_NAME=S1380

FlagModAlbedoBasedOnTempSGS=true
coupledRunSkipSCyc=false
#---------

PBS_O_WORKDIR=${EXPDIR}

module load netcdf
module load netcdf-fortran
module load hdf5
module load pnetcdf

#export XOS_MMM_L_ARENA_FREE=1
export FORT90L="-Wl,-T"
#export OMPI_MCA_plm_ple_memory_allocation_policy=bind_local
export PLE_MPI_STD_EMPTYFILE="off"
export PARALLEL=${atm_THREADS_NUM}
export OMP_NUM_THREADS=${atm_THREADS_NUM}
#export fu11bf=1

## End of setting *******************************************************************************

source ${EXPDIR}/../common/APEI07CoupleExp_job_inc_Odyssey.sh
