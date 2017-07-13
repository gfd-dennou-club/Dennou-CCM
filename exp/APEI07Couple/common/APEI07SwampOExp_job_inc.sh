#------------------------------------------------------------------------------------------
# Copyright (c) 2016-2017 Yuta Kawai. All rights reserved.
#-------------------------------------------------------------------------------------------
# * Dennou-CCM launcher script
#   A shell script used to perform climate numerical experiments with a coupled model. 
#
#   The coupled system is composed of atmopheric general circulation model, ocean general
#   circulation model, sea ice model and others. This script assumes that coupling 'Jcup'
#   is used to couple their models. 
#   
#   The coupled system is integrated temporally in the  following steps. 
#   1) First, the coupled model  is integrated for a short period (ex. 6 months). After that,
#      Some ruby scripts generate data files necessary to specify initial and boundary
#      condition in next OGCM standalone run. 
#   2) Ocean model (often run with sea ice model) alone is integrated for much longer period
#      than that of 1) (ex. 10 years). After that, some ruby scripts generate data files necessary
#      to specify initial condition of ocean in next coupled run.
#   3) Go back to step 1. 
#   
#********************************************************************************************


${coupledRunSkipSCyc:=false}
${FlagVerticalFilter:=false}
${HDEFoldTimeHour:=3}
${KMINGSVF:=2}
${KMAXGSVF:=31}
${FlagPRCPPC:=true}

#--------------------------------------------------------------------------------------------

## Definition of some functions ##############################

function create_dir() {
    dirPath=$1
    if [ ! -e $dirPath ]; then
	echo "Create directory '${dirPath}' .."
	mkdir $dirPath
#	chown ykawai:ykawai $dirPath
#	chmod g+w $dirPath
    else
	echo "Skip mkdir operation because '${dirPath}' already exist."
    fi
}

### Main parts ##############################################

# Prepare directories to save output data.

echo "Create some directories to save data.."
for ((n=StartCycleNum; n<=nCycle; n++)) ; do
    create_dir "${atm_wdir}/cycle${n}-couple"
done

cd $PBS_O_WORKDIR

#- Perform temporal integration of coupled system -------------------------------

coupledRunRestartTime=$(((StartCycleNum-1)*coupledTimeIntrvPerCycle))
for ((n=StartCycleNum; n<=nCycle; n++)) ; do

    ######################################################################
    # Run coupled model
    ######################################################################
    
    atmDirPath="${atm_wdir}/cycle${n}-couple"
    coupledRunEndTime=$((coupledRunRestartTime + coupledTimeIntrvPerCycle))
    
    echo "-- cycle=${n} (AGCM run) -- time range =${coupledRunRestartTime} - ${coupledRunEndTime} [day]"
    
    echo "** Create configuration file for AGCM **"

    sedArgs=`cat <<EOF 
     s!#restart_file_io_nml_InputFile#!${atm_wdir}/cycle$((n-1))-couple/rst.nc!g;
     s!#restart_file_io_nml_IntValue#!${RestartIntValue}.0!g;
     s!#restart_surftemp_io_nml_InputFile#!${atm_wdir}/cycle$((n-1))-couple/rst_sst.nc!g;
     s!#timeset_nml_RestartTimeValue#!${coupledRunRestartTime}!g;
     s!#timeset_nml_InitYear#!2000!g; 
     s!#timeset_nml_EndYear#!2000!g;
     s!#timeset_nml_EndDay#!$((coupledRunEndTime+1))!g;
     s!#timeset_nml_DelTimeMin#!${DelTimeMin}.0!g;
     s!#gtool_historyauto_nml_IntValue#!${HistIntValueDay}.0!g; 
     s!#rad_DennouAGCM_nml_RstInputFile#!${atm_wdir}/cycle$((n-1))-couple/rst_rad.nc!g;
     s!#rad_DennouAGCM_nml_SolarConst#!${SolarConst}.0!g;
     s!#dynamics_hspl_vas83_nml_HDEFoldTimeValue#!${HDEFoldTimeHour}!g;
     s!#dynamics_hspl_vas83_nml_FlagVertFilter#!${FlagVerticalFilter}!g;
     s!#dynamics_hspl_vas83_nml_KMINGSVF#!${KMINGSVF}!g;
     s!#dynamics_hspl_vas83_nml_KMAXGSVF#!${KMAXGSVF}!g;
     s!#cloud_none_nml_FlagPRCPPC#!${FlagPRCPPC}!g;
EOF
    ` 
    atm_nml=${atmDirPath}/${atm_nml_template##*/}
    sed -e "${sedArgs}" ${atm_nml_template} > ${atm_nml}


    #
    if [ $n -eq $StartCycleNum ] && $coupledRunSkipSCyc ; then
	echo "skip coupled run .."
    else
        echo "** Execute AGCM  ******************************"

	${MPIRUN}                                                   \
	-wdir ${atmDirPath} -env OMP_NUM_THREADS ${atm_THREADS_NUM} \
	-n ${atm_PE_NUM} ${atm_pe} -N=${atm_nml}                    \
	1> Stdout_${exp_name} 2>Stderr_${exp_name}
    fi

    coupledRunRestartTime=${coupledRunEndTime}
done

