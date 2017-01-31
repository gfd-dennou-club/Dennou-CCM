#------------------------------------------------------------------------------------------
# Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
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

cp    ${TOPDIR}/bin/atm_driver ${atm_wdir}
cp    ${TOPDIR}/bin/ocn_driver ${ocn_wdir}
cp    ${ocn_standalone_pedir}/${ocn_standalone_pename} ${ocn_wdir}/ocn_standalone
cp -r ${ocn_standalone_libdir} ${ocn_wdir}

echo "Create some directories to save data.."
for ((n=1; n<=nCycle; n++)) ; do
    create_dir "${atm_wdir}/cycle${n}-couple"
    create_dir "${ocn_wdir}/cycle${n}-couple"
    create_dir "${ocn_wdir}/cycle${n}-standalone"
done

cd $PBS_O_WORKDIR

#- Perform temporal integration of coupled system -------------------------------

coupledRunRestartTime=$(((StartCycleNum-1)*coupledTimeIntrvPerCycle))
for ((n=StartCycleNum; n<=nCycle; n++)) ; do

    ######################################################################
    # Run coupled model
    ######################################################################
    
    atmDirPath="${atm_wdir}/cycle${n}-couple"
    ocnDirPath="${ocn_wdir}/cycle${n}-couple"
    ocnDirPath_standalone="${ocn_wdir}/cycle${n}-standalone"

    coupledRunEndTime=$((coupledRunRestartTime + coupledTimeIntrvPerCycle))
    
    echo "-- cycle=${n} (coupled AOGCM run) -- time range =${coupledRunRestartTime} - ${coupledRunEndTime} [day]"
    
    echo "** Create configuration file for AGCM **"

    sedArgs=`cat <<EOF 
     s!#restart_file_io_nml_InputFile#!${atm_wdir}/cycle$((n-1))-couple/rst.nc!g;
     s!#restart_file_io_nml_IntValue#!730.0!g;
     s!#timeset_nml_RestartTimeValue#!${coupledRunRestartTime}!g;
     s!#timeset_nml_InitYear#!2000!g; 
     s!#timeset_nml_EndYear#!2000!g;
     s!#timeset_nml_EndDay#!$((coupledRunEndTime+1))!g;
     s!#gtool_historyauto_nml_IntValue#!146.0!g; 
     s!#rad_DennouAGCM_nml_RstInputFile#!${atm_wdir}/cycle$((n-1))-couple/rst_rad.nc!g;
     s!#rad_DennouAGCM_nml_SolarConst#!${SolarConst}!g;
EOF
    ` 
    atm_nml=${atmDirPath}/${atm_nml_template##*/}
    sed -e "${sedArgs}" ${atm_nml_template} > ${atm_nml}

    echo "** Create configuration file for OGCM **"
    OcnRestartInFile=""
    SIceRestartInFile=""

    sedArgs=`cat << EOF
     s!#gtool_historyauto_nml_IntValue#!146.0!g;
     s!#gtool_historyauto_nml_OriginValue#!${coupledRunRestartTime}!g;
     s!#gtool_historyauto_nml_TerminusValue#!${coupledRunEndTime}!g;
     s!#OcnRestartFile_nml_InputFileName#!${OcnRestartInFile}!g; 
     s!#OcnRestartFile_nml_OutputFileName#!RestartOcnData.nc!g;
     s!#OcnRestartFile_nml_IntValue#!730.0!g;
     s!#SIceRestartFile_nml_InputFileName#!${SIceRestartInFile}!g; 
     s!#SIceRestartFile_nml_OutputFileName#!RestartSIceData.nc!g;
     s!#SIceRestartFile_nml_IntValue#!730.0!g;
     s!#TemporalInteg_nml_DelTimeHour#!${coupleODelTimeHour}!g;
     s!#TemporalInteg_nml_RestartTimeVal#!${coupledRunRestartTime}!g;
     s!#TemporalInteg_nml_InitYear#!2000!g; 
     s!#TemporalInteg_nml_EndYear#!2000!g; s!#TemporalInteg_nml_EndDay#!$((coupledRunEndTime+1))!g;
     s!#BoundaryCondition_nml_ThermBCSurface#!PrescFlux!g;
     s!#BoundaryCondition_nml_SaltBCSurface#!PrescFlux!g;
     s!#Exp_APECoupleClimate_nml_RunCycle#!${n}!g;
     s!#Exp_APECoupleClimate_nml_RunTypeName#!Coupled!g;
     s!#Exp_APECoupleClimate_nml_SfcBCDataDir#!${ocn_wdir}/cycle$((n-1))-couple/!g;
     s!#Exp_APECoupleClimate_nml_SfcBCMeanInitTime#!${coupledRunRestartTime}.0!g;
     s!#Exp_APECoupleClimate_nml_SfcBCMeanEndTime#!${coupledRunRestartTime}.0!g;
EOF
    `
    sedArgs2=""
    if [ $standaloneTimeIntrvPerCyc -gt 0 ] ; then
	sedArgs2=`cat << EOF
     s!#Exp_APECoupleClimate_nml_RestartDataDir#!${ocn_wdir}/cycle$((n-1))-standalone/!g;
     s!#Exp_APECoupleClimate_nml_RestartMeanInitTime#!$((standaloneTimeIntrvPerCycle))!g;
     s!#Exp_APECoupleClimate_nml_RestartMeanEndTime#!${standaloneTimeIntrvPerCycle}.0!g;
EOF
    `
    else
	sedArgs2=`cat << EOF
     s!#Exp_APECoupleClimate_nml_RestartDataDir#!${ocn_wdir}/cycle$((n-1))-couple/!g;
     s!#Exp_APECoupleClimate_nml_RestartMeanInitTime#!${coupledRunRestartTime}.0!g;
     s!#Exp_APECoupleClimate_nml_RestartMeanEndTime#!${coupledRunRestartTime}.0!g;
EOF
    `	
    fi
    
    ocn_nml=${ocnDirPath}/${ocn_nml_template##*/}
    sed -e "${sedArgs}" ${ocn_nml_template} | sed -e "${sedArgs2}" > ${ocn_nml}

    #
    if [ $n -eq $StartCycleNum ] && $coupledRunSkipSCyc ; then
	echo "skip coupled run .."
    else
        echo "** Execute Dennou-OGCM  ******************************"

	cp    ${EXPDIR}/DCCM_AtmT21.conf ${atmDirPath}/DCCM.conf
	cp    ${EXPDIR}/DCCM_AtmT21.conf ${ocnDirPath}/DCCM.conf
	
	${MPIRUN}                                                   \
	-wdir ${atmDirPath} -env OMP_NUM_THREADS ${atm_THREADS_NUM} \
	-n ${atm_PE_NUM} ${atm_pe} -N=${atm_nml} :                  \
        -wdir ${ocnDirPath} -env OMP_NUM_THREADS ${ocn_THREADS_NUM} \
        -env LD_LIBRARY_PATH ${ocn_wdir}/lib                        \
        -n ${ocn_PE_NUM} ${ocn_pe} --N=${ocn_nml}                   \
	1> Stdout_couple_${exp_name} 2>Stderr_couple_${exp_name}

	if [ $? -ne 0 ]; then
	  echo "Exit stauts is 0.  Fail to run DCPCM. Exit.."; exit
	fi
        coupledRunEndTimeSec=`echo "$coupledRunEndTime*86400" | bc`
    fi

    coupledRunRestartTime=${coupledRunEndTime}
	
    #########################################################################
    # Run standalone ocean model with sea-ice model
    ########################################################################

    if [ $standaloneTimeIntrvPerCyc -gt 0 ] ; then
    
	echo "-- cycle=${n} (OGCM stadalone run) -- ${standaloneTimeIntrvPerCycle} [day]"

	sedArgs=`cat << EOF
      s!#gtool_historyauto_nml_IntValue#!1825.0!g;
      s!#gtool_historyauto_nml_OriginValue#!0.0!g;
      s!#gtool_historyauto_nml_TerminusValue#!${standaloneTimeIntrvPerCycle}!g;
      s!#OcnRestartFile_nml_InputFileName#!!g; 
      s!#OcnRestartFile_nml_OutputFileName#!RestartOcnData.nc!g;
      s!#OcnRestartFile_nml_IntValue#!9125.0!g;
      s!#SIceRestartFile_nml_InputFileName#!!g; 
      s!#SIceRestartFile_nml_OutputFileName#!RestartSIceData.nc!g;
      s!#SIceRestartFile_nml_IntValue#!9125.0!g;
      s!#TemporalInteg_nml_DelTimeHour#!${standaloneODelTimeHour}!g;
      s!#TemporalInteg_nml_RestartTimeVal#!0.0!g;
      s!#TemporalInteg_nml_InitYear#!2000!g;
      s!#TemporalInteg_nml_EndYear#!$((2000 + standaloneTimeIntrvPerCycle/365))!g;
      s!#TemporalInteg_nml_EndDay#!10!g;
      s!#BoundaryCondition_nml_ThermBCSurface#!PrescFlux_Han1984!g;
      s!#BoundaryCondition_nml_SaltBCSurface#!PrescFlux!g;
      s!#Exp_APECoupleClimate_nml_RunCycle#!${n}!g;
      s!#Exp_APECoupleClimate_nml_RunTypeName#!Standalone!g;
      s!#Exp_APECoupleClimate_nml_SfcBCDataDir#!${ocn_wdir}/cycle$((n))-couple/!g;
      s!#Exp_APECoupleClimate_nml_SfcBCMeanInitTime#!$((coupledRunEndTime - 438)).0!g;
      s!#Exp_APECoupleClimate_nml_SfcBCMeanEndTime#!${coupledRunEndTime}.0!g;
      s!#Exp_APECoupleClimate_nml_RestartDataDir#!${ocn_wdir}/cycle$((n))-couple/!g;
      s!#Exp_APECoupleClimate_nml_RestartMeanInitTime#!${coupledRunEndTime}.0!g;
      s!#Exp_APECoupleClimate_nml_RestartMeanEndTime#!${coupledRunEndTime}.0!g;
EOF
    `
	ocn_nml=${ocnDirPath_standalone}/${ocn_nml_template##*/}
	sed -e "${sedArgs}" ${ocn_nml_template} > ${ocn_nml}

	#
	${MPIRUN}                                                                             \
	    -wdir ${ocnDirPath_standalone} -env OMP_NUM_THREADS ${ocn_standalone_THREADS_NUM} \
	    -env LD_LIBRARY_PATH ${LD_LIBRARY_PATH}                                           \
	    -n ${ocn_standalone_PE_NUM}                                                       \
	    ${ocn_standalone_pe} --N=${ocn_nml}                                               \
	    1> Stdout_standalone_${exp_name} 2>Stderr_standalone_${exp_name}
	
	if [ $? -ne 0 ]; then
	    echo "Exit stauts is 0.  Fail to run Dennou-OGCM(stand-alone mode). Exit.."; exit
	fi
    fi

done
