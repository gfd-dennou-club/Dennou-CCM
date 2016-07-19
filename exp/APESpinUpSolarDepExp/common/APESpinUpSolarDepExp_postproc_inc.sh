#!/bin/bash
#------------------------------------------------------------------------------------------
# Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
#-------------------------------------------------------------------------------------------

RUBYCMD=ruby

#----------------------------------

TOPDIR=`pwd`
POSTPROC_STAGE1_RB=/home/ykawai/workspace/Dennou-CCM/tool/postproc/Dennou-CCM_postproc_stage1.rb
POSTPROC_STAGE2_RB=/home/ykawai/workspace/Dennou-CCM/tool/postproc/Dennou-CCM_postproc_stage2.rb
POSTPROC_STAGE3_RB=/home/ykawai/workspace/Dennou-CCM/tool/postproc/Dennou-CCM_postproc_stage3.rb

#----------------------------------
cd ${TOPDIR}

echo "SolarConst=${SolarConst}"
echo "TOPDIR=${TOPDIR}"

if [ ${STAGE1} = "T" ]; then
    echo "Dennou-CCM postprocess stage1 --------------------"
    $RUBYCMD ${POSTPROC_STAGE1_RB} \
	     --topdir=${TOPDIR}/ --atmdir=run_S${SolarConst}/atm --ocndir=run_S${SolarConst}/ocn \
	     --cycles=${CycStart}:${CycEnd}
fi

if [ ${STAGE2} = "T" ]; then
    echo "Dennou-CCM postprocess stage2 --------------------"
    $RUBYCMD ${POSTPROC_STAGE2_RB} \
	     --topdir=${TOPDIR}/ --atmdir=run_S${SolarConst}/atm --ocndir=run_S${SolarConst}/ocn \
	     --cycles=${CycStart}:${CycEnd}
fi

if [ ${STAGE3} = "T" ]; then
    echo "Dennou-CCM postprocess stage3 --------------------"
    $RUBYCMD ${POSTPROC_STAGE3_RB} \
          --topdir=${TOPDIR}/ --atmdir=run_S${SolarConst}/atm --ocndir=run_S${SolarConst}/ocn \
          --cyclesTrans=${CycStart}:${CycEnd} \
          --cyclesMean=${MeanCycStart}:${MeanCycEnd}
fi

