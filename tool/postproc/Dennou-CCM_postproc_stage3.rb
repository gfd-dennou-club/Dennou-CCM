#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------
# * Run merge_ncf for DCPAM output data
#   and diagvar for Dennou-OGCM output data. 
#
# Configuration **********************************

require 'optparse'

TOPDIR = `pwd`

OPTS = ARGV.getopts('', 'topdir:', 'atmdir:', 'ocndir:', 'cyclesTrans:', 'cyclesMean:').inject({}) {
  |hash,(k,v)| hash[k.to_sym] = v; hash
}

p OPTS

cycles = OPTS[:cyclesTrans].split(":")
CycleBegin_Transit = cycles[0].to_i
CycleEnd_Transit = cycles[1].to_i

cycles = OPTS[:cyclesMean].split(":")
CycleBegin_Mean = cycles[0].to_i
CycleEnd_Mean = cycles[1].to_i

ATMDATADIR = OPTS[:topdir] + OPTS[:atmdir]
OCNDATADIR = OPTS[:topdir] + OPTS[:ocndir]


DCPCM_POSTPROC_DIR = File.expand_path(File.dirname(__FILE__))
DCPCM_POSTPROC_TOOLSDIR = File.expand_path(File.dirname(__FILE__)) + "/tools"

RUBY = "ruby"

###############################################

#-----------------------------------------------

################################################

require 'fileutils'

MERGE_COUPLERUNDATA_DISTDIRNAME = "transition"
MEAN_COUPLERUNDATA_DISTDIRNAME = "mean_state"
MEAN_STANDALONERUNDATA_DISTDIRNAME = "mean_state_standalone"

#-- Atmoshere --------------------------------------------------

MERGE_COPLEDRUNDATA_PROG = DCPCM_POSTPROC_DIR + "/common/merge_couplerun_data.rb"
MERGE_COUPLERUNDATA_LIST_ATM = \
{ "GlobalMeanQuants.nc" => "o2d_SfcTemp,o2d_SfcAlbedo,OLR,mOSR,RadTOA", 
  "TotEngy.nc"=>"TotEngy", "IntEngy.nc"=>"IntEngy", "KinEngy.nc"=>"KinEngy", 
  "PotEngy.nc"=>"PotEngy", "LatEngy.nc"=>"LatEngy" }

MEAN_COPLEDRUNDATA_PROG = DCPCM_POSTPROC_DIR + "/common/mean_couplerun_data.rb"
MEAN_STANDALONERUNDATA_PROG = DCPCM_POSTPROC_DIR + "/common/mean_standalonerun_data.rb"

MEAN_XT_COUPLERUNDATA_LIST_ATM = {
  "U.nc" => "U", "V.nc"=>"V", "Temp.nc"=>"Temp", "QH2OVap.nc"=>"QH2OVap", "Diagnose.nc"=>"MSF", 
  "o2d_SfcTemp.nc"=>"o2d_SfcTemp", "o2d_SfcAlbedo.nc"=>"o2d_SfcAlbedo",
  "PRCP.nc"=>"PRCP", "PotTemp.nc"=>"PotTemp" }
MEAN_XT_COUPLERUNDATA_LIST2_ATM = {
  "CloudCoverforRad.nc"=>"CloudCoverforRad" }
MEAN_T_COUPLERUNDATA_LIST_ATM = {  "EngyFlx.nc"=>"totStatEnFlxLat,dryStatEnFlxLat,moistStatEnFlxLat" }

#-- Ocean ------------------------------------------------------

MERGE_COUPLERUNDATA_LIST_OCN = { 
  "GlobalMeanQuants.nc" => "PTemp,Salt", 
  "GlobalMeanQuants_SfcFlx.nc" => "SfcHFlxO,FreshWtFlxS", 
}

MEAN_XT_COUPLERUNDATA_LIST_OCN = {
  "U.nc" => "U", "V.nc"=>"V", "PTemp.nc"=>"PTemp", "Salt.nc"=>"Salt", 
  "a2o_WindStressX.nc"=>"a2o_WindStressX",
  "a2o_WindStressY.nc"=>"a2o_WindStressY",
  "a2o_SfcAirTemp.nc"=>"a2o_SfcAirTemp",  
  "SfcHFlxO.nc"=>"SfcHFlxO", "FreshWtFlxS.nc"=>"FreshWtFlxS"
}

MEAN_T_COUPLERUNDATA_LIST_OCN = {
  "MassStreamFunc.nc"=>"MassStreamFunc", 
  "SGSEddyMixAnalysis.nc" \
  =>  "EngyFlxLatTot,EngyFlxLatEuler,EngyFlxLatBolus,EngyFlxLatIPDiff,EngyFlxLatTotHDiv,PT_TendVINT_IPDiff,PT_TendVINT_GM,PT_TendVINT_ADV,PT_TendVINT_HDiff,PT_TendVINT_VDiff,PT_TendVINT_TOT",
}

MEAN_XT_STANDALONERUNDATA_LIST_OCN = {
  "SfcHFlxO.nc" => "SfcHFlxO", 
  "FreshWtFlxS.nc" => "FreshWtFlxS", 
 }
MEAN_T_STANDALONERUNDATA_LIST_OCN = {
  "SGSEddyMixAnalysis.nc" \
  =>  "EngyFlxLatTot,EngyFlxLatEuler,EngyFlxLatBolus,EngyFlxLatIPDiff,EngyFlxLatTotHDiv,PT_TendVINT_IPDiff,PT_TendVINT_GM,PT_TendVINT_ADV,PT_TendVINT_HDiff,PT_TendVINT_VDiff,PT_TendVINT_TOT",
}

#-- Sea ice ------------------------------------------------------

SICE_NCFNAME="history_sice.nc"

MERGE_COUPLERUNDATA_LIST_SICE = {
  "GlobalMeanQuants_SIce.nc" => "IceThick,SnowThick",
}

MEAN_XT_COUPLERUNDATA_LIST_SICE = {
  "history_sice.nc"=>"SIceSfcTemp,IceThick,SnowThick,SIceCon",
}

MEAN_XT_COUPLERUNDATA_LIST2_SICE = {
  "history_sice.nc"=>"Wice", 
}

MEAN_T_COUPLERUNDATA_LIST_SICE = {}


# Main part #############################################################################################

def mkdir(newdirPath)
  if !File.exist?(newdirPath) then
    puts "Create directory: #{newdirPath}"
    Dir.mkdir(newdirPath)
  end
end

def merge_CoupledRunData(dirPath, cycleBegin, cycleEnd, nc, varList, distDirPath)
  puts "Perform merge_CoupledRunData .. (Dir=#{dirPath})"
  puts "Cycle=#{cycleBegin}..#{cycleEnd}"
  puts "VarList=#{varList}"
  Dir.chdir(dirPath){
    varList.split(",").each{|var|
      `#{RUBY} #{MERGE_COPLEDRUNDATA_PROG} #{cycleBegin} #{cycleEnd} #{nc} #{var} . #{distDirPath}`
    }
  }  
end

def merge_CoupledStandAloneRunData(dirPath, cycleBegin, cycleEnd, nc, varList, distDirPath)
  puts "Perform merge_CoupledRunData .. (Dir=#{dirPath})"
  puts "Cycle=#{cycleBegin}..#{cycleEnd}"
  puts "VarList=#{varList}"
  Dir.chdir(dirPath){
    varList.split(",").each{|var|
      `#{RUBY} #{MERGE_COPLEDRUNDATA_PROG} #{cycleBegin} #{cycleEnd} #{nc} #{var} . #{distDirPath}`
    }
  }  
end

def mean_CoupledRunData(dirPath, cycleBegin, cycleEnd, nc, varList, meanAxisNames, distDirPath)
  puts "Perform mean_CoupledRunData .. (Dir=#{dirPath})"
  puts "Cycle=#{cycleBegin}..#{cycleEnd}"
  puts "VarList=#{varList}"
  puts "Axis=#{meanAxisNames}"
  
  Dir.chdir(dirPath){
    varList.split(",").each{|var|
      `#{RUBY} #{MEAN_COPLEDRUNDATA_PROG} #{cycleBegin} #{cycleEnd} #{nc} #{var} #{meanAxisNames} . #{distDirPath}`
    }
  }  
end

def mean_StandaloneRunData(dirPath, cycleBegin, cycleEnd, nc, varList, meanAxisNames, distDirPath)
  puts "Perform mean_StandaloneRunData .. (Dir=#{dirPath})"
  puts "Cycle=#{cycleBegin}..#{cycleEnd}"
  puts "VarList=#{varList}"
  puts "Axis=#{meanAxisNames}"

  Dir.chdir(dirPath){
    varList.split(",").each{|var|
      p var
      `#{RUBY} #{MEAN_STANDALONERUNDATA_PROG} #{cycleBegin} #{cycleEnd} #{nc} #{var} #{meanAxisNames} . #{distDirPath}`
    }
  }  
end

#--- Atmospheric component ----------------------------------------------------------------------------------

distDirPath = ATMDATADIR+"/#{MERGE_COUPLERUNDATA_DISTDIRNAME}"
MERGE_COUPLERUNDATA_LIST_ATM.each{|nc,varList|
  mkdir(distDirPath)
  merge_CoupledRunData(ATMDATADIR, CycleBegin_Transit, CycleEnd_Transit, nc, varList, distDirPath)
}


#
distDirPath = ATMDATADIR+"/#{MEAN_COUPLERUNDATA_DISTDIRNAME}"
MEAN_XT_COUPLERUNDATA_LIST_ATM.each{|nc,varList|
  mkdir(distDirPath)
  mean_CoupledRunData(ATMDATADIR, CycleBegin_Mean, CycleEnd_Mean, nc, varList, "time,lon", distDirPath)
}
MEAN_T_COUPLERUNDATA_LIST_ATM.each{|nc,varList|
  mkdir(distDirPath)
  mean_CoupledRunData(ATMDATADIR, CycleBegin_Mean, CycleEnd_Mean, nc, varList, "time", distDirPath)
}

#--- Ocean component ----------------------------------------------------------------------------------

p MERGE_COUPLERUNDATA_LIST_OCN
#
distDirPath = OCNDATADIR+"/#{MERGE_COUPLERUNDATA_DISTDIRNAME}"

MERGE_COUPLERUNDATA_LIST_OCN.each{|nc,varList|
  mkdir(distDirPath)
  merge_CoupledRunData(OCNDATADIR, CycleBegin_Transit, CycleEnd_Transit, nc, varList, distDirPath)
}

distDirPath = OCNDATADIR+"/#{MEAN_COUPLERUNDATA_DISTDIRNAME}"
MEAN_XT_COUPLERUNDATA_LIST_OCN.each{|nc,varList|
  mkdir(distDirPath)
  mean_CoupledRunData(OCNDATADIR, CycleBegin_Mean, CycleEnd_Mean, nc, varList, "time,lon", distDirPath)
}
=begin
MEAN_T_COUPLERUNDATA_LIST_OCN.each{|nc,varList|
  mkdir(distDirPath)
  mean_CoupledRunData(OCNDATADIR, CycleBegin_Mean, CycleEnd_Mean, nc, varList, "time", distDirPath)
}

distDirPath = OCNDATADIR+"/#{MEAN_STANDALONERUNDATA_DISTDIRNAME}"
MEAN_T_STANDALONERUNDATA_LIST_OCN.each{|nc,varList|
  mkdir(distDirPath)
  mean_StandaloneRunData(OCNDATADIR, CycleBegin_Mean, CycleEnd_Mean, nc, varList, "time", distDirPath)
}
=end

#--- Sea-ice component ----------------------------------------------------------------------------------

p MERGE_COUPLERUNDATA_LIST_SICE
#
distDirPath = OCNDATADIR+"/#{MERGE_COUPLERUNDATA_DISTDIRNAME}"

MERGE_COUPLERUNDATA_LIST_SICE.each{|nc,varList|
  mkdir(distDirPath)
  merge_CoupledRunData(OCNDATADIR, CycleBegin_Transit, CycleEnd_Transit, nc, varList, distDirPath)
}

distDirPath = OCNDATADIR+"/#{MEAN_COUPLERUNDATA_DISTDIRNAME}"
MEAN_XT_COUPLERUNDATA_LIST_SICE.each{|nc,varList|
  mkdir(distDirPath)
  mean_CoupledRunData(OCNDATADIR, CycleBegin_Mean, CycleEnd_Mean, nc, varList, "time,lon", distDirPath)
}
MEAN_T_COUPLERUNDATA_LIST_SICE.each{|nc,varList|
  mkdir(distDirPath)
  mean_CoupledRunData(OCNDATADIR, CycleBegin_Mean, CycleEnd_Mean, nc, varList, "time", distDirPath)
}
