#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------
# * Run merge_ncf for DCPAM output data
#   and diagvar for Dennou-OGCM output data. 
#
# Configuration **********************************

require 'optparse'

TOPDIR = `pwd`

OPTS = ARGV.getopts('', 'topdir:', 'atmdir:', 'ocndir:', 'cycles:').inject({}) {
  |hash,(k,v)| hash[k.to_sym] = v; hash
}

p OPTS

cycles = OPTS[:cycles].split(":")
CycleBegin = cycles[0].to_i
CycleEnd = cycles[1].to_i

ATMDATADIR = OPTS[:topdir] + OPTS[:atmdir]
OCNDATADIR = OPTS[:topdir] + OPTS[:ocndir]


DCPCM_POSTPROC_DIR = File.expand_path(File.dirname(__FILE__))
DCPCM_POSTPROC_TOOLSDIR = File.expand_path(File.dirname(__FILE__)) + "/tools"

RUBY = "ruby"

###############################################

#-----------------------------------------------

DIAGATM_PROG =DCPCM_POSTPROC_DIR + "/atm_analysis/Diagnose.rb"
EngyFlxLatAtm_PROG = DCPCM_POSTPROC_DIR + "/atm_analysis/EngyFlxLat.rb"
EngyChkAtm_PROG = DCPCM_POSTPROC_DIR + "/atm_analysis/EngyCheck.rb"
SurfFlx_PROG = DCPCM_POSTPROC_DIR + "/ocn_analysis/SurfFlx.rb"
EngyFlxLatOcn_PROG = DCPCM_POSTPROC_DIR + "/ocn_analysis/EngyFlxLat.rb"
ATM_GlobalMeanQuants_PROG = DCPCM_POSTPROC_DIR + "/atm_analysis/GlobalMeanQuants.rb"
OCN_GlobalMeanQuants_PROG = DCPCM_POSTPROC_DIR + "/ocn_analysis/GlobalMeanQuants.rb"


################################################

require 'fileutils'

################################################

def mkdir(newdirPath)
  if !File.exist?(newdirPath) then
    puts "Create directory: #{newdirPath}"
    Dir.mkdir(newdirPath)
  end
end

def diagnose_ATMNC(dirPath)
  puts "Perform diagnose .. (Dir=#{dirPath})"
#  diagvar_nml = dirPath+"/diagVarConfig.nml"
#  FileUtils.cp(DIAGVAR_NML_ORI, diagvar_nml)
  Dir.chdir(dirPath){
    `#{RUBY} #{DIAGATM_PROG}`
  }
end

def run_AtmGlobalMeanQuants(dirPath)
  puts "Perform GlobalMeanQuants.rb .. (Dir=#{dirPath})"
  Dir.chdir(dirPath){
      `#{RUBY} #{ATM_GlobalMeanQuants_PROG}`
  }
end

def run_OcnGlobalMeanQuants(dirPath)
  puts "Perform GlobalMeanQuants.rb .. (Dir=#{dirPath})"
  Dir.chdir(dirPath){
      `#{RUBY} #{OCN_GlobalMeanQuants_PROG}`
  }
end

def run_EngyFlxLatAtm(dirPath)
  puts "Perform EngyFlxLat.rb .. (Dir=#{dirPath})"
  Dir.chdir(dirPath){
    `#{RUBY} #{EngyFlxLatAtm_PROG}`
  }
end

def run_EngyChkAtm(dirPath)
  puts "Perform EngyCheck.rb .. (Dir=#{dirPath})"
  Dir.chdir(dirPath){
    `#{RUBY} #{EngyChkAtm_PROG}`
  }
end

def run_SurfFlx(dirPath)
  puts "Perform SurfFlx.rb .. (Dir=#{dirPath})"
  Dir.chdir(dirPath){
    `#{RUBY} #{SurfFlx_PROG}`
  }
end

def run_EngyFlxLatOcn(dirPath)
  puts "Perform EngyFlxLat.rb .. (Dir=#{dirPath})"
  Dir.chdir(dirPath){
    `#{RUBY} #{EngyFlxLatOcn_PROG}`
  }
end

for n in CycleBegin..CycleEnd
  atmdir_couplerun = ATMDATADIR + "/cycle#{n}-couple"
  ocndir_couplerun = OCNDATADIR + "/cycle#{n}-couple"
  ocndir_standlonerun = OCNDATADIR + "/cycle#{n}-standalone"

  diagnose_ATMNC(atmdir_couplerun)
  run_AtmGlobalMeanQuants(atmdir_couplerun)
  run_EngyFlxLatAtm(atmdir_couplerun)
  run_EngyChkAtm(atmdir_couplerun)

  run_OcnGlobalMeanQuants(ocndir_couplerun)
#  run_SurfFlx(ocndir_couplerun)
#  run_EngyFlxLatOcn(ocndir_couplerun)
end

