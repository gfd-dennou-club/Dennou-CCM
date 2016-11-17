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

################################################

require 'fileutils'
require File.expand_path(File.dirname(__FILE__) + "/ocn_analysis/DennouOGCMUtil.rb")

MERGE_NCF_PROG = DCPCM_POSTPROC_TOOLSDIR + "/merge_ncf"
#DIAGVAROCN_PROG = DCPCM_POSTPROC_TOOLSDIR + "/diagVar_axisym"
GPCAT_PROG = DCPCM_POSTPROC_TOOLSDIR + "/gpcat"

OcnVarDef = DennouOGCMUtil::VarNameDef

################################################

def mkdir(newdirPath)
  if !File.exist?(newdirPath) then
    puts "Create directory: #{newdirPath}"
    Dir.mkdir(newdirPath)
  end
end

def merge_ATMNC(dirPath)
  puts "Perform merge_ncf.. (Dir=#{dirPath})"
  
  merge_nml = dirPath + "/../merge.nml"
  FileUtils.cp(merge_nml, dirPath)
  Dir.chdir(dirPath){
    `#{MERGE_NCF_PROG}`
  }
end

def diagVar_OCNNC(dirPath, suffix="")
  puts "Perform diagvar .. (Dir=#{dirPath})"
  diagvar_nml = dirPath + "/../diagVarConfig#{suffix}.nml"
  Dir.chdir(dirPath){
    puts "#{diagvar_nml}" 
    `#{DIAGVAROCN_PROG} --N=#{diagvar_nml}`
  }
end

def gpcat_GMSchemeVars_OCNNC(dirPath)
  puts "gpcat .. (Dir=#{dirPath})"

  Dir.chdir(dirPath){
    [ OcnVarDef::BolusV,"PTempIsopycDiffTend" ].each{|var|
      ncfname = "#{var}.nc"
      if File.exist?(ncfname) then
        puts "rm #{ncfname} .."
        FileUtils.rm(ncfname)
      end
      `#{GPCAT_PROG} -v #{var} -o #{var}.nc GMSchemeAnalysis.nc`
    }
  }
end

for n in CycleBegin..CycleEnd
  atmdir_couplerun = ATMDATADIR + "/cycle#{n}-couple"
  ocndir_couplerun = OCNDATADIR + "/cycle#{n}-couple"
  ocndir_standalonerun = OCNDATADIR + "/cycle#{n}-standalone"

   merge_ATMNC(atmdir_couplerun)
  
  # diagVar_OCNNC(ocndir_couplerun, "-couple")
  # gpcat_GMSchemeVars_OCNNC(ocndir_couplerun)  
  # diagVar_OCNNC(ocndir_standalonerun, "-standalone")
end


