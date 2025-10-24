#!/usr/env ruby

############################

require "numru/ggraph"
require "optparse"
require "fileutils"
require "parallel"
include NumRu

opt = OptionParser.new
options = {}
opt.on("-c", "--cyc_range_diagnoe <param>",  "the range of cycle"){|v| options[:cyc_range] = v}
opt.on("-m", "--cyc_range_merge <param>",  "the range of cycle"){|v| options[:cyc_range_merge] = v}
opt.on("-s", "--solar_const  <param>",  "the location of experimental directory "){|v| options[:solar_const] = v}
opt.on("-e", "--exp_suffix  <param>",  "the suffix of experimetal directory"){|v| options[:exp_suffix] = v}
opt.on("-i", "--interval_cyc <param>", "the interval of each cycle. (day)"){|v| options[:interval_cyc] = v}
opt.on("-l", "--list_analysis <param>", "array of the name of analysis"){|v| options[:list_analysis] = v}
#opt.on("-v", "--var_name <param>",  "variable name of surface temperature"){|v| options[:var_name] = v}
opt.parse(ARGV)

cyc_start = 1
cyc_end   = 1
if options[:cyc_range] != nil then
  cyc_start = options[:cyc_range].split(":")[0].to_i
  cyc_end   = options[:cyc_range].split(":")[1].to_i
end

cyc_start_merge=cyc_start
cyc_end_merge=cyc_end
if options[:cyc_range_merge] != nil then
  cyc_start_merge = options[:cyc_range_merge].split(":")[0].to_i
  cyc_end_merge   = options[:cyc_range_merge].split(":")[1].to_i
end

intrv_cyc_couple     = 730.0
intrv_cyc_standalone = 50.0 * 365.0
if options[:interval_cyc] != nil then
  intrv_cyc_couple = options[:interval_cyc].split(",")[0].to_f
  intrv_cyc_standalone = options[:interval_cyc].split(",")[1].to_f
end

list_analysis = ["MSF", "Iceline", "EngyGlMean"]
if options[:list_analysis] != nil then
  list_analysis = options[:list_analysis].split(",")
end
LIST_ANALYSIS=list_analysis

SolarConst = (options[:solar_const]==nil) ? 1380 : options[:solar_const].to_i
dir = (options[:exp_suffix]==nil) ? "S#{SolarConst}" : "S#{SolarConst}#{options[:exp_suffix]}"
TargetDir="./#{dir}/atm"
BeginCycDiag=cyc_start
EndCycDiag=cyc_end

DCPCM_TOOL_DIR="/data/ra000005/a04028/DCPCM/src/dcpcm_base_v1.0/tool"
CMD_GLOBALMEAN="#{DCPCM_TOOL_DIR}/postproc/atm_analysis/GlobalMeanQuants.rb"
CMD_CHECKTOT="#{DCPCM_TOOL_DIR}/postproc/atm_analysis/EngyCheck.rb"
CMD_ENGYFLXLAT="#{DCPCM_TOOL_DIR}/postproc/atm_analysis/EngyFlxLat.rb"
CMD_DIAGNOSE="#{DCPCM_TOOL_DIR}/postproc/atm_analysis/Diagnose.rb"
CMD_MERGENCF="/data/ra000005/a04028/DCPCM/src/util_merge-2011-03-28-2/merge_ncf"
CMD_ICELINE="../../../script/iceline_lat_tseries.rb"
CMD_ATMXTMEANSTATE="../../../script/atm_xt_mean_state.rb"

MERGENCF_CNFFILE = "../../../../common/dcpam/merge-4proc.nml"

#---------------------------------------------

BeginCycMerge=cyc_start_merge
EndCycMerge=cyc_end_merge

CoupledCycIntDay = UNumeric[intrv_cyc_couple, "day"]
StandaloneCycIntDay = UNumeric[intrv_cyc_standalone, "day"]

varListAtm=["PTemp", "Salt"]
varListSfcFlx=["SfcHFlxO", "FreshWtFlxS"]
varListIce=["IceThick", "SnowThick", "SIceEn", "SIceCon"]


NProc=Parallel.processor_count
p "NProc=#{NProc}"

###########################

require "fileutils"
require "open3"
require "numru/ggraph"
include NumRu

def exec_cmd(cmd, cyc=-1)
  lines_o = "";  lines_e = ""
  Open3.popen3(cmd) do |i, o, e, w|
    o.each do |line| lines_o <<  "#{cyc}:"+line  end
    e.each do |line| lines_e <<  "#{cyc}:"+line  end
  end
  puts "#{lines_o}";   puts "#{lines_e}"  
end

def combine_ncfile(beginCyc, endCyc, ofname, ifname, varname)

  nCyc = endCyc - beginCyc + 1
  nSubCyc = 10
  nBlock = nCyc / nSubCyc
  nBlock += 1 if (nCyc % nSubCyc != 0)
  
  index_info = Hash.new
  subfnames = []
  for i in 0..nBlock-1
    beginSubCyc = i*nSubCyc + beginCyc
    endSubCyc = [beginSubCyc + nSubCyc - 1, endCyc].min
    index_info[i] = { "begin_cyc"=>beginSubCyc, "end_cyc"=>endSubCyc }
    subfnames.push( "tmp#{i}-#{ofname}" )
  end

  ret = Parallel.map(0..nBlock-1, :in_processes => NProc){|i|
    index = index_info[i]
    fnames = []
    p "block:#{i}"
    (index["begin_cyc"]..index["end_cyc"]).each{|i|
      fnames.push("cycle#{i}-couple/#{ifname}")
    }

    FileUtils.rm_f(subfnames[i]) if File.exist?(subfnames[i])    
    exec_cmd( "gpcat -o #{subfnames[i]} -v #{varname} #{fnames.join(" ")}",
              "#{index["begin_cyc"]}..#{index["end_cyc"]}"            )  
  }
  
  #
  FileUtils.rm_f(ofname) if File.exist?(ofname)
  exec_cmd("gpcat -o #{ofname} -v #{varname} #{subfnames.join(" ")}")  

  #--
  for i in 0..nBlock-1  
    FileUtils.rm_f(subfnames[i]) if File.exist?(subfnames[i])    
  end
end


#------------------------------------------------------------

cycles = (BeginCycDiag..EndCycDiag).map{|i| i}
ret = Parallel.map(cycles, :in_processes =>NProc){|i|
  Dir.chdir("#{TargetDir}/cycle#{i}-couple"){
    
    p "run merge_ncf dir=#{Dir::pwd} )"
    `cp #{MERGENCF_CNFFILE} merge.nml`
    exec_cmd(CMD_MERGENCF, i) 

    if LIST_ANALYSIS.include?("EngyGlMean") then            
      p "run chceck total energy dir=#{Dir::pwd} )"
      exec_cmd("ruby #{CMD_CHECKTOT}", i)
    end

    if LIST_ANALYSIS.include?("GlobalMean") then                
      p "run GlobalMeanQuants.rb( dir=#{Dir::pwd} )"
      exec_cmd("ruby #{CMD_GLOBALMEAN}")
    end
    
    if LIST_ANALYSIS.include?("MSF") then    
      p "run Diagnose.rb dir=#{Dir::pwd} )"
      exec_cmd("ruby #{CMD_DIAGNOSE} --const_util ConstUtil_INTH07", i)
    end

    if LIST_ANALYSIS.include?("EngyFlxLat") then    
      p "run EngyFlxLat.rb dir=#{Dir::pwd} )"
      exec_cmd("ruby #{CMD_ENGYFLXLAT} --const_util ConstUtil_INTH07", i)
    end
  }
}

Dir.chdir(TargetDir){

  if LIST_ANALYSIS.include?("EngyGlMean") then        
    for var in ["TotEngy", "IntEngy", "KinEngy", "PotEngy", "LatEngy", "Mass"]
      combine_ncfile(BeginCycMerge, EndCycMerge, "#{var}GlMean.nc", "#{var}.nc", var)
    end
  end
  if LIST_ANALYSIS.include?("GlobalMean") then        
    for var in ["OLR", "mOSR"]
      combine_ncfile(BeginCycMerge, EndCycMerge, "#{var}GlMean.nc", "GlobalMeanQuants.nc", var)
    end
  end
  if LIST_ANALYSIS.include?("Iceline") then      
    exec_cmd("ruby #{CMD_ICELINE} -c #{BeginCycMerge}:#{EndCycMerge} -r #{BeginCycDiag}:#{EndCycDiag}  -v o2d_SfcTemp -i #{CoupledCycIntDay}:#{StandaloneCycIntDay}")
  end

  if LIST_ANALYSIS.include?("MeanState") then
    FileUtils.mkdir_p("./mean_state") unless FileTest.exist?("./mean_stae")

    list_addition = "EngyFlxLat"
    cmd="ruby #{CMD_ATMXTMEANSTATE} -c #{BeginCycMerge}:#{EndCycMerge} -s ./ -d ./mean_state -l #{list_addition} -v o2d_SfcTemp"
    p cmd
    exec_cmd(cmd)
  end  
}

#for var in ["U", "Temp"]
#  combine_ncfile(BeginCyc, EndCyc, "#{var}.nc","#{var}_rank00*.nc", var)
#end

#for var in ["totStatEnFlxLat", "dryStatEnFlxLat", "moistStatEnFlxLat"]
#  combine_ncfile(BeginCyc, EndCyc, "#{var}.nc","EngyFlx.nc", var)
#end

#for var in varListAtm
#  combine_ncfile(BeginCyc, EndCyc, "#{var}GlMean.nc", "GlobalMeanQuants.nc", var)
#end

#for var in varListSfcFlx
#  combine_ncfile(BeginCyc, EndCyc, "#{var}GlMean.nc", "GlobalMeanQuants_SfcFlx.nc", var)
#end

#for var in varListIce
#  combine_ncfile(BeginCyc, EndCyc, "#{var}GlMean.nc", "GlobalMeanQuants_SIce.nc", var)
#end
