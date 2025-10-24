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

list_analysis = ["MSF"]
if options[:list_analysis] != nil then
  list_analysis = options[:list_analysis].split(",")
end
LIST_ANALYSIS=list_analysis

DCPCM_TOOL_DIR="/data/ra000005/a04028/DCPCM/src/dcpcm_base_v1.0/tool"
CMD_MSF="#{DCPCM_TOOL_DIR}/postproc/ocn_analysis/MSF.rb"
CMD_OCNXTMEANSTATE="../../../script/ocn_xt_mean_state.rb"
CMD_GLOBALMEAN_TSERISE="#{DCPCM_TOOL_DIR}/postproc/ocn_analysis/collect-glmean-tseries.rb"
CMD_OCNDIAG="/home/ykawai/lib/Dennou-OGCM/bin/ocndiag"

CMD_SICEXTMEANSTATE="../../../script/sice_xt_mean_state.rb"

#---------------------------------------------

SolarConst = (options[:solar_const]==nil) ? 1380 : options[:solar_const].to_i
dir = (options[:exp_suffix]==nil) ? "S#{SolarConst}" : "S#{SolarConst}#{options[:exp_suffix]}"

TargetDir="./#{dir}/ocn"

BeginCycDiag=cyc_start
EndCycDiag=cyc_end

BeginCycMerge=cyc_start_merge
EndCycMerge=cyc_end_merge

INTRV_COUPLE = intrv_cyc_couple
INTRV_STANDALONE = intrv_cyc_standalone

CoupledCycIntDay = UNumeric[intrv_cyc_couple, "day"]
StandaloneCycIntDay = UNumeric[intrv_cyc_standalone, "day"]

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

def gen_ocndiag_nmlfile(ofname, dogcm_conf_name, cyc)

  if File.exist?(ofname) then
    FileUtils.rm_f(ofname)
  end

  tstart = (cyc-1)*INTRV_COUPLE
  tend   = tstart + INTRV_COUPLE
  outputIntDay = 146.0
  
  f = File.open(ofname, "w")
  f.print <<-END
 &ocndiag_nml
  configNmlFileDOGCM = "#{dogcm_conf_name}",
  TimeStart          = #{tstart.to_s},
  TimeEnd            = #{tend.to_s},
  TimeUnits          = 'day',
  TimeInt            = #{outputIntDay.to_s}, 
 /
 &gtool_historyauto_nml
   IntValue      = #{outputIntDay.to_s},  
   IntUnit       = 'day',                 
   OriginValue   = #{tstart.to_s},
   OriginUnit    = 'day',                             
   TerminusValue = #{tend.to_s}, 
   TerminusUnit  = 'day',                             
   FilePrefix    = 'diag/',
  /
 &gtool_historyauto_nml
   Name = 'DensPot, BVFreq', 
   Precision = 'float'
  /
 &gtool_historyauto_nml
   Name = 'EulerHT, BolusHT, IsoDiffHT, TotHT, MSF_GM, OcnHT, OcnHT_conv, NumDiffTend', 
   Precision = 'float'
  /
  END
  f.close
end

#------------------------------------------------------------


if LIST_ANALYSIS.include?("GlMeanTSerise") then
  Dir.chdir("#{TargetDir}"){
    cmd =   "ruby #{CMD_GLOBALMEAN_TSERISE}"            \
            + " -c #{BeginCycMerge}:#{EndCycMerge}"       \
            + " -r #{BeginCycDiag}:#{EndCycDiag}"         \
            + " -i #{INTRV_COUPLE},#{INTRV_STANDALONE} -p -P"
    exec_cmd(cmd)
  }
end

ret = Parallel.map(BeginCycDiag..EndCycDiag, :in_processes =>NProc){|i|
#(BeginCycDiag..EndCycDiag).each{|i|

  Dir.chdir("#{TargetDir}/cycle#{i}-couple"){
    currentDir = Dir::pwd
    if LIST_ANALYSIS.include?("MSF") then    
      p "run MSF.rb dir=#{currentDir} )"
      exec_cmd("ruby #{CMD_MSF}")
    end

    if LIST_ANALYSIS.include?("OcnDiagDo") then    
      p "run ocndiag dir=#{currentDir} )"
      ocn_conf = currentDir+"/dogcm_APEI07Couple_Pl42L60Opt.conf"
      ocn_conf = currentDir+"/dogcm_APEI07Couple_Pl64L60_I07SfcAlbMod.conf" if !File.exist?(ocn_conf)
      
      gen_ocndiag_nmlfile("ocndiag.conf", ocn_conf, i)
      if !File.exist?("diag") then
        p "create directory.."
        FileUtils.mkdir("diag")
      end
      
      exec_cmd("#{CMD_OCNDIAG} --N=./ocndiag.conf")
    end
  }
  
}

Dir.chdir(TargetDir){

=begin  
  if LIST_ANALYSIS.include?("EngyGlMean") then        
    for var in ["TotEngy", "IntEngy", "KinEngy", "PotEngy", "LatEngy", "Mass"]
      combine_ncfile(BeginCycMerge, EndCycMerge, "#{var}GlMean.nc", "#{var}.nc", var)
    end
  end
=end
  
  if LIST_ANALYSIS.include?("MeanState") then
    FileUtils.mkdir_p("./mean_state") unless FileTest.exist?("./mean_stae")

    list_addition = "EngyFlxLat,Stability"
    cmd="ruby #{CMD_OCNXTMEANSTATE} -c #{BeginCycMerge}:#{EndCycMerge} -s ./ -d ./mean_state -l #{list_addition} -v o2d_SfcTemp"
    p cmd
    exec_cmd(cmd)
  end  

  if LIST_ANALYSIS.include?("MeanStateSIce") then
    FileUtils.mkdir_p("./mean_state") unless FileTest.exist?("./mean_stae")

    list_addition = "EngyFlxLat,Stability"
    cmd="ruby #{CMD_SICEXTMEANSTATE} -c #{BeginCycMerge}:#{EndCycMerge} -s ./ -d ./mean_state -l #{list_addition} -v o2d_SfcTemp"
    p cmd
    exec_cmd(cmd)
  end  
  
}