#!/usr/env ruby
require "numru/ggraph"
require "optparse"
require "fileutils"
require "parallel"
require "benchmark"
include NumRu

opt = OptionParser.new
options = {}
opt.on("-c", "--cyc_range <param>",  "the range of cycle"){|v| options[:cyc_range] = v}
opt.on("-s", "-source_dir  <param>",  "the location of experimental directory "){|v| options[:exp_path] = v}
opt.on("-d", "--dist_dir  <param>",  "the distination for output data "){|v| options[:dist_path] = v}
opt.on("-l", "--list_addition  <param>",  "the list of additional target data "){|v| options[:list_add_target] = v}
opt.on("-v", "--var_name <param>",  "variable name of surface temperature"){|v| options[:var_name] = v}

opt.parse(ARGV)

cyc_start = 1
cyc_end   = 1
if options[:cyc_range] != nil then
  cyc_start = options[:cyc_range].split(":")[0].to_i
  cyc_end   = options[:cyc_range].split(":")[1].to_i
end

target_data_list=[""]
if options[:list_add_target] != nil then
  target_data_list = options[:list_add_target].split(",")
end

sfctempvar_name = "SurfTemp"
if options[:var_name] != nil then
  sfctempvar_name = options[:var_name]
end

TargetDir = (options[:exp_path] == nil) ? Dir.pwd : options[:exp_path]
DistDir = (options[:dist_path] == nil) ? Dir.pwd : options[:dist_path]
TargetDataList=target_data_list

BeginCombineCyc=cyc_start
EndCombineCyc=cyc_end

Grav = UNumeric[9.8, "m/s2"]

SfcTempVarName = sfctempvar_name
NProc=Parallel.processor_count

p "TargetDir: "+TargetDir
p "DistDir:"+DistDir
p "Period: #{cyc_start}:#{cyc_end}"
p "ListAddition:"+TargetDataList.join(",")

def merge_ncfile(ncname, varname, ovarname, cutOpt={}, meanOpt=[], rmissFlag=false, ofilename=ovarname+".nc")

  p "read NetCDF files (#{ncname}).."

  ncpathList = Hash.new
  (BeginCombineCyc..EndCombineCyc).each{|i|
    ncpathList[i] = "#{DistDir}/tmp#{i}_#{ofilename}"
  }
  p ncpathList
 
  ret = Parallel.map(BeginCombineCyc..EndCombineCyc, :in_processes =>NProc ){|i|  
#  (BeginCombineCyc..EndCombineCyc).each{|i|
    p "#{i}." if i%10 == 0
    p "#{TargetDir}\/cycle#{i}-couple\/#{ncname}"
    gphys = GPhys::IO.open("#{TargetDir}/cycle#{i}-couple/#{ncname}", varname)
#    gphys = GPhys::IO.open("#{TargetDir}\/cycle#{i}-couple\/#{ncname}", varname)
#    p gphys

    gphys = gphys.cut(cutOpt) if cutOpt.length > 0

    ps_read_flag = false    
    meanOpt.each{|opt|
      case opt      
      when "time"
        if rmissFlag then
          ocn_mask = (gphys < -1e20).where
          rmiss = NArray.float(1).fill(DCL::glrget('RMISS'))
          gphys[ocn_mask] = rmiss
          gphys.set_att('missing_value', rmiss)
        end
      when "sig"
        gp_sig_weight = GPhys::IO.open("#{TargetDir}/cycle#{i}-couple/U.nc", "sig_weight") 
        gp_ps = GPhys::IO.open(/#{TargetDir}\/cycle#{i}-couple\/H.nc/, "H")
        gphys = gp_ps/Grav*(gphys*gp_sig_weight).sum("sig")
      else
        gphys = gphys.mean(opt)
      end
    }
    
    ofile = NetCDF::create(ncpathList[i])
    tlen = gphys.axis("time").length
    GPhys::IO.write(ofile, gphys[false,1..tlen-1].copy.rename(ovarname))
    ofile.close
  }

  p "Output.."

  ofile = NetCDF::create("#{DistDir}/#{ofilename}")
  gphys = GPhys::IO.open(ncpathList.values, ovarname)
  gphys = gphys.mean("time") if meanOpt.include?("time")
  
  GPhys::IO.write(ofile, gphys)
#  NetCDF_Conventions.add_history(ofile, File.basename($0)+" "+ncpathList[0])  
  ofile.close

  p "remove tmp files.."
  FileUtils::rm_f(ncpathList.values)
end


=begin
merge_ncfile("history_sice.nc", "IceThick", "IceThick", {}, ["lon","time"])
merge_ncfile("history_sice.nc", "SnowThick", "SnowThick", {}, ["lon","time"])
merge_ncfile("history_sice.nc", "SIceV", "SIceV", {}, ["lon","time"])
merge_ncfile("history_sice.nc", "SIceTemp", "SIceTemp", {}, ["lon","time"], true)
merge_ncfile("history_sice.nc", "SIceSfcTemp", "SIceSfcTemp", {}, ["lon","time"], true)
=end
 merge_ncfile("history_sice.nc", "SfcHFlxAI", "SfcHFlxAI_mlon", {}, ["lon"])
#merge_ncfile("history_sice.nc", "SfcHFlxAI", "SfcHFlxAI", {}, ["lon","time"], true)

# merge_ncfile("SfcHFlxO.nc", "SfcHFlxO", "SfcHFlxO", {}, ["lon", "time"])
