#!/usr/env ruby
require "numru/ggraph"
require "optparse"
require "fileutils"
include NumRu

opt = OptionParser.new
options = {}
opt.on("-c", "--cyc_range <param>",  "the range of cycle"){|v| options[:cyc_range] = v}
opt.on("-s", "-source_dir  <param>",  "the location of experimental directory "){|v| options[:exp_path] = v}
opt.on("-d", "--dist_dir  <param>",  "the distination for output data "){|v| options[:dist_path] = v}
opt.on("-l", "--list_addition  <param>",  "the list of additional target data "){|v| options[:list_add_target] = v}
#opt.on("-v", "--var_name <param>",  "variable name of surface temperature"){|v| options[:var_name] = v}
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

TargetDir = (options[:exp_path] == nil) ? Dir.pwd : options[:exp_path]
DistDir = (options[:dist_path] == nil) ? Dir.pwd : options[:dist_path]
TargetDataList=target_data_list

BeginCombineCyc=cyc_start
EndCombineCyc=cyc_end

Grav = UNumeric[9.8, "m/s2"]

p "TargetDir: "+TargetDir
p "DistDir:"+DistDir
p "Period: #{cyc_start}:#{cyc_end}"
p "ListAddition:"+TargetDataList.join(",")

def merge_ncfile(ncname, varname, ovarname, cutOpt={}, meanOpt=[], ofilename=ovarname+".nc")
  ncpathList = []
  p "read NetCDF files (#{ncname}).."

  for i in BeginCombineCyc..EndCombineCyc
    p "#{i}." if i%10 == 0
    p "#{TargetDir}\/cycle#{i}-couple\/#{ncname}"
    gphys = GPhys::IO.open(/#{TargetDir}\/cycle#{i}-couple\/#{ncname}/, varname)
#    gphys = GPhys::IO.open("#{TargetDir}\/cycle#{i}-couple\/#{ncname}", varname)
#    p gphys

    gphys = gphys.cut(cutOpt) if cutOpt.length > 0

    ps_read_flag = false    
    meanOpt.each{|opt|
      case opt      
      when "time"
      when "sig"
        gp_sig_weight = GPhys::IO.open("#{TargetDir}/cycle#{i}-couple/U.nc", "sig_weight") 
        gp_ps = GPhys::IO.open(/#{TargetDir}\/cycle#{i}-couple\/Ps_rank(\d\d\d\d\d\d).nc/, "Ps")
        gphys = gp_ps/Grav*(gphys*gp_sig_weight).sum("sig")
      else
        gphys = gphys.mean(opt)
      end
    }
    
    ncpath = "#{DistDir}/tmp#{i}_#{ofilename}"
    ncpathList.push(ncpath)
    ofile = NetCDF::create(ncpath)
    tlen = gphys.axis("time").length
    GPhys::IO.write(ofile, gphys[false,1..tlen-1].copy.rename(ovarname))
    ofile.close
  end

  p "Output.."

  ofile = NetCDF::create("#{DistDir}/#{ofilename}")
  gphys = GPhys::IO.open(ncpathList, ovarname)
  gphys = gphys.mean("time") if meanOpt.include?("time")
  
  GPhys::IO.write(ofile, gphys)
#  NetCDF_Conventions.add_history(ofile, File.basename($0)+" "+ncpathList[0])  
  ofile.close

  p "remove tmp files.."
  FileUtils::rm_f(ncpathList)
end

#=begin
merge_ncfile("OLRA_rank(\\d\\d\\d\\d\\d\\d).nc", "OLRA", "OLRA", {}, ["lon","time"])
merge_ncfile("OSRA_rank(\\d\\d\\d\\d\\d\\d).nc", "OSRA", "OSRA", {}, ["lon","time"])
merge_ncfile("SLRA_rank(\\d\\d\\d\\d\\d\\d).nc", "SLRA", "SLRA", {}, ["lon","time"])
merge_ncfile("SSRA_rank(\\d\\d\\d\\d\\d\\d).nc", "SSRA", "SSRA", {}, ["lon","time"])

merge_ncfile("EvapA_rank(\\d\\d\\d\\d\\d\\d).nc", "EvapA", "LatHFlxA", {}, ["lon","time"])
merge_ncfile("SensA_rank(\\d\\d\\d\\d\\d\\d).nc", "SensA", "SenHFlxA", {}, ["lon","time"])

merge_ncfile("TauX_rank(\\d\\d\\d\\d\\d\\d).nc", "TauX", "TauX", {}, ["lon","time"])
merge_ncfile("TauY_rank(\\d\\d\\d\\d\\d\\d).nc", "TauY", "TauY", {}, ["lon","time"])

merge_ncfile("PRCP_rank(\\d\\d\\d\\d\\d\\d).nc", "PRCP", "PRCP", {}, ["lon","time"])

merge_ncfile("Ps_rank(\\d\\d\\d\\d\\d\\d).nc", "Ps", "Ps", {}, ["lon","time"])
merge_ncfile("SurfTemp_rank(\\d\\d\\d\\d\\d\\d).nc", "SurfTemp", "SfcTemp", {}, ["lon","time"])

merge_ncfile("^Temp_rank(\\d\\d\\d\\d\\d\\d).nc", "Temp", "Temp", {}, ["lon","time"])
merge_ncfile("U_rank(\\d\\d\\d\\d\\d\\d).nc", "U", "U", {}, ["lon","time"])
merge_ncfile("V_rank(\\d\\d\\d\\d\\d\\d).nc", "V", "V", {}, ["lon","time"])
merge_ncfile("QH2OVap_rank(\\d\\d\\d\\d\\d\\d).nc", "QH2OVap", "QH2OVap", {}, ["lon","time"])
merge_ncfile("SigDot_rank(\\d\\d\\d\\d\\d\\d).nc", "SigDot", "SigDot", {}, ["lon","time"])
merge_ncfile("Diagnos(\\w).nc", "MSF", "MSF", {}, ["time"])

#=end
#merge_ncfile("QH2OVap_rank(\\d\\d\\d\\d\\d\\d).nc", "QH2OVap", "PWV", {}, ["sig", "lon","time"])

if TargetDataList.include?("EngyFlxLat") then
  merge_ncfile("EngyFl(\\w).nc", "dryStatEnFlxLat", "dryStatEnFlxLat", {}, ["time"])
  merge_ncfile("EngyFl(\\w).nc", "moistStatEnFlxLat", "moistStatEnFlxLat", {}, ["time"])
  merge_ncfile("EngyFl(\\w).nc", "latentEnFlxLat", "latentEnFlxLat", {}, ["time"])
end

