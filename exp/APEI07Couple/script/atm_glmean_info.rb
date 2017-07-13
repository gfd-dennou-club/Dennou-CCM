#!/usr/env ruby
require "numru/ggraph"
require "optparse"
include NumRu

opt = OptionParser.new
options = {}

opt.on("-l", "--loc_exp_data  <param>",  "the location of experimental directory "){|v| options[:exp_path] = v}
opt.on("-w", "--weight_file <param>",  "the weight for integration"){|v| options[:weight_file] = v}
opt.on("-v", "--variable <param>",  "the name of variable"){|v| options[:var_name] = v}
opt.parse(ARGV)

TargetDir = (options[:exp_path] == nil) ? Dir.pwd : options[:exp_path]
IntegWeightNCFile = (options[:weight_file] == nil) ? "U.nc" : options[:weight_file]
VarName = (options[:var_name] == nil) ? "U" : options[:var_name]

p "weight file: #{IntegWeightNCFile}"
p "Target dir: #{TargetDir}"
p "Target NetCDF data: #{VarName}.nc@#{VarName}"

weight_hash = Hash.new
weight_hash["lon"] = GPhys::IO.open(IntegWeightNCFile, "lon_weight")
weight_hash["lat"] = GPhys::IO.open(IntegWeightNCFile, "lat_weight")
weight_hash["sig"] = sig_weight = GPhys::IO.open(IntegWeightNCFile, "sig_weight")

gp_var = GPhys::IO.open("#{TargetDir}/#{VarName}.nc", VarName)

gp = gp_var
int_wt = UNumeric[1.0, "1"]
for axname in gp_var.axnames
  p axname
  wt = weight_hash[axname]
  gp = (gp*wt).mean(axname)
  p wt.sum(axname)
  int_wt = int_wt * wt.sum(axname)
end

p gp / int_wt





