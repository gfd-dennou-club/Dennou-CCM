#!/usr/env ruby
# coding: utf-8
require "numru/ggraph"
require "optparse"
require "fileutils"
include NumRu

opt = OptionParser.new
options = {}
opt.on("-l", "--loc_exp_data  <param>",  "the location of experimental directory "){|v| options[:exp_path] = v}
opt.on("-d", "--dist_dir  <param>",  "the distination for output data "){|v| options[:dist_path] = v}
opt.on("-o", "--output_fig", "flag to set whether the figures are output"){|v| options[:flag_output]=true}
opt.on("-p", "--prefix_imgfile <param>", "prefix"){|v| options[:imgfile_prefix]=v}
opt.on("-c", "--climate_state <param>", "climate state"){|v| options[:climate_state]=v}
#opt.on("-v", "--var_name <param>",  "variable name of surface temperature"){|v| options[:var_name] = v}
opt.parse(ARGV)

CLIMATE_SNOWBALL = "SnowBallIce"
CLIMATE_PARTIALICE = "PartialIce"
CLIMATE_RUNAWAY = "Runaway"
CLIMATE_WARM    = "Warm"
CLIMATE_PARTICE_COLD = "PartIceCold"
CLIMATE_PARTICE_LARGE = "PartIceLarge"

TargetDir = (options[:exp_path] == nil) ? Dir.pwd : options[:exp_path]
DistDir = (options[:dist_path] == nil) ? Dir.pwd : options[:dist_path]
FlagOutputIMG = (options[:flag_output] == nil) ? false : options[:flag_output]
PrefixOutputIMG = (options[:imgfile_prefix] == nil) ? "" : options[:imgfile_prefix]
ClimateState = (options[:climate_state] == nil) ? CLIMATE_PARTIALICE : options[:climate_state]
PlanetName = "Earth"

LatentHeat = UNumeric[2.4253e6, "J/kg"]
E3 = UNumeric[5.2e3, "m"]
TotDepth = UNumeric[5.2e3, "m"]

p "TargetDir: "+TargetDir
p "DistDir:"+DistDir
p "Is output:"+FlagOutputIMG.to_s
p "Prefix:"+PrefixOutputIMG if FlagOutputIMG
p "ClimateState:"+ClimateState

#------------------------------------------

require File.expand_path(File.dirname(__FILE__) + "/../../../tool/postproc/common/ConstUtil.rb")
require File.expand_path(File.dirname(__FILE__) + "/../../../tool/postproc/common/DCModelIOUtil.rb")
require File.expand_path(File.dirname(__FILE__) + "/../../../tool/postproc/atm_analysis/DCPAMUtil.rb")
eval("include ConstUtil::#{PlanetName}")
include DCModelIOUtil
include NumRu

#------------------------------------------

DCL::swlset('lwnd', false) if FlagOutputIMG


def get_gp(comp, varname)
  targetDir = TargetDir.sub("/?/", "/#{comp}/")
  p "get_gp: #{targetDir}/#{varname}.nc@#{varname} .." 
  return GPhys::IO.open(targetDir + "/" + varname+".nc", varname)
end

def prep_dcl(iws, itr=1, clrmap=10)
#  DCLExt.sg_set_params('ifont'=>2)
  DCL.sgscmn(clrmap)           
  DCL.gropn(iws)              
  DCL.sgpset('isub', 96)     # 下付き添字を表す制御文字を '_' から '`' に
  DCL.glpset('lmiss',true)   # DCLの欠損値処理を on に．

  GGraph.set_fig('itr'=>itr) 
  GGraph.set_fig('viewport'=>[0.13,0.87,0.22,0.79])  
  GGraph.set_axes('ytitle'=>'sigma')
end

def rename_pngfile(fbasename)
  File.rename(DistDir+"/dcl_0001.png", DistDir+"/#{PrefixOutputIMG}#{fbasename}.png")
end

#-------------------

def meridional_heat_flux_fig(ocnHT, atmHT, itr=1)
  p "meridional_heat_flux_fig.."

  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    htflx_min = -1.0; htflx_max = 1.0
  when CLIMATE_PARTICE_COLD, CLIMATE_PARTICE_LARGE then
    htflx_min = -2.0; htflx_max = 2.0    
  when CLIMATE_WARM then
    htflx_min = -5.0; htflx_max = 5.0
  when CLIMATE_RUNAWAY
    htflx_min = -20.0; htflx_max = 20.0    
  else
    htflx_min = -3.0; htflx_max = 3.0
  end

  atmHT_ = atmHT.interpolate(ocnHT.axis("lat").pos)
  totHT = atmHT_ + ocnHT
  GGraph.set_axes("ytitle"=>"meridional heat flux", "yunits"=>"PW")
  GGraph.line(totHT, true, "titl"=>"meridional heat flux", "min"=>htflx_min, "max"=>htflx_max, "index"=>13,
              "legend"=>"tot", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(atmHT_, false,  "index"=>23,
              "legend"=>"atm", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(ocnHT, false, "index"=>33,
              "legend"=>"ocn", "legend_dx"=>0.024, "legend_size"=>0.018)
  DCL.grcls
  
  rename_pngfile("HeatFluxLat") if FlagOutputIMG   
end

atmHT = get_gp("atm", "moistStatEnFlxLat")
ocnHT = get_gp("ocn", "TotHT")

meridional_heat_flux_fig( ocnHT, atmHT )

