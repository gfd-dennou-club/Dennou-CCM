# coding: utf-8
#!/usr/env ruby

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

opt.parse(ARGV)

CLIMATE_SNOWBALL = "SnowBallIce"
CLIMATE_PARTIALICE = "PartialIce"
CLIMATE_RUNAWAY = "Runaway"
CLIMATE_WARM = "Warm"
CLIMATE_PARTICE_COLD = "PartIceCold"
CLIMATE_PARTICE_LARGE = "PartIceLarge"

TargetDir = (options[:exp_path] == nil) ? Dir.pwd : options[:exp_path]
DistDir = (options[:dist_path] == nil) ? Dir.pwd : options[:dist_path]
FlagOutputIMG = (options[:flag_output] == nil) ? false : options[:flag_output]
PrefixOutputIMG = (options[:imgfile_prefix] == nil) ? "" : options[:imgfile_prefix]
ClimateState = (options[:climate_state] == nil) ? CLIMATE_PARTIALICE : options[:climate_state]
PlanetName = "Earth"

LatentHeat = UNumeric[2.4253e6, "J/kg"]

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

def get_gp(ncname, varname=ncname.sub(".nc",""))
  return GPhys::IO.open(TargetDir + "/" + ncname, varname)
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

def prep_dcl_tserise(iws, itr=1, clrmap=1, row=1, col=1)

  DCL.swpset("iwidth", 650)
  DCL.swpset("iheight", 550)

  DCL.sgscmn(clrmap)           
  DCL.gropn(iws)              
  DCL.sgpset('isub', 96)     # 下付き添字を表す制御文字を '_' から '`' に
  DCL.glpset('lmiss',true)   # DCLの欠損値処理を on に．

  GGraph.set_fig('itr'=>itr) 

  DCL.sldiv('t', row, col)
  DCL.sgpset("lfull", true)
  DCL.uzfact(0.75)
  
end

def rename_pngfile(fbasename)
  File.rename(DistDir+"/dcl_0001.png", DistDir+"/#{PrefixOutputIMG}#{fbasename}.png")
end

def sicethick_siceEn_glmean_tserise_fig(gmSIceThick, gmSIceEn, itr=1)

  case ClimateState
  when CLIMATE_SNOWBALL then
    sicethick_min = 0; sicethick_max = 400
    siceen_min = -5e11; siceen_max = 0.0
  when CLIMATE_PARTICE_LARGE then
    sicethick_min = 0; sicethick_max = 80
    siceen_min = -5e10; siceen_max = 0.0
  when CLIMATE_PARTICE_COLD then
    sicethick_min = 0; sicethick_max = 55
    siceen_min = -2e10; siceen_max = 0.0
  when CLIMATE_WARM then
    sicethick_min = 0; sicethick_max = 5
    siceen_min = -1e9; siceen_max = 0.0
  when CLIMATE_RUNAWAY then
    sicethick_min = 0; sicethick_max = 10
    siceen_min = -5e8; siceen_max = 0.0
  else
    sicethick_min = 0; sicethick_max = 10
    siceen_min = -3e9; siceen_max = 0.0
  end

  prep_dcl_tserise(1, itr, 1, 1, 2)

  GGraph.set_fig('viewport'=>[0.15,0.9,0.05,0.34])
  GGraph.set_axes("ytitle"=>"ice thickness", "yunits"=>"m")
  GGraph.line(gmSIceThick, true, {  "title"=>"time serise of global-mean field", "max"=>sicethick_max, "min"=>sicethick_min,
                                "index"=>13})

  GGraph.set_fig('viewport'=>[0.15,0.9,0.1,0.38])      
  GGraph.set_axes("ytitle"=>"sea ice energy", "yunits"=>"J")
  GGraph.line(gmSIceEn, true, {  "title"=>"", "max"=>siceen_max, "min"=>siceen_min,
                               "index"=>13})
  DCL.grcls
  
  rename_pngfile("SIceThickSIceEnGlMean_tserise") if FlagOutputIMG       
end

#---

icethick = get_gp("IceThickGlMean-standalone.nc", "IceThick")
siceEn = get_gp("SIceEnGlMean-standalone.nc", "SIceEn")

#----

sicethick_siceEn_glmean_tserise_fig(icethick, siceEn)
