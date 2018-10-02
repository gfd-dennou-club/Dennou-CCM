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

opt.parse(ARGV)

CLIMATE_SNOWBALL = "SnowBallIce"
CLIMATE_PARTIALICE = "PartialIce"
CLIMATE_RUNAWAY = "Runaway"
CLIMATE_WARM = "Warm"

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

def rename_pngfile(fbasename)
  File.rename(DistDir+"/dcl_0001.png", DistDir+"/#{PrefixOutputIMG}#{fbasename}.png")
end

def iceline_lat_tserise_fig(iceline, itr=1)
  prep_dcl(1, itr, 8)
  GGraph.set_fig('viewport'=>[0.1,0.9,0.22,0.79])  

  GGraph.set_axes("ytitle"=>"ice-line latitude", "yunits"=>"degrees")
  GGraph.line(iceline, true, {"titl"=>"time serise of iceline latitude", "max"=>95.0, "min"=>-5.0,
                              "index"=>13})
  DCL.grcls
  
  rename_pngfile("IcelineLat_tserise") if FlagOutputIMG       
end

def sfctemp_glmean_tserise_fig(gmSfcTemp, itr=1)
  prep_dcl(1, itr, 8)
  GGraph.set_fig('viewport'=>[0.1,0.9,0.22,0.79])  

  case ClimateState
  when CLIMATE_SNOWBALL then
    sfctemp_min = 180; sfctemp_max = 320
  when CLIMATE_RUNAWAY then
    sfctemp_min = 250; sfctemp_max = 430
  else
    sfctemp_min = 180; sfctemp_max = 320
  end
  
  GGraph.set_axes("ytitle"=>"global-mean surface temperature", "yunits"=>"K")
  GGraph.line(gmSfcTemp, true, {"titl"=>"time serise of global-mean surface temperature",
                                "max"=>sfctemp_max, "min"=>sfctemp_min,
                                "index"=>13})
  DCL.grcls
  
  rename_pngfile("SfcTempGlMean_tserise") if FlagOutputIMG       
end

def energy_glmean_tserise_fig(totEn, intEn, potEn, latEn, kinEn, itr=1)
  prep_dcl(1, itr, 8)
  GGraph.set_fig('viewport'=>[0.1,0.9,0.22,0.79])  

  case ClimateState
  when CLIMATE_RUNAWAY then
    normalized_en_max = 3.5; latEnFac=1.0; 
  else
    normalized_en_max = 1.5; latEnFac=10.0
  end

  tot_En_mean = totEn.mean("time")
  
  legend_common = {"legend_dx"=>0.024, "legend_vx"=> 0.65, "legend_size"=>0.018}
  GGraph.set_axes("ytitle"=>"normalized global-mean energy", "yunits"=>"1")
  GGraph.line( totEn/tot_En_mean, true, {"titl"=>"time serise of normalized energy", "max"=>normalized_en_max, "min"=>0,
                                         "index"=>13, "legend"=>"total", "legend_vy"=> 0.75}.merge(legend_common) )
  GGraph.line( intEn/tot_En_mean, false,              
               {"index"=>23, "legend"=>"internal"}.merge(legend_common) )

  GGraph.line( potEn/tot_En_mean, false,              
               {"index"=>33, "legend"=>"potential"}.merge(legend_common) )

  GGraph.line( latEn/tot_En_mean*latEnFac, false,              
               {"index"=>43, "legend"=>"latent(#{latEnFac.to_i} times)"}.merge(legend_common) )
  
  GGraph.line( kinEn/tot_En_mean*5e2, false,              
               {"index"=>53, "legend"=>"kinetic(500 times)"}.merge(legend_common) )
  DCL.grcls
  
  rename_pngfile("EngyGlMean_tserise") if FlagOutputIMG     
end

def taxis_unitconv(gp)
  taxis = gp.axis("time")
  if taxis.pos.units.to_s == "day" then
    taxis_new = Axis.new.set_pos(taxis.pos/UNumeric[365, "day/year"])
    grid_new = gp.copy.grid.change_axis(gp.dim_index("time"), taxis_new)
    gp = GPhys.new(grid_new, gp.data)
  end
  return gp
end

#---

totEn = taxis_unitconv( get_gp("TotEngyGlMean.nc", "TotEngy") )
intEn = taxis_unitconv( get_gp("IntEngyGlMean.nc", "IntEngy") )
latEn = taxis_unitconv( get_gp("LatEngyGlMean.nc", "LatEngy") )
potEn = taxis_unitconv( get_gp("PotEngyGlMean.nc", "PotEngy") )
kinEn = taxis_unitconv( get_gp("KinEngyGlMean.nc", "KinEngy") )

iceline = get_gp("iceline_lat.nc", "iceline_lat")
gmSfcTemp = get_gp("sfctemp_glmean.nc", "SfcTemp")

#----

energy_glmean_tserise_fig(totEn, intEn, potEn, latEn, kinEn)
sfctemp_glmean_tserise_fig(gmSfcTemp)
iceline_lat_tserise_fig(iceline)

p totEn
