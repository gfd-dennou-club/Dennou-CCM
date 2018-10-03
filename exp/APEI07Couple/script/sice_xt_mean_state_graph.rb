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


def get_gp(varname)
  p "get_gp: #{TargetDir}/#{varname}.nc@#{varname} .." 
  return GPhys::IO.open(TargetDir + "/" + varname+".nc", varname)
end

def prep_dcl(iws, itr=1, clrmap=10)
#  DCLExt.sg_set_params('ifont'=>2)
  DCL.sgscmn(clrmap)           
  DCL.sgpset('isub', 96)     # 下付き添字を表す制御文字を '_' から '`' に
  DCL.glpset('LMISS',true)   # DCLの欠損値処理を on に．
  DCL.gropn(iws)              

  GGraph.set_fig('itr'=>itr) 
  GGraph.set_fig('viewport'=>[0.13,0.87,0.22,0.79])  
  GGraph.set_axes('ytitle'=>'sigma')
end

def prep_dcl_split(iws, itr=1, clrmap=1, row=1, col=1)

  DCL.swpset("iwidth", 650)
  DCL.swpset("iheight", 550)

  DCL.sgscmn(clrmap)           
  DCL.sgpset('isub', 96)     # 下付き添字を表す制御文字を '_' から '`' に
  DCL.glpset('LMISS',true)   # DCLの欠損値処理を on に．
  DCL.gropn(iws)              

  GGraph.set_fig('itr'=>itr) 

  DCL.sldiv('t', row, col)
  DCL.sgpset("lfull", true)
  DCL.uzfact(0.75)
  
end

def rename_pngfile(fbasename)
  File.rename(DistDir+"/dcl_0001.png", DistDir+"/#{PrefixOutputIMG}#{fbasename}.png")
end

def meridional_heat_flux_fig(totHT, eulerHT, isoDiffHT, bolusHT, itr=1)
  p "meridional_heat_flux_fig.."

  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    htflx_min = -1.0; htflx_max = 1.0
  when CLIMATE_PARTICE_LARGE then
    htflx_min = -2.0; htflx_max = 2.0    
  when CLIMATE_PARTICE_COLD then
    htflx_min = -2.0; htflx_max = 2.0    
  when CLIMATE_WARM then
    htflx_min = -5.0; htflx_max = 5.0
  when CLIMATE_RUNAWAY
    htflx_min = -20.0; htflx_max = 20.0    
  else
    htflx_min = -3.0; htflx_max = 3.0
  end

  GGraph.set_axes("ytitle"=>"meridional heat flux", "yunits"=>"PW")
  GGraph.line(totHT, true, "titl"=>"meridional heat flux", "min"=>htflx_min, "max"=>htflx_max, "index"=>13,
              "legend"=>"tot", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(eulerHT, false,  "index"=>23,
              "legend"=>"euler", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(bolusHT, false, "index"=>33,
              "legend"=>"bolus", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(isoDiffHT, false, "index"=>43,
              "legend"=>"isopyc", "legend_dx"=>0.024, "legend_size"=>0.018)
  DCL.grcls
  
  rename_pngfile("HeatFluxLat") if FlagOutputIMG   
end

def sicetemp_fig(sfctemp_, icetemp_, itr=1)
  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    sfctemp_min = -130; sfctemp_max = 0.0
    icetemp_min = -50.0; icetemp_max = 0.0
  when CLIMATE_RUNAWAY then
    sfctemp_min = -10; sfctemp_max = 0.0
    icetemp_min = -10.0; icetemp_max = 0.0
  when CLIMATE_WARM then
    sfctemp_min = -70; sfctemp_max = 0.0
    icetemp_min = -10.0; icetemp_max = 0.0
  when CLIMATE_PARTICE_COLD then
    sfctemp_min = -100; sfctemp_max = 0.0
    icetemp_min = -30.0; icetemp_max = 0.0
  when CLIMATE_PARTICE_LARGE then
    sfctemp_min = -100; sfctemp_max = 0.0
    icetemp_min = -30.0; icetemp_max = 0.0    
  else
    sfctemp_min = -80; sfctemp_max = 0.0
    icetemp_min = -20.0; icetemp_max = 0.0
  end

  rm = NArray.float(1).fill(DCL::glrget('RMISS'))
  sfctemp = GPhys.new(sfctemp_.grid, VArray.new(sfctemp_.val,
                                                {'missing_value'=>rm, 'units'=>'degC'}, sfctemp_.name))
  icetemp = GPhys.new(icetemp_.grid, VArray.new(icetemp_.val,
                                                {'missing_value'=>rm, 'units'=>'degC'}, icetemp_.name))

#  ocn_mask = (sfctemp < -1e10).where
#  sfctemp[ocn_mask] = rm
#  ocn_mask = (icetemp < -1e10).where
#  icetemp[ocn_mask] = rm
  #  p sfctemp.val[10..20]

  
  GGraph.set_axes("ytitle"=>"temperature", "yunits"=>"degC")
  GGraph.line(sfctemp, true,  "title"=>"sea-ice temperature", "min"=>sfctemp_min, "max"=>sfctemp_max, "index"=>13,
              "legend"=>"surface", "legend_vx"=>-0.2, "legend_vy"=>-0.45, "legend_dx"=>0.025, "legend_size"=>0.018 )
  GGraph.line(icetemp.cut('sig2'=>-0.25), false,  "index"=>23,
              "legend"=>"upper layer", "legend_vx"=>-0.2, "legend_dx"=>0.025, "legend_size"=>0.018)
  GGraph.line(icetemp.cut('sig2'=>-0.75), false,  "index"=>33, 
              "legend"=>"lower layer", "legend_vx"=>-0.2, "legend_dx"=>0.025, "legend_size"=>0.018)  
  DCL.grcls
  
  rename_pngfile("SIceTemp_xtmean") if FlagOutputIMG
end

def sicethick_fig(ice, snow, itr=1)
  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    thick_min = 0; thick_max = 600
  when CLIMATE_RUNAWAY then
    thick_min = 0; thick_max = 3
  when CLIMATE_WARM then
    thick_min = 0; thick_max = 20
  when CLIMATE_PARTICE_COLD then
    thick_min = 0; thick_max = 100
  when CLIMATE_PARTICE_LARGE then
    thick_min = 0; thick_max = 300
  else
    thick_min = 0; thick_max = 30
  end
  
  GGraph.set_axes("ytitle"=>"thickness", "yunits"=>"m")
  GGraph.line(ice, true,  "title"=>"sea-ice thickness", "min"=>thick_min, "max"=>thick_max, "index"=>13, "type"=>3,
              "legend"=>"ice", "legend_vx"=>-0.18, "legend_dx"=>0.025, "legend_size"=>0.018)
  GGraph.line(ice+snow, false,  "index"=>13, "type"=>1, 
              "legend"=>"ice+snow", "legend_vx"=>-0.18, "legend_dx"=>0.025, "legend_size"=>0.018)

  DCL.grcls
  
  rename_pngfile("SIceThick_xtmean") if FlagOutputIMG
end

def sfchflx_fig(sfchflx, itr=1)
  prep_dcl(1, itr, 8)
  
  case ClimateState
  when CLIMATE_SNOWBALL then
    sfchflx_min = -30.0; sfchflx_max = 30.0
  when CLIMATE_RUNAWAY then
    sfchflx_min = -30.0; sfchflx_max = 30.0
  else
    sfchflx_min = -50.0; sfchflx_max = 25.0
  end
  
  GGraph.set_axes("ytitle"=>"surface flux", "yunits"=>"K")
  GGraph.line(sfchflx, true,  "title"=>"surface heat flux(ocn<->atm/sice)",
              "min"=>sfchflx_min, "max"=>sfchflx_max,
              "index"=>13)
  DCL.grcls

  rename_pngfile("SfcHFlxO_xtmean") if FlagOutputIMG
end

#------------------------------------------

icethick = get_gp("IceThick")
snowthick = get_gp("SnowThick")                  
sfctemp = get_gp("SIceSfcTemp")
sicetemp = get_gp("SIceTemp")

#------------------------------------------

sicethick_fig( icethick, snowthick, 1 )
sicetemp_fig(sfctemp, sicetemp, 1)
