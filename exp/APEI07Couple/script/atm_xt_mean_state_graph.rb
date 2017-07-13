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


def get_gp(varname)
  p "get_gp: #{TargetDir}/#{varname}.nc@#{varname} .." 
  return GPhys::IO.open(TargetDir + "/" + varname+".nc", varname)
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

def u_temp_fig(u, temp, itr)
  p "u_temp_fig: itr=#{itr}"
  prep_dcl(1,itr,10)

  case  ClimateState
  when CLIMATE_RUNAWAY
    temp_min = 110; temp_max = 450; temp_int = 10
  else
    temp_min = 60; temp_max = 400; temp_int = 10
  end
  
  u_min = -120; u_max = 120; u_int = 10
  
  GGraph.tone( temp, true, "titl"=>"U, T", "int"=>temp_int, "min"=>temp_min, "max"=>temp_max )
  GGraph.contour( u, false, "titl"=>"U, T", "int"=>u_int, "min"=>u_min, "max"=>u_max )
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.25)
  DCL.grcls

  rename_pngfile("U-T_xtmean_itr#{itr.to_i}") if FlagOutputIMG
end

def msf_qvap_fig(msf, qvap, itr)
  p "msf_qvap_fig: itr=#{itr}"
  prep_dcl(1, itr, 28)

  qvap_min = 1e-9; qvap_max = 1.0; qvap_int = 1e-2
  msf_min = -120; msf_max = 120; msf_int = 10

  case ClimateState
  when CLIMATE_SNOWBALL then
    qvap_min = 1e-10; qvap_max = 5e-3
  when CLIMATE_WARM then
    qvap_max = 1e-1
  when CLIMATE_RUNAWAY
    qvap_min = 1e-7; qvap_max = 2.0
  else
    qvap_max = 1e-1
  end

  GGraph.tone( qvap, true, "titl"=>"MSF, QVap", "int"=>qvap_int, "log"=>true, "log_cycle"=>3,
               "min"=>qvap_min, "max"=>qvap_max, "clr_min"=>10, "clr_max"=>99, "annotate"=>false )
  sig_pos = qvap.axis("sig").pos.rename("sigm")
  GGraph.contour( msf.interpolate(sig_pos), false, "titl"=>"MSF, QVap", "int"=>msf_int, "min"=>msf_min, "max"=>msf_max )
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.25, "log"=>true)
  DCL.grcls

  rename_pngfile("MSF-QH2OVap_xtmean_itr#{itr.to_i}") if FlagOutputIMG
end

def energy_flux_fig(olr, osr, slr, ssr, rain, evap, sens, itr=1)
  p "energy_flux_fig.."
  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    energyflx_min = -50; energyflx_max = 250;
  when CLIMATE_RUNAWAY
    energyflx_min = -50; energyflx_max = 600;
  else
    energyflx_min = -50; energyflx_max = 500;    
  end

  GGraph.set_axes("ytitle"=>"heat flux", "yunits"=>"W/m2")
  GGraph.line(rain, true, "titl"=>"heat flux", "min"=>energyflx_min, "max"=>energyflx_max, "index"=>13,
              "legend"=>"PRCP", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(olr, false, "index"=>23,
              "legend"=>"OLR", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(-osr, false, "index"=>33,
              "legend"=>"-OSR", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(slr, false, "index"=>43,
              "legend"=>"SLR", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(evap, false, "index"=>53,
              "legend"=>"Evap", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(sens, false, "index"=>63,
              "legend"=>"Sens", "legend_dx"=>0.024, "legend_size"=>0.018)
  DCL.grcls
  
  rename_pngfile("EnergyFlux_xtmean") if FlagOutputIMG
end

def wind_stress_x_fig(taux, itr=1)
  p "wind_stress_x_fig.."

  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    taux_min = -0.01; taux_max = 0.02
  else
    taux_min = -0.05; taux_max = 0.1
  end
  
  GGraph.set_axes("ytitle"=>"wind stress (lon)", "yunits"=>"N/m2")
  GGraph.line(taux, true,  "title"=>"wind stress (lon)", "min"=>taux_min, "max"=>taux_max,
              "index"=>13)
  DCL.grcls
  
  rename_pngfile("TauX_xtmean") if FlagOutputIMG  
end

def wind_stress_y_fig(tauy, itr=1)
  p "wind_stress_y_fig.."
  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    tauy_min = -0.01; tauy_max = 0.01
  else
    tauy_min = -0.05; tauy_max = 0.05
  end
  
  GGraph.set_axes("ytitle"=>"wind stress (lat)", "yunits"=>"N/m2")
  GGraph.line(tauy, true, "title"=>"wind stress (lat)", "min"=>tauy_min, "max"=>tauy_max,
             "index"=>13)
  DCL.grcls
  
  rename_pngfile("TauY_xtmean") if FlagOutputIMG  
end

def meridional_heat_flux_fig(dryStatEn, latentEn, moistStatEn, itr=1)
  p "meridional_heat_flux_fig.."

  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    htflx_min = -2.0; htflx_max = 2.0
  when CLIMATE_WARM then
    htflx_min = -11.0; htflx_max = 11.0
  when CLIMATE_RUNAWAY
    htflx_min = -40.0; htflx_max = 40.0    
  else
    htflx_min = -6.0; htflx_max = 6.0
  end

  GGraph.set_axes("ytitle"=>"meridional heat flux", "yunits"=>"PW")
  GGraph.line(moistStatEn, true, "titl"=>"meridional heat flux", "min"=>htflx_min, "max"=>htflx_max, "index"=>13,
              "legend"=>"moist", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(dryStatEn, false,  "index"=>23,
              "legend"=>"dry", "legend_dx"=>0.024, "legend_size"=>0.018)
  GGraph.line(latentEn, false, "index"=>33,
              "legend"=>"latent", "legend_dx"=>0.024, "legend_size"=>0.018)
  DCL.grcls
  
  rename_pngfile("HeatFluxLat") if FlagOutputIMG   
end

def sfctemp_fig(sfctemp, itr=1)
  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    sfctemp_min = 160.0; sfctemp_max = 300.0
  when CLIMATE_RUNAWAY then
    sfctemp_min = 280.0; sfctemp_max = 500.0
  else
    sfctemp_min = 180.0; sfctemp_max = 320.0
  end
  
  GGraph.set_axes("ytitle"=>"surface temperature", "yunits"=>"K")
  GGraph.line(sfctemp, true,  "title"=>"surface temperature", "min"=>sfctemp_min, "max"=>sfctemp_max,
              "index"=>13)
  DCL.grcls
  
  rename_pngfile("SfcTemp_xtmean") if FlagOutputIMG
end

def prcp_fig(prcp, itr=1)
  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    prcp_min = 0.0; prcp_max = 500.0
  when CLIMATE_PARTICE_COLD then
    prcp_min = 0.0; prcp_max = 800.0
  when CLIMATE_WARM then
    prcp_min = 0.0; prcp_max = 8000.0
  when CLIMATE_RUNAWAY then
    prcp_min = 0.0; prcp_max = 8000.0
  else
    prcp_min = 0.0; prcp_max = 4000.0
  end


  prcp_mm_yr = prcp * (UNumeric[86400*365.0,'s'] * UNumeric[1e3, "mm/m"] / UNumeric[1e3, "kg/m3"])
  GGraph.set_axes("ytitle"=>"precipitation", "yunits"=>"mm/year")
  GGraph.line(prcp_mm_yr, true, 
              "title"=>"precipitation", "min"=>prcp_min, "max"=>prcp_max,
              "index"=>13)
  DCL.grcls
  
  rename_pngfile("PRCP_xtmean") if FlagOutputIMG  
end

#------------------------------------------

u = get_gp("U")
temp = get_gp("Temp")

msf = get_gp("MSF")
qvap = get_gp("QH2OVap")

olr = get_gp("OLRA")
slr = get_gp("SLRA")
osr = get_gp("OSRA")
ssr = get_gp("SSRA")

evap = get_gp("LatHFlxA")
sens = get_gp("SenHFlxA")

taux = get_gp("TauX")
tauy = get_gp("TauY")

prcp = get_gp("PRCP")
sfctemp = get_gp("SfcTemp")

dryEnFlx = get_gp("dryStatEnFlxLat")
latentEnFlx = get_gp("latentEnFlxLat")
moistEnFlx = get_gp("moistStatEnFlxLat")

#------------------------------------------

u_temp_fig( u, temp, 1)
u_temp_fig( u, temp, 2)

msf_qvap_fig( msf, qvap, 1)
msf_qvap_fig( msf, qvap, 2)

energy_flux_fig(olr, osr, slr, ssr, prcp*LatentHeat, evap, sens)

wind_stress_x_fig(-taux)
wind_stress_y_fig(-tauy)

meridional_heat_flux_fig(dryEnFlx, latentEnFlx, moistEnFlx)

prcp_fig(prcp)
sfctemp_fig(sfctemp)
