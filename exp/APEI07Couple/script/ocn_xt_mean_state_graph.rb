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

def z_to_figz(z, upperLyrDepth)
    return (z > -upperLyrDepth) ? z :
                - upperLyrDepth + upperLyrDepth*(upperLyrDepth + z)/(TotDepth - upperLyrDepth)
end

def sig_to_figZgrid(grid, itr)
  gp_sig = grid.axis('sig').to_gphys
  nz = gp_sig.shape[0]
  na_z = NArray.sfloat(nz)

  upperLyrDepth = UNumeric[1e3,"m"]
  for k in 0..nz-1
    z = E3*gp_sig.val[k]
    na_z[k] = (itr == 2)? z_to_figz(z, upperLyrDepth) : z
#    p na_z[k]
  end

  ry1 = (0..26).map{|i|
    z = -i*200.0
    (itr == 2)? z_to_figz(z, upperLyrDepth) : z
  }
  
  cy = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0', '2.0', '3.0', '4.0', '5.0']
  ry2 = cy.map{|z|
    z_ = -z.to_f*1e3
    (itr == 2)? z_to_figz(z_, upperLyrDepth) : z_
  }
  va_z = VArray.new(na_z, {"long_name"=>"ocean depth", "units"=>"m"}, "depth")
  return grid.change_axis(1, Axis.new.set_pos(va_z)), ry1, ry2, cy
end

def u_ptemp_fig(u_, ptemp_, itr)
  p "u_ptemp_fig: itr=#{itr}"
  prep_dcl(1,1,40)

  u_min = -0.8; u_max = 0.3; u_int = 0.025
  
  case  ClimateState
  when CLIMATE_SNOWBALL then
    temp_min = 271; temp_max = 280; temp_int = 0.25
  when CLIMATE_PARTICE_LARGE  then
    temp_min = 271; temp_max = 290; temp_int = 0.5
  when CLIMATE_PARTICE_COLD then
    temp_min = 271; temp_max = 297; temp_int = 0.5
  when CLIMATE_WARM then
    temp_min = 271; temp_max = 313; temp_int = 2
  when CLIMATE_RUNAWAY
    temp_min = 280; temp_max = 360; temp_int = 5
    u_min = -1.0; u_max = 1.0; u_int = 0.05    
  else
    temp_min = 271; temp_max = 304; temp_int = 1
  end
  


  fig_z_grid, ry1, ry2, cy = sig_to_figZgrid(u_.grid,itr)  
  u = GPhys.new(fig_z_grid, u_.data)  
  ptemp = GPhys.new(fig_z_grid, ptemp_.data)
  
  label = "U, " + DCL::csgi(135)
  GGraph.next_axes('yside'=>'u')
  GGraph.tone( ptemp, true, "titl"=>label, "int"=>temp_int, "min"=>temp_min, "max"=>temp_max )
  GGraph.contour( u, false, "titl"=>label, "int"=>u_int, "min"=>u_min, "max"=>u_max )
  DCL::uyaxlb('L', ry1, ry2, cy, 4)
  DCL::uyaxlb('R', ry1, ry2, cy, 4)
  DCL::uysttl('L', 'depth (km)', 0.0)
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.5)

  DCL.grcls
  

  rename_pngfile("U-PTemp_xtmean_itr#{itr.to_i}") if FlagOutputIMG
end

def gen_levels(min_val, max_val, interval, is_exceed_val=false)
  nlevels = ((max_val - min_val)/interval).to_i + 1
  levels = (0..nlevels-1).map{|i|
    min_val + i*interval
  }
  if is_exceed_val then
    rmiss = DCL.glpget('rmiss')
    levels.unshift(rmiss); levels.push(rmiss)
  end
  return levels
end

def msf_salt_fig(msf_, salt_, itr)
  p "msf_salt_fig: itr=#{itr}"
  prep_dcl(1, 1, 1)


  salt_levels = []
  case ClimateState
  when CLIMATE_SNOWBALL then
    salt_levels = gen_levels(36.0, 37.5, 0.25) 
    msf_min = -50; msf_max = 50; msf_int = 5
  when CLIMATE_PARTICE_LARGE then
    salt_levels = [33.0].concat( gen_levels(33.4, 36.4, 0.2) )
    msf_min = -50; msf_max = 50; msf_int = 5    
  when CLIMATE_PARTICE_COLD then
    salt_levels = [33.0].concat( gen_levels(33.4, 36.4, 0.2) )
    msf_min = -50; msf_max = 50; msf_int = 5
  when CLIMATE_WARM then
    salt_levels = [33.0].concat( gen_levels(33.5, 37.25, 0.25).concat([37.5,38.0,39.0,40.0]) )
    msf_min = -50; msf_max = 50; msf_int = 5
  when CLIMATE_RUNAWAY
    salt_levels = [20.0,30.0,33.0].concat( gen_levels(33.5, 37.25, 0.25).concat([37.5,38.0,39.0]) )
    msf_min = -50; msf_max = 50; msf_int = 5
  else
    salt_levels = [33.0].concat( gen_levels(33.4, 36.8, 0.2).concat([37.2, 38.0]) )
    msf_min = -50; msf_max = 50; msf_int = 5
  end

  p salt_levels
  
  sig_pos = salt_.axis("sig").pos.rename("sigm")  
  fig_z_grid, ry1, ry2, cy = sig_to_figZgrid(salt_.grid,itr)  
  salt = GPhys.new(fig_z_grid, salt_.data)  
  msf = GPhys.new(fig_z_grid, msf_.interpolate(sig_pos).data)

  GGraph.next_axes('yside'=>'u')  
  GGraph.tone( salt, true, "titl"=>"MSF, Salt", "levels"=>salt_levels)
  GGraph.contour( msf, false, "titl"=>"MSF, Salt", "int"=>msf_int, "min"=>msf_min, "max"=>msf_max )
  DCL::uyaxlb('L', ry1, ry2, cy, 4)
  DCL::uyaxlb('R', ry1, ry2, cy, 4)
  DCL::uysttl('L', 'depth (km)', 0.0)
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.5)
  DCL.grcls

  rename_pngfile("MSF-Salt_xtmean_itr#{itr.to_i}") if FlagOutputIMG
end

def stratification_fig(densPot_, bvFreq_, itr=1)
  p "densPot_bvFreq_fig: itr=#{itr}"
  prep_dcl(1, 1, 1)


  case ClimateState
  when CLIMATE_SNOWBALL then
  when CLIMATE_PARTICE_COLD, CLIMATE_PARTICE_LARGE then
    bvFreq_min = -1e-5; bvFreq_max = 2e-4; bvFreq_int = 1e-5
    bvFreq_levels = [-2e-6, -8e-7, -5e-8, -2e-8, 0.0, 2e-7, 5e-7, 8e-7, 2e-6, 5e-6, 8e-6, 2e-5, 5e-5, 8e-5, 2e-4, 5e-4, 8e-4]
    densPot_min = -4; densPot_max = 2; densPot_int = 0.2
  when CLIMATE_WARM then
    bvFreq_min = -1e-5; bvFreq_max = 2e-4; bvFreq_int = 1e-5
    bvFreq_levels = [-2e-6, -8e-7, -5e-8, -2e-8, 0.0, 2e-7, 5e-7, 8e-7, 2e-6, 5e-6, 8e-6, 2e-5, 5e-5, 8e-5, 2e-4, 5e-4, 8e-4]
    densPot_min = -4; densPot_max = 2; densPot_int = 0.2
  when CLIMATE_RUNAWAY
  else
    bvFreq_min = -1e-5; bvFreq_max = 2e-4; bvFreq_int = 1e-5
    bvFreq_levels = [-2e-6, -8e-7, -5e-8, -2e-8, 0.0, 2e-7, 5e-7, 8e-7, 2e-6, 5e-6, 8e-6, 2e-5, 5e-5, 8e-5, 2e-4, 5e-4, 8e-4]
    densPot_min = -4; densPot_max = 2; densPot_int = 0.2
  end

  fig_z_grid, ry1, ry2, cy = sig_to_figZgrid(bvFreq_.grid,itr)
  bvFreq = GPhys.new(fig_z_grid, bvFreq_.data)
  densPot = GPhys.new(fig_z_grid, densPot_.data) 

  GGraph.next_axes('yside'=>'u')  
  GGraph.tone( bvFreq, true, "titl"=>"DensPot, BvFreq", "levels"=>bvFreq_levels )#"int"=>bvFreq_int, "min"=>bvFreq_min, "max"=>bvFreq_max )
  GGraph.contour( densPot, false, "titl"=>"DensPot, BvFreq", "int"=>densPot_int, "min"=>densPot_min, "max"=>densPot_max )
  DCL::uyaxlb('L', ry1, ry2, cy, 4)
  DCL::uyaxlb('R', ry1, ry2, cy, 4)
  DCL::uysttl('L', 'depth (km)', 0.0)
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.5, "constwidth"=>true, "chval_fmt"=>"B")
  DCL.grcls

  rename_pngfile("DensPot-BvFreq_xtmean_itr#{itr.to_i}") if FlagOutputIMG
end

def meridional_heat_flux_fig(totHT, eulerHT, isoDiffHT, bolusHT, itr=1)
  p "meridional_heat_flux_fig.."

  prep_dcl(1, itr, 8)

  case ClimateState
  when CLIMATE_SNOWBALL then
    htflx_min = -1.0; htflx_max = 1.0
  when CLIMATE_PARTICE_COLD, CLIMATE_PARTICE_LARGE then
    htflx_min = -2.0; htflx_max = 2.0    
  when CLIMATE_WARM then
    htflx_min = -4.0; htflx_max = 4.0
  when CLIMATE_RUNAWAY
    htflx_min = -20.0; htflx_max = 20.0    
  else
    htflx_min = -2.0; htflx_max = 2.0
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

def sfchflx_fig(sfchflx, itr=1)
  prep_dcl(1, itr, 8)
  
  case ClimateState
  when CLIMATE_SNOWBALL then
    sfchflx_min = -30.0; sfchflx_max = 30.0
  when CLIMATE_RUNAWAY then
    sfchflx_min = -30.0; sfchflx_max = 30.0
  when CLIMATE_WARM then
    sfchflx_min = -60.0; sfchflx_max = 40.0
  else
    sfchflx_min = -50.0; sfchflx_max = 25.0
  end

  GGraph.set_axes("ytitle"=>"surface heat flux", "yunits"=>"K")
  GGraph.line(sfchflx, true,  "title"=>"surface heat flux(ocn<->atm/sice)",
              "min"=>sfchflx_min, "max"=>sfchflx_max,
              "index"=>13)
  DCL.grcls

  rename_pngfile("SfcHFlxO_xtmean") if FlagOutputIMG
end

#------------------------------------------


u = get_gp("U")
ptemp = get_gp("PTemp")

msf = get_gp("MSF")
salt = get_gp("Salt")

totHT = get_gp("TotHT")
eulerHT = get_gp("EulerHT")
bolusHT = get_gp("BolusHT")
isoDiffHT = get_gp("IsoDiffHT")

bvFreq = get_gp("BVFreq")
densPot = get_gp("DensPot")

sfchflxo = get_gp("SfcHFlxO")

#------------------------------------------

p "u_ptemp fig"
u_ptemp_fig( u, ptemp, 1)
u_ptemp_fig( u, ptemp, 2)

p "msf_salt fig"
msf_salt_fig( msf, salt, 1)
msf_salt_fig( msf, salt, 2)

p "stratification fig"
stratification_fig( densPot, bvFreq, 1)
stratification_fig( densPot, bvFreq, 2)

p "meridional heat flux fig"
meridional_heat_flux_fig(totHT, eulerHT, isoDiffHT, bolusHT)

p "sfchflx fig"
sfchflx_fig(sfchflxo, 1)


