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
opt.on("-r", "--cyc_range <param>",  "the range of cycle"){|v| options[:cyc_range] = v}
opt.on("-v", "--var_name <param>",  "variable name of surface temperature"){|v| options[:var_name] = v}

opt.parse(ARGV)

CLIMATE_SNOWBALL = "SnowBallIce"
CLIMATE_PARTIALICE = "PartialIce"
CLIMATE_RUNAWAY = "Runaway"
CLIMATE_WARM = "Warm"

cyc_start = 1
cyc_end   = 1
if options[:cyc_range] != nil then
  cyc_start = options[:cyc_range].split(":")[0].to_i
  cyc_end   = options[:cyc_range].split(":")[1].to_i
end

sfctempvar_name = "SurfTemp"
if options[:var_name] != nil then
  sfctempvar_name = options[:var_name]
end

BeginCombineCyc=cyc_start
EndCombineCyc=cyc_end

SfcTempVarName = sfctempvar_name

TargetDir = (options[:exp_path] == nil) ? Dir.pwd : options[:exp_path]
DistDir = (options[:dist_path] == nil) ? Dir.pwd : options[:dist_path]
FlagOutputIMG = (options[:flag_output] == nil) ? false : options[:flag_output]
PrefixOutputIMG = (options[:imgfile_prefix] == nil) ? "" : options[:imgfile_prefix]
ClimateState = (options[:climate_state] == nil) ? CLIMATE_PARTIALICE : options[:climate_state]

PlanetName = "Earth"

DCL::swlset('lwnd', false) if FlagOutputIMG

def prep_dcl(iws, itr=1, clrmap=10)
#  DCLExt.sg_set_params('ifont'=>2)
  DCL.sgscmn(clrmap)           
  DCL.gropn(iws)              
  DCL.sgpset('isub', 96)     # 下付き添字を表す制御文字を '_' から '`' に
  DCL.glpset('lmiss',true)   # DCLの欠損値処理を on に．

  GGraph.set_fig('itr'=>itr) 
  GGraph.set_fig('viewport'=>[0.13,0.87,0.22,0.79])  
end

def rename_pngfile(fbasename)
  File.rename(DistDir+"/dcl_0001.png", DistDir+"/#{PrefixOutputIMG}#{fbasename}.png")
end

def merge_ncfile(ncname, varname, ovarname, cutOpt={}, meanOpt=[], ofilename=ovarname+".nc")
  ncpathList = []
  p "read NetCDF files (#{ncname}).."

  for i in BeginCombineCyc..EndCombineCyc
    p "#{i}." if i%10 == 0
    p "#{TargetDir}\/cycle#{i}-couple\/#{ncname}"
#    gphys = GPhys::IO.open(/#{TargetDir}\/cycle#{i}-couple\/#{ncname}/, varname)
    gphys = GPhys::IO.open("#{TargetDir}/cycle#{i}-couple/#{ncname}", varname)
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

def combine_ncfile_xmean_open(ncvarname, varname=ncvarname, ovarname=ncvarname, cutOpt={})
  p "combine & xmean (var=#{varname}).."
  merge_ncfile("#{ncvarname}.nc", varname, ovarname, cutOpt, ["lon"])

  p "open #{varname}.nc@#{varname} .."
  return GPhys::IO.open(ovarname+".nc", ovarname)
end

def sfctemp_xmean_tserise_fig(sfctemp,itr=1)
  prep_dcl(1,itr,10)

  case ClimateState
  when CLIMATE_SNOWBALL then
    sfctemp_min = 180; sfctemp_max = 320; sfctemp_int = 5
  when CLIMATE_RUNAWAY then
    sfctemp_min = 250; sfctemp_max = 450; sfctemp_int = 5
  else
    sfctemp_min = 180; sfctemp_max = 320; sfctemp_int = 5
  end

  GGraph::tone_and_contour( sfctemp, true, "titl"=>"SfcTemp",
                            "int"=>sfctemp_int, "max"=>sfctemp_max, "min"=>sfctemp_min, "exchange"=>true )
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.25)
  DCL.grcls

  rename_pngfile("SfcTempXMean_tserise") if FlagOutputIMG
end

def tempLowerLyr_xmean_tserise_fig(temp,itr=1)
  prep_dcl(1,itr,10)

  case ClimateState
  when CLIMATE_SNOWBALL then
    temp_min = 271; temp_max = 275; temp_int = 0.1
  when CLIMATE_RUNAWAY then
    temp_min = 271; temp_max = 300; temp_int = 1
  else
    temp_min = 271; temp_max = 280; temp_int = 0.2
  end
  
  GGraph::tone_and_contour( temp, true, "titl"=>"PTemp (depth 5 km)",
                            "int"=>temp_int, "max"=>temp_max, "min"=>temp_min, "exchange"=>true )
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.25)
  DCL.grcls

  rename_pngfile("PTempSig1.0XMean_tserise") if FlagOutputIMG
end

def tempMiddleLyr_xmean_tserise_fig(temp,itr=1)
  prep_dcl(1,itr,10)

  case ClimateState
  when CLIMATE_SNOWBALL then
    temp_min = 271; temp_max = 280; temp_int = 0.25
  when CLIMATE_RUNAWAY then
    temp_min = 271; temp_max = 320; temp_int = 1
  else
    temp_min = 271; temp_max = 290; temp_int = 0.5
  end
  
  GGraph::tone_and_contour( temp, true, "titl"=>"PTemp (depth 2.5km)",
                            "int"=>temp_int, "max"=>temp_max, "min"=>temp_min, "exchange"=>true )
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.25)
  DCL.grcls

  rename_pngfile("PTempSig0.5XMean_tserise") if FlagOutputIMG
end

def tempUpperLyr_xmean_tserise_fig(temp,itr=1)
  prep_dcl(1,itr,10)

  case ClimateState
  when CLIMATE_SNOWBALL then
    temp_min = 271; temp_max = 290; temp_int = 0.5
  when CLIMATE_RUNAWAY then
    temp_min = 271; temp_max = 350; temp_int = 2
  else
    temp_min = 271; temp_max = 310; temp_int = 1
  end
  
  GGraph::tone_and_contour( temp, true, "titl"=>"PTemp (depth 0.5km)",
                            "int"=>temp_int, "max"=>temp_max, "min"=>temp_min, "exchange"=>true )
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.25)
  DCL.grcls

  rename_pngfile("PTempSig0.1XMean_tserise") if FlagOutputIMG
end

tempLowerLyr = combine_ncfile_xmean_open("PTemp", "PTemp", "PTempSig1.0XMean", {"sig"=>-1.0})
tempLowerLyr_xmean_tserise_fig(tempLowerLyr)

tempMiddleLyr = combine_ncfile_xmean_open("PTemp", "PTemp", "PTempSig0.5XMean", {"sig"=>-0.5})
tempMiddleLyr_xmean_tserise_fig(tempMiddleLyr)

tempUpperLyr = combine_ncfile_xmean_open("PTemp", "PTemp", "PTempSig0.1XMean", {"sig"=>-0.1})
tempUpperLyr_xmean_tserise_fig(tempUpperLyr)

["PTempSig0.9", "PTempSig0.5", "PTempSig0.1"].each{|ovarname|
  ncpath = DistDir+"/#{ovarname}XMean.nc"
#  FileUtils::rm_f(ncpath) if FileTest::exist?(ncpath)
}
