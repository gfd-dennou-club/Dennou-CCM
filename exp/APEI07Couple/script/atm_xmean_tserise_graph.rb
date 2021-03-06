#!/usr/env ruby
# coding: utf-8
require "numru/ggraph"
require "optparse"
require "fileutils"
require "parallel"

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

NProc=8 #Parallel.processor_count

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
  p "read NetCDF files (#{ncname}).."

  nCyc = EndCombineCyc - BeginCombineCyc + 1
  nSubCyc = 10
  nBlock = nCyc / nSubCyc
  nBlock += 1 if (nCyc % nSubCyc != 0)

  subcyc_fname = Hash.new
  index_info = Hash.new
  for i in 0..nBlock-1
    beginSubCyc = i*nSubCyc + 1
    endSubCyc = [beginSubCyc + nSubCyc - 1, EndCombineCyc].min
    index_info[i] = { "begin_cyc"=>beginSubCyc, "end_cyc"=>endSubCyc }
    subcyc_fname[i] = "#{DistDir}/tmp#{beginSubCyc}-#{endSubCyc}_#{ofilename}"
#    p index_info[i]
  end
  
  ret = Parallel.map(0..nBlock-1, :in_processes =>NProc ){|b|

    id_start = index_info[b]["begin_cyc"]
    id_end = index_info[b]["end_cyc"]

    gphys_list = []
    tmp_file_list = []
    
    p "block=#{b} (cyc=#{id_start}:#{id_end}) #{TargetDir}\/cycle*-couple\/#{ncname}"
    (id_start..id_end).each{|i|
      gphys = GPhys::IO.open(/#{TargetDir}\/cycle#{i}-couple\/#{ncname}/, varname)

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
      tlen = gphys.axis("time").length
      gphys_list.push(gphys[false,1..tlen-1])
    }

    ofile = NetCDF::create(subcyc_fname[b])    
    GPhys::IO.write(ofile, GPhys.join(gphys_list), ovarname)
    ofile.close
  }

  p "Output.."
  ofile = NetCDF::create("#{DistDir}/#{ofilename}")
  gphys = GPhys::IO.open(subcyc_fname.values, ovarname)
  gphys = gphys.mean("time") if meanOpt.include?("time")
  
  GPhys::IO.write(ofile, gphys)
#  NetCDF_Conventions.add_history(ofile, File.basename($0)+" "+ncpathList[0])  
  ofile.close

  p "remove tmp files.."
  FileUtils::rm_f(subcyc_fname.values)
end

def combine_ncfile_xmean_open(ncvarname, varname=ncvarname, ovarname=ncvarname, cutOpt={})
  p "combine & xmean (var=#{varname}).."
  merge_ncfile("#{ncvarname}_rank(\\d\\d\\d\\d\\d\\d).nc", varname, ovarname, cutOpt, ["lon"])

  p "open #{varname}.nc@#{varname} .."
  gp = GPhys::IO.open(ovarname+".nc", ovarname)

  if (EndCombineCyc - BeginCombineCyc + 1 >= 1000) then
    time = gp.axis("time").pos
    tlen = time.length
    gp = GPhys::IO.open_gturl("#{ovarname}.nc@#{ovarname},time=#{time.val[0]}:#{time.val[tlen-1]}:2")
  end
  return gp
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
    temp_min = 180; temp_max = 320; temp_int = 5
  when CLIMATE_RUNAWAY then
    temp_min = 250; temp_max = 450; temp_int = 5
  else
    temp_min = 180; temp_max = 320; temp_int = 5
  end
  
  GGraph::tone_and_contour( temp, true, "titl"=>"Temp (sigma=0.9)",
                            "int"=>temp_int, "max"=>temp_max, "min"=>temp_min, "exchange"=>true )
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.25)
  DCL.grcls

  rename_pngfile("TempSig0.9XMean_tserise") if FlagOutputIMG
end

def tempUpperLyr_xmean_tserise_fig(temp,itr=1)
  prep_dcl(1,itr,10)

  case ClimateState
  when CLIMATE_SNOWBALL then
    temp_min = 80; temp_max = 220; temp_int = 5
  when CLIMATE_RUNAWAY then
    temp_min = 250; temp_max = 450; temp_int = 5
  else
    temp_min = 150; temp_max = 290; temp_int = 5
  end
  
  GGraph::tone_and_contour( temp, true, "titl"=>"Temp (sigma=0.3)",
                            "int"=>temp_int, "max"=>temp_max, "min"=>temp_min, "exchange"=>true )
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.25)
  DCL.grcls

  rename_pngfile("TempSig0.3XMean_tserise") if FlagOutputIMG
end

sfctemp = combine_ncfile_xmean_open(SfcTempVarName, SfcTempVarName, "SurfTempXMean")
sfctemp_xmean_tserise_fig(sfctemp)

tempLowerLyr = combine_ncfile_xmean_open("^Temp", "Temp", "TempSig0.9XMean", {"sig"=>0.9})
tempLowerLyr_xmean_tserise_fig(tempLowerLyr)

tempUpperLyr = combine_ncfile_xmean_open("^Temp", "Temp", "TempSig0.3XMean", {"sig"=>0.3})
tempUpperLyr_xmean_tserise_fig(tempUpperLyr)

["SurfTemp", "TempSig0.9", "TempSig0.3"].each{|ovarname|
  ncpath = DistDir+"/#{ovarname}XMean.nc"
#  FileUtils::rm_f(ncpath) if FileTest::exist?(ncpath)
}
