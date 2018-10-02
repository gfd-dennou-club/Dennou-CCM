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
CLIMATE_PARTICE_COLD = "PartIceCold"
CLIMATE_PARTICE_LARGE = "PartIceLARGE"


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
E3 = UNumeric[5.2e3, "m"]
TotDepth = UNumeric[5.2e3, "m"]
CoupledCycIntDay    = 730.0
StandaloneCycIntDay = 18250.0
TAxisThinIntrv      = 2

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
  GGraph.set_axes('ytitle'=>'')
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
    beginSubCyc = i*nSubCyc + BeginCombineCyc
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
    (id_start..id_end).each{|n|
      i = n - id_start
      ["couple", "standalone"].each{|mode|
        fname = "#{TargetDir}/cycle#{n}-#{mode}/#{ncname}"
        gphys = GPhys::IO.open(fname, varname)

        gphys = gphys.cut(cutOpt) if cutOpt.length > 0

        meanOpt.each{|opt|
          case opt      
          when "time"
          else
            gphys = gphys.mean(opt)
          end
        }

        gphys = gphys.copy
        time_pos = gphys.axis("time").pos
        tlen = time_pos.length

        time_day = (time_pos.val[1..tlen-1] - time_pos.val[0]) + (CoupledCycIntDay + StandaloneCycIntDay)*(n-1)
        time_day = time_day + CoupledCycIntDay if mode == "standalone"
        time_year = time_day /UNumeric[365.0, 'day/year']        
        taxis = Axis.new.set_pos( VArray.new(time_year, {'long_name'=>'time', 'units'=>'year'}, 'time') )

        grid = gphys.grid
        time_dim = grid.dim_index("time")
        grid.delete_axes(time_dim)
        grid.set_axis(time_dim, taxis)

        gphys_list.push( GPhys.new(grid, gphys.data[false,1..tlen-1]) )
      }
    }

    ofile = NetCDF::create(subcyc_fname[b])
    GPhys::IO.write(ofile, GPhys.join(gphys_list), ovarname)
    ofile.close
  }

  p "Output.."
  ofile = NetCDF::create("#{DistDir}/#{ofilename}")
  gphys = GPhys::IO.open(subcyc_fname.values, ovarname)
  GPhys::IO.write(ofile, gphys)
#  NetCDF_Conventions.add_history(ofile, File.basename($0)+" "+ncpathList[0])  
  ofile.close

  p "remove tmp files.."
  FileUtils::rm_f(subcyc_fname.values)
end

def combine_ncfile_xmean_open(ncvarname, varname=ncvarname, ovarname=ncvarname, cutOpt={})
  p "combine & xmean (var=#{varname}).."
  merge_ncfile("#{ncvarname}.nc", varname, ovarname, cutOpt, ["lon"])

  p "open #{varname}.nc@#{varname} .."
  gp = GPhys::IO.open(ovarname+".nc", ovarname)
  if (EndCombineCyc - BeginCombineCyc + 1 >= 200) then
    time = gp.axis("time").pos
    tlen = time.length
    p "#{ovarname}.nc@#{ovarname},time=#{time.val[0]}:#{time.val[tlen-1]}:#{TAxisThinIntrv}"
    gp = GPhys::IO.open_gturl("#{ovarname}.nc@#{ovarname},time=#{time.val[0]}:#{time.val[tlen-1]}:#{TAxisThinIntrv}")
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

def z_to_figz(z, upperLyrDepth)
    return (z > -upperLyrDepth) ? z :
                - upperLyrDepth + upperLyrDepth*(upperLyrDepth + z)/(TotDepth - upperLyrDepth)
end

def sig_to_figZgrid(grid, itr, sig_name='sig')
  gp_sig = grid.axis(sig_name).to_gphys
  nz = gp_sig.shape[0]
  na_z = NArray.sfloat(nz)

  upperLyrDepth = UNumeric[1e3,"m"]
#  upperLyrDepth = UNumeric[0.5e3,"m"]

  for k in 0..nz-1
    z = E3*gp_sig.val[k]
    na_z[k] = (itr == 2)? z_to_figz(z, upperLyrDepth) : z
  end

  ry1 = (0..26).map{|i|
    z = -i*200.0
    (itr == 2)? z_to_figz(z, upperLyrDepth) : z
  }
  
#  cy = ['0.0', '0.2', '0.4', '0.5', '1.0', '2.0', '3.0', '4.0', '5.0']
  cy = ['0.0', '0.2', '0.4', '0.6', '0.8', '1.0', '2.0', '3.0', '4.0', '5.0']
  ry2 = cy.map{|z|
    z_ = -z.to_f*1e3
    (itr == 2)? z_to_figz(z_, upperLyrDepth) : z_
  }
  va_z = VArray.new(na_z, {"long_name"=>"ocean depth", "units"=>"m"}, "depth")
  return grid.change_axis(0, Axis.new.set_pos(va_z)), ry1, ry2, cy
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

def tempZ_xmean_tserise_fig(temp_,itr=1,lat)
  prep_dcl(1,1,40)

  case ClimateState
  when CLIMATE_SNOWBALL then
    temp_min = 271; temp_max = 290; temp_int = 0.5
  when CLIMATE_RUNAWAY then
    temp_min = 271; temp_max = 350; temp_int = 2
  else
    temp_min = 271; temp_max = 300; temp_int = 1
  end

  fig_z_grid, ry1, ry2, cy = sig_to_figZgrid(temp_.grid,itr)
  temp = GPhys.new(fig_z_grid, temp_.data).convert_units("degC")
  
  GGraph.set_axes('yside'=>'u')

  rmiss = DCL.glpget('rmiss')  
  temp_levels = [rmiss,-1.0].concat( gen_levels(1.0, 29.0, 2.0))#.concat([rmiss])
  temp_levels_s = temp_levels.map{|level| level.to_i.to_s }  
  GGraph::tone_and_contour( temp, true, "titl"=>"PTemp (#{lat} degrees_north)",
                            "levels"=>temp_levels, "label"=>temp_levels_s, "exchange"=>true, "index"=>[3,2] )

  DCL::uyaxlb('L', ry1, ry2, cy, 4)
  DCL::uyaxlb('R', ry1, ry2, cy, 4)
  DCL::uysttl('L', 'depth (km)', 0.0)  
  GGraph.color_bar("charfact"=>0.75, "vlength"=>0.5, "units"=>"deg C")
  DCL.grcls

  rename_pngfile("PTempLat#{lat}XMean_tserise") if FlagOutputIMG
end

#tempLowerLyr = combine_ncfile_xmean_open("PTemp", "PTemp", "PTempSig1.0XMean", {"sig"=>-1.0})
#tempLowerLyr_xmean_tserise_fig(tempLowerLyr)

#tempMiddleLyr = combine_ncfile_xmean_open("PTemp", "PTemp", "PTempSig0.5XMean", {"sig"=>-0.5})
#tempMiddleLyr_xmean_tserise_fig(tempMiddleLyr)

#tempUpperLyr = combine_ncfile_xmean_open("PTemp", "PTemp", "PTempSig0.1XMean", {"sig"=>-0.1})
#tempUpperLyr_xmean_tserise_fig(tempUpperLyr)

temp_lat70 = combine_ncfile_xmean_open("PTemp", "PTemp", "PTempLat70XMean", {"lat"=>70})
tempZ_xmean_tserise_fig(temp_lat70, 2, 0)

temp_lat30 = combine_ncfile_xmean_open("PTemp", "PTemp", "PTempLat30XMean", {"lat"=>30})
tempZ_xmean_tserise_fig(temp_lat30, 2, 0)

temp_lat0 = combine_ncfile_xmean_open("PTemp", "PTemp", "PTempLat0XMean", {"lat"=>0})
tempZ_xmean_tserise_fig(temp_lat0, 2, 0)

["PTempSig1.0", "PTempSig0.5", "PTempSig0.1",
 "PTempLat0XMean"].each{|ovarname|
  ncpath = DistDir+"/#{ovarname}XMean.nc"
  FileUtils::rm_f(ncpath) if FileTest::exist?(ncpath)
}
