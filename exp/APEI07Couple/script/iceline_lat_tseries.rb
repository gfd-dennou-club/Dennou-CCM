#!/usr/env ruby
require "numru/ggraph"
require "optparse"
require "parallel"
include NumRu

opt = OptionParser.new
options = {}
opt.on("-c", "--cyc_range <param>",  "the range of cycle"){|v| options[:cyc_range] = v}
opt.on("-r", "--run_rb_cyc_range <param>",  "the range of cycle over which the searching of iceline latitude is performed."){|v| options[:run_rb_cyc_range] = v}
opt.on("-v", "--var_name <param>",  "variable name of surface temperature"){|v| options[:var_name] = v}
opt.on("-p", "--periodic_coupling",  "the periodic coupling moode is used."){|v| options[:is_periodic_couple] = v}
opt.on("-i", "--interval_cyc <param>", "the interval of each cycle. (day)"){|v| options[:interval_cyc] = v}
opt.parse(ARGV)

cyc_start = 1
cyc_end   = 1
if options[:cyc_range] != nil then
  cyc_start = options[:cyc_range].split(":")[0].to_i
  cyc_end   = options[:cyc_range].split(":")[1].to_i
end

run_rb_cyc_start = cyc_start
run_rb_cyc_end   = cyc_end
if options[:run_rb_cyc_range] != nil then
  run_rb_cyc_start = options[:run_rb_cyc_range].split(":")[0].to_i
  run_rb_cyc_end   = options[:run_rb_cyc_range].split(":")[1].to_i
end

intrv_cyc_couple     = 2.0 * 365.0
intrv_cyc_standalone = 50.0 * 365.0
if options[:interval_cyc] != nil then
  intrv_cyc_couple = options[:interval_cyc].split(",")[0].to_f
  intrv_cyc_standalone = options[:interval_cyc].split(",")[1].to_f
end

var_name = "SurfTemp"
if options[:var_name] != nil then
  var_name = options[:var_name]
end

TARGET_DIR=Dir.pwd

BeginRBCyc=run_rb_cyc_start
EndRBCyc=run_rb_cyc_end
BeginCombineCyc=cyc_start
EndCombineCyc=cyc_end

ICELINE_TEMP=UNumeric[263.0, "K"]
SFCTEMP_VARNAME=var_name

CoupledCycIntDay = UNumeric[intrv_cyc_couple, "day"]
StandaloneCycIntDay = UNumeric[intrv_cyc_standalone, "day"]

p "cyc_range: #{cyc_start}--#{cyc_end}"
p "the period of coupled/stadalone run: #{CoupledCycIntDay}:#{StandaloneCycIntDay} [day]"

#---------------------------------------------

NProc=Parallel.processor_count

p "NProc=#{NProc}"

###########################

require "fileutils"
require "open3"

def exec_cmd(cmd)
  Open3.popen3(cmd) do |i, o, e, w|
    o.each do |line| p line end
  end
end

def combine_ncfile(beginCyc, endCyc, ofname, ifname, varname, modes)

  if File.exist?(ofname) then
    FileUtils.rm_f(ofname)
  end

  
  fnames = []
  data_ary   = []
  taxis_ary  = []
  tot_tlen = 0
  gp0 = nil
  
  (beginCyc..endCyc).each{|i|
    modes.each{|mode|
      fname = TARGET_DIR+"/cycle#{i}-couple/iceline_lat.nc" # "cycle#{i}-#{mode}/#{ifname}"
      if !File.exist?(fname) then
        p "#{fname} is not found. Check!"
      end
      fnames.push(fname)
      gp = GPhys::IO.open(fname, varname)
      time_pos = gp.axis("time").pos
      tlen = time_pos.length

      tot_tlen = tot_tlen + tlen - 1
      gp = gp.cut('time'=>time_pos.val[1]..time_pos.val[tlen-1])
      gp0 = gp if i == beginCyc

      time_day = (time_pos.val[1..tlen-1] - time_pos.val[0]) + (CoupledCycIntDay + StandaloneCycIntDay)*(i-1)
      time_day = time_day + CoupledCycIntDay if mode == "standalone"
      time_year = time_day /UNumeric[365.0, 'day/year']        
      #    gp.axis('time').set_pos(VArray.new(time_year, nil, 'year'))
      data_ary.push( gp.val[true].to_a )
      taxis_ary.push( time_year.to_a )
    }
  }
  p fnames.join(" ")

  va = VArray.new( NArray.to_na(data_ary.flatten),
                   {'long_name'=>gp0.long_name, 'units'=>gp0.units.to_s},
                   gp0.name )
  taxis = Axis.new.set_pos( VArray.new(NArray.to_na(taxis_ary.flatten), {'long_name'=>'time', 'units'=>'year'}, 'time') )

  gp_combine = GPhys.new(Grid.new(taxis), va)
  ofile = NetCDF.create(ofname)
  GPhys::IO.write(ofile, gp_combine)
  ofile.close
end

def extract_iceline(cycle, varname)
  ncname = TARGET_DIR + "/cycle#{cycle}-couple/#{varname}_rank(\\d\\d\\d\\d\\d\\d).nc"
  p "ncname=#{ncname}"
  gp_SfcTemp = GPhys::IO.open(/#{ncname}/, varname)
  gp_Tg_meanlon = gp_SfcTemp.mean("lon")
  ncname = TARGET_DIR + "/cycle#{cycle}-couple/U.nc"
  gp_lat_weight = GPhys::IO.open(ncname, "lat_weight")
  
  time_axis = gp_Tg_meanlon.axis("time")
  tlen = time_axis.pos.length
  lat_axis = gp_Tg_meanlon.axis("lat")
  latlen = lat_axis.pos.length

  gp_Tg_dev     = gp_Tg_meanlon - ICELINE_TEMP
  gp_sign_check = gp_Tg_dev[1..latlen-1,true]*gp_Tg_dev[0..latlen-2,true]
    
  ofile = NetCDF.create(TARGET_DIR+"/cycle#{cycle}-couple/iceline_lat.nc")

  gp_iceline = GPhys.new(
    Grid.new(time_axis),
    VArray.new(NArray.sfloat(tlen), {"long_name"=> "ice line latitude", "units"=>"degree_north"}, "iceline_lat")
  )

  gp_iceline_nh = GPhys.new(
    Grid.new(time_axis),
    VArray.new(NArray.sfloat(tlen), {"long_name"=> "ice line latitude (north hemisphere)", "units"=>"degree_north"}, "iceline_lat_nh")
  )
  gp_iceline_sh = GPhys.new(
    Grid.new(time_axis),
    VArray.new(NArray.sfloat(tlen), {"long_name"=> "ice line latitude (south hemisphere)", "units"=>"degree_north"}, "iceline_lat_sh")
  )
  
  gp_Tg_glmean = GPhys.new(
    Grid.new(time_axis),
    VArray.new(NArray.sfloat(tlen), {"long_name"=> "global mean temperature", "units"=>"degree_north"}, "SfcTemp")
  )

  for n in 0..tlen-1

    icelat_nh = 90.0
    icelat_sh = 90.0

    # p gp_Tg_dev.val[true,n]
    #    p gp_sign_check.val[0..latlen-2,n]

    for j in 0..latlen-2
      if (gp_sign_check[j,n] <= 0.0) then
        tmp_icelat = ( (gp_Tg_dev.val[j+1,n].abs*lat_axis.pos[j] + gp_Tg_dev.val[j,n].abs*lat_axis.pos[j+1]) / (gp_Tg_dev.val[j,n] - gp_Tg_dev.val[j+1,n]).abs ).abs

        if lat_axis.pos.val[j] > 0.0 then
          if (gp_Tg_dev[j,n] < gp_Tg_dev[j+1,n]) and (icelat_sh == 90.0 and icelat_nh == 90.0) then
            icelat_sh = 0.0
#            p "accross eq.."
          else
            icelat_nh = tmp_icelat
          end
        else
          if (gp_Tg_dev[j,n] > gp_Tg_dev[j+1,n]) and (icelat_sh != 90.0) then
            icelat_nh = 0.0
#            p "accross eq.."
          else
            icelat_sh = tmp_icelat
          end
        end
      end
    end
    
    if (icelat_nh == 90.0 and icelat_sh == 90.0) and gp_Tg_dev[latlen/2,n] < 0.0 then
      gp_iceline[n] = 0.0
      gp_iceline_nh[n] = 0.0
      gp_iceline_sh[n] = 0.0
    else
      gp_iceline_nh[n] = icelat_nh
      gp_iceline_sh[n] = icelat_sh
      gp_iceline[n] = 0.5*(icelat_nh + icelat_sh)      
    end

    gp_Tg_glmean[n] = ((gp_Tg_meanlon[true,n]*gp_lat_weight).sum * 0.5).to_f
    p "cyc=#{cycle}: n=#{n} iceline=#{gp_iceline.val[n]},  sfctemp=#{gp_Tg_glmean.val[n]}"
    p "nh=", icelat_nh, "sh=", icelat_sh
  end
  
  GPhys::IO.write(ofile, gp_iceline)
  GPhys::IO.write(ofile, gp_iceline_nh)
  GPhys::IO.write(ofile, gp_iceline_sh)
  GPhys::IO.write(ofile, gp_Tg_glmean)
  ofile.close
end

cycles = (BeginRBCyc..EndRBCyc).map{|i| i}
ret = Parallel.map(cycles, :in_processes => NProc){|i|
  extract_iceline(i, SFCTEMP_VARNAME)
}

[""].each{|suffix|
#["", "_nh", "_sh"].each{|suffix|
  ncname = "iceline_lat#{suffix}.nc"
  if File.exist?(ncname) then
    FileUtils::rm_f(ncname)
  end
  combine_ncfile(BeginCombineCyc, EndCombineCyc, "iceline_lat#{suffix}.nc", "tmp*-iceline_lat.nc", "iceline_lat#{suffix}", ["couple"])
}
#exec_cmd("gpcat -o iceline_lat.nc -v iceline_lat tmp*-iceline_lat.nc")

if File.exist?("sfctemp_glmean.nc") then
  FileUtils::rm_f("sfctemp_glmean.nc")  
end
#exec_cmd("gpcat -o sfctemp_glmean.nc -v SfcTemp tmp*-iceline_lat.nc")
combine_ncfile(BeginCombineCyc, EndCombineCyc, "sfctemp_glmean.nc", "tmp*-iceline_lat.nc", "SfcTemp", ["couple"])

#FileUtils::rm_f(Dir.glob("tmp*-iceline_lat.nc"))
