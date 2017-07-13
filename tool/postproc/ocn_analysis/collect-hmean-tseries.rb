#!/usr/env ruby
require "numru/ggraph"
require "optparse"
include NumRu

############################

varListOcn=["PTemp"] #, "Salt"]
varListSfcFlx=[] #["SfcHFlxO", "FreshWtFlxS"]
varListIce=[] #"IceThick", "SnowThick", "SIceEn", "SIceCon"]

opt = OptionParser.new
options = {}
opt.on("-c", "--cyc_range <param>",  "the range of cycle"){|v| options[:cyc_range] = v}
opt.on("-r", "--run_rb_cyc_range <param>",  "the range of cycle over which GlobalMeanQuants.rb is run."){|v| options[:run_rb_cyc_range] = v}
opt.on("-p", "--periodic_coupling",  "the periodic coupling moode is used."){|v| options[:is_periodic_couple] = v}
opt.on("-i", "--interval_cyc <param>", "the interval of each cycle. (year)"){|v| options[:interval_cyc] = v}
opt.parse(ARGV)

cyc_start = 1
cyc_end   = 1
if options[:cyc_range] != nil then
  cyc_start = options[:cyc_range].split(":")[0].to_i
  cyc_end   = options[:cyc_range].split(":")[1].to_i
end

run_rb_cyc_start = 1
run_rb_cyc_end   = 1
if options[:run_rb_cyc_range] != nil then
  run_rb_cyc_start = options[:run_rb_cyc_range].split(":")[0].to_i
  run_rb_cyc_end   = options[:run_rb_cyc_range].split(":")[1].to_i
end

intrv_cyc_couple     = 2.0 * 365.0
intrv_cyc_standalone = 50.0 * 365.0
if options[:is_periodic_couple] then
  intrv_cyc_standalone = options[:interval_cyc].split(",")[1].to_f * 365.0
end

p "cyc_range: #{cyc_start}--#{cyc_end}"
p "run_rb_cyc_range: #{run_rb_cyc_start}--#{run_rb_cyc_end}"

#---------------------------------------------

BeginGlMeanCyc=run_rb_cyc_start
EndGlMeanCyc=run_rb_cyc_end
BeginCombineCyc=cyc_start
EndCombineCyc=cyc_end

CMD_HORIZONMEAN="/home/ykawai/workspace/Dennou-CCM/tool/postproc/ocn_analysis/HorizontalMeanQuants.rb"
CMD_MSF="/home/ykawai/workspace/Dennou-CCM/tool/postproc/ocn_analysis/GlobalMeanQuants.rb"

#---------------------------------------------

CoupledCycIntDay = UNumeric[intrv_cyc_couple, "day"]
StandaloneCycIntDay = UNumeric[intrv_cyc_standalone, "day"]

###########################

require "fileutils"
require "open3"

def exec_cmd(cmd)
  Open3.popen3(cmd) do |i, o, e, w|
    o.each do |line| p line end
  end
end

def horizontal_mean()
  gp = GPhys::IO.open("PTemp.nc", "PTemp")
  ofile = NetCDF.create("HorizontalMeanQuants.nc")

  lon_w = GPhys::IO.open("PTemp.nc", "lon_weight")
  lat_w = GPhys::IO.open("PTemp.nc", "lat_weight")
  [0.0, 0.2, 0.4, 0.6, 0.8, 1.0].each{|sig|
    gp_sig = gp.cut("sig"=>-1.0*sig)
    # horizontal mean
    gp_sig_hmean = (gp_sig*lon_w*lat_w).sum("lon").sum("lat")/(4.0*NMath::PI)
    gp_sig_hmean.units = gp.units
    GPhys::IO.write(ofile, gp_sig_hmean.rename(gp.name+"_sig#{sig.to_s}"))

    [0, 20, 40, 60, 80].each{|lat|
      gp_sig_local = (gp_sig.cut("lat"=>lat)*lon_w ).sum("lon")/(2.0*NMath::PI)
      gp_sig_local.units = gp.units
      GPhys::IO.write(ofile, gp_sig_local.rename(gp.name+"_sig#{sig.to_s}_lat#{lat.to_s}"))      
    }
  }
  ofile.close
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
      fname = "cycle#{i}-#{mode}/#{ifname}"
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

#------------------------------------------------------------

modes = ["couple"]
modes.push("standalone") if options[:is_periodic_couple]
p modes
p BeginCombineCyc
p EndGlMeanCyc

for i in BeginGlMeanCyc..EndGlMeanCyc
  modes.each{|mode|
    Dir.chdir("cycle#{i}-#{mode}"){
      p "run GlobalMeanQuants.rb( dir=#{Dir::pwd} )"
#      exec_cmd(CMD_GLOBALMEAN)
      horizontal_mean()      
    }
  }
end

for var in varListOcn
  [0.0, 0.2, 0.6, 1.0].each{|sig|
    var_ = var + "_sig#{sig.to_s}"
    combine_ncfile(BeginCombineCyc, EndCombineCyc, "#{var_}HorMean-standalone.nc", "HorizontalMeanQuants.nc", var_, modes)
    [0, 20, 40, 60, 80].each{|lat|
      var_lat = var + "_sig#{sig.to_s}_lat#{lat.to_s}"
      combine_ncfile(BeginCombineCyc, EndCombineCyc, "#{var_lat}-standalone.nc", "HorizontalMeanQuants.nc", var_lat, modes)
    }
  }
end

for var in varListSfcFlx
  combine_ncfile(BeginCombineCyc, EndCombineCyc, "#{var}GlMean-standalone.nc", "GlobalMeanQuants_SfcFlx.nc", var, modes)
end

for var in varListIce
  combine_ncfile(BeginCombineCyc, EndCombineCyc, "#{var}GlMean-standalone.nc", "GlobalMeanQuants_SIce.nc", var, modes)
end
