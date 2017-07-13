#!/usr/env ruby
require "numru/ggraph"
require "optparse"
require "parallel"
require "benchmark"
include NumRu

############################

varListOcn=["PTemp", "Salt"]
varListSfcFlx=["SfcHFlxO", "FreshWtFlxS"]
varListIce=["IceThick", "SnowThick", "SIceEn", "SIceCon"]

opt = OptionParser.new
options = {}
opt.on("-c", "--cyc_range <param>",  "the range of cycle"){|v| options[:cyc_range] = v}
opt.on("-r", "--run_rb_cyc_range <param>",  "the range of cycle over which GlobalMeanQuants.rb is run."){|v| options[:run_rb_cyc_range] = v}
opt.on("-p", "--periodic_coupling",  "the periodic coupling moode is used."){|v| options[:is_periodic_couple] = v}
opt.on("-i", "--interval_cyc <param>", "the interval of each cycle. (year)"){|v| options[:interval_cyc] = v}
opt.on("-P", "--ParallelMode",  "flag for parallel or serial mode"){|v| options[:is_parallel_mode] = v}

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
  intrv_cyc_standalone = options[:interval_cyc].split(",")[1].to_f #* 365.0
end

nproc = 1
if options[:is_parallel_mode] then
  nproc = Parallel.processor_count
end

p "cyc_range: #{cyc_start}--#{cyc_end}"
p "run_rb_cyc_range: #{run_rb_cyc_start}--#{run_rb_cyc_end}"
p "NProc=#{nproc}"

#---------------------------------------------

BeginGlMeanCyc=run_rb_cyc_start
EndGlMeanCyc=run_rb_cyc_end
BeginCombineCyc=cyc_start
EndCombineCyc=cyc_end

CMD_GLOBALMEAN="/home/ykawai/workspace/Dennou-CCM/tool/postproc/ocn_analysis/GlobalMeanQuants.rb"
CMD_MSF="/home/ykawai/workspace/Dennou-CCM/tool/postproc/ocn_analysis/GlobalMeanQuants.rb"

NProc = nproc
#---------------------------------------------

CoupledCycIntDay = UNumeric[intrv_cyc_couple, "day"]
StandaloneCycIntDay = UNumeric[intrv_cyc_standalone, "day"]

###########################

require "fileutils"
require "open3"

def exec_cmd(cmd, cyc=-1)
  lines_o = "";  lines_e = ""
  Open3.popen3(cmd) do |i, o, e, w|
    o.each do |line| lines_o <<  "#{cyc}:"+line  end
    e.each do |line| lines_e <<  "#{cyc}:"+line  end
  end
  puts "#{lines_o}";   puts "#{lines_e}"  
end

def combine_ncfile(beginCyc, endCyc, ofname, ifname, varname, modes)

  if File.exist?(ofname) then
    FileUtils.rm_f(ofname)
  end

  nCycle = endCyc - beginCyc + 1
  
  data_ary   = Array.new(nCycle)
  taxis_ary  = Array.new(nCycle)

  (beginCyc..endCyc).each{|n|
      i = n - beginCyc
      modes.each{|mode|
        fname = "cycle#{n}-#{mode}/#{ifname}"
        if !File.exist?(fname) then
          p "#{fname} is not found. Check!"
        end
        gp = GPhys::IO.open(fname, varname)
        time_pos = gp.axis("time").pos
        tlen = time_pos.length

        gp = gp.cut('time'=>time_pos.val[1]..time_pos.val[tlen-1])

        time_day = (time_pos.val[1..tlen-1] - time_pos.val[0]) + (CoupledCycIntDay + StandaloneCycIntDay)*(n-1)
        time_day = time_day + CoupledCycIntDay if mode == "standalone"
        time_year = time_day /UNumeric[365.0, 'day/year']        
        #    gp.axis('time').set_pos(VArray.new(time_year, nil, 'year'))

        data_ary[i] =  gp.val[true].to_a
        taxis_ary[i] =  time_year.to_a
      }
  }
  
  #  p "fname=#{fnames.join(" ")}"
  #  p data_ary

  puts "cobine ncfile: #{varname} period: #{beginCyc}..#{endCyc}"

  gp0 = GPhys::IO.open("cycle#{beginCyc}-couple/#{ifname}", varname)
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

bench_result = Benchmark.realtime do
  ret = Parallel.map(BeginGlMeanCyc..EndGlMeanCyc, :in_processes =>NProc){|i|
    #(BeginGlMeanCyc..EndGlMeanCyc).each{|i|
    modes.each{|mode|
      Dir.chdir("cycle#{i}-#{mode}"){
        p "run GlobalMeanQuants.rb( dir=#{Dir::pwd} )"
        exec_cmd(CMD_GLOBALMEAN)
      
        #p "run MSF.rb( dir=#{Dir::pwd} )"
        #exec_cmd(CMD_MSF)
      }
    }
  }
end
puts "elapse time (cmd_globalmean) : #{bench_result} s"

bench_result = Benchmark.realtime do
  varList = []
  var2ncfile = Hash.new
  varListOcn.each{|var| varList.push(var); var2ncfile[var] = "GlobalMeanQuants.nc"}
  varListSfcFlx.each{|var| varList.push(var); var2ncfile[var] = "GlobalMeanQuants_SfcFlx.nc"}
  varListIce.each{|var| varList.push(var); var2ncfile[var] = "GlobalMeanQuants_SIce.nc"}
  
  ret = Parallel.map(varList, :in_processes =>NProc){|var|
    combine_ncfile(BeginCombineCyc, EndCombineCyc, "#{var}GlMean-standalone.nc", var2ncfile[var], var, modes)    
  }  
end
puts "elapse time (combine_ncfile) : #{bench_result} s"
