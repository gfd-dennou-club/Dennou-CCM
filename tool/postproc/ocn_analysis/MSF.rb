#!/usr/bin/env ruby

#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------

# Configuration **********************************

OutputNCName = 'MSF.nc'
CurrentDir=Dir::pwd
PlanetName = 'Earth'

#**************************************************

#-----------------------------------------------------------------

require "numru/ggraph"
require "optparse"
require File.expand_path(File.dirname(__FILE__) + "/../common/ConstUtil.rb")
require File.expand_path(File.dirname(__FILE__) + "/../common/DCModelIOUtil.rb")
require File.expand_path(File.dirname(__FILE__) + "/DennouOGCMUtil.rb")

eval("include ConstUtil::#{PlanetName}")
include DCModelIOUtil
include NumRu

#-----------------------------------------------------------------

opt = OptionParser.new
options = {}
opt.on("-t", "--time_range <param>",  "time range"){|v| options[:time_range] = v}
opt.parse(ARGV)

tstart = 0.0
tend   = -1.0
if options[:time_range] != nil then
  tstart = options[:time_range].split(":")[0].to_f
  tend   = options[:time_range].split(":")[1].to_f
  p "time range: #{tstart}--#{tend}"
end

#-----------------------------------------------------------------

VarDef = DennouOGCMUtil::VarNameDef
AxisDef = DennouOGCMUtil::AxisNameDef

#-----------------------------------------------------------------

puts "CurrentDir=#{CurrentDir} .."
@dogcmUtil = DennouOGCMUtil.new(PlanetName, "#{CurrentDir}/#{VarDef::U}.nc")

varList = [VarDef::V, 'H']
gp_V, gp_H = GPhysUtil.get_GPhysObjs(varList)
if tend > 0.0 then
  gp_V = gp_V.cut("time"=>tstart..tend)
  gp_H = gp_H.cut("time"=>tstart..tend)
end

ofile = NetCDF::create(OutputNCName)
GPhys::IO.each_along_dims_write( 
  [gp_V, gp_H], ofile, AxisDef::Time){
  |v, h|

  time = v.axis("time").pos
  puts "time=#{time.val[0]} [#{time.units}] .."

  msf = @dogcmUtil.calc_MSF(v, h)

  [ \
    GPhysUtil.redef_GPhysObj( msf/1e9,                  \
                 "MSF", "mass stream function(note: We define 1Sv=10^9 [kg/s].", "Sv")
  ]
}

ofile.close
