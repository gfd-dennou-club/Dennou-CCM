#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------

# Configuration **********************************

OutputNCName = 'Diagnose.nc'
CurrentDir=Dir::pwd
PlanetName = 'Earth'

#**************************************************

require "numru/ggraph"
require "optparse"

opt = OptionParser.new
options = {}
opt.on("-c", "--const_util  <param>",  "name of module for physical constants"){|v| options[:const_util_name] = v}
opt.parse(ARGV)

CONST_UTIL_NAME = (options[:const_util_name] == nil) ? "ConstUtil" : options[:const_util_name]
p "const_util: #{CONST_UTIL_NAME}"

#-----------------------------------------------------------------

require File.expand_path(File.dirname(__FILE__) + "/../common/#{CONST_UTIL_NAME}.rb")
require File.expand_path(File.dirname(__FILE__) + "/../common/DCModelIOUtil.rb")
require File.expand_path(File.dirname(__FILE__) + "/DCPAMUtil.rb")

eval("include #{CONST_UTIL_NAME}::#{PlanetName}")
include DCModelIOUtil
include NumRu

#-----------------------------------------------------------------

@dcpamUtil = nil
VarDef = DCPAMUtil::VarNameDef
AxisDef = DCPAMUtil::AxisNameDef

#-----------------------------------------------------------------

puts "CurrentDir=#{CurrentDir} .."
@dcpamUtil = DCPAMUtil.new(PlanetName, "#{CurrentDir}/#{VarDef::U}.nc")

#varList = [VarDef::U, VarDef::V, VarDef::Temp, VarDef::Ps, VarDef::QH2OVap, VarDef::Height]
#gp_U, gp_V, gp_Temp, gp_Ps, gp_QVap, gp_Height = GPhysUtil.get_GPhysObjs(varList, CurrentDir, nil, "_rank*")
varList = [VarDef::V, VarDef::Ps, VarDef::Temp, VarDef::QH2OVap]
gp_V, gp_Ps, gp_Temp, gp_QH2OVap = GPhysUtil.get_GPhysObjs(varList, CurrentDir, nil, "_rank*")


ofile = NetCDF::create(OutputNCName)
GPhys::IO.each_along_dims_write( 
#  [gp_U, gp_V, gp_Temp, gp_Ps, gp_QVap, gp_Height], ofile, AxisDef::Time){
#  |u, v, temp, ps, qvap, height|
  [gp_V, gp_Ps, gp_Temp, gp_QH2OVap], ofile, AxisDef::Time){
  |v, ps, temp, qvap|

  time = v.axis("time").pos
  puts "time=#{time.val[0]} [#{time.units}] .."

  msf = @dcpamUtil.calc_MSF(v, ps)

#=begin  
  press = @dcpamUtil.calc_Pressure(ps)
  e_sat = UNumeric[611.0, "Pa"]
  qvap_sat = Atm::EpsV * e_sat * (LatentHeatV/Atm::GasRWat*(1.0/UNumeric[273.0,"K"] - 1.0/temp)).exp / press
  relhumd = qvap/qvap_sat * 100.0
#=end
  
  [ \
    GPhysUtil.redef_GPhysObj( msf/1e9,                  \
                 "MSF", "mass stream function(note: We define 1Sv=10^9 [kg/s].", "Sv"), \
    GPhysUtil.redef_GPhysObj( relhumd,                  \
                 "RelHumd", "relative humidity.", "%")
  ]
}

ofile.close
