#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------

# Configuration **********************************

OutputNCName = 'Diagnose.nc'
CurrentDir=Dir::pwd
PlanetName = 'Earth'

#**************************************************

#-----------------------------------------------------------------

require "numru/ggraph"
require File.expand_path(File.dirname(__FILE__) + "/../common/ConstUtil.rb")
require File.expand_path(File.dirname(__FILE__) + "/../common/DCModelIOUtil.rb")
require File.expand_path(File.dirname(__FILE__) + "/DCPAMUtil.rb")

eval("include ConstUtil::#{PlanetName}")
include DCModelIOUtil
include NumRu

#-----------------------------------------------------------------

@dcpamUtil = nil
VarDef = DCPAMUtil::VarNameDef
AxisDef = DCPAMUtil::AxisNameDef

#-----------------------------------------------------------------

puts "CurrentDir=#{CurrentDir} .."
@dcpamUtil = DCPAMUtil.new(PlanetName, "#{CurrentDir}/#{VarDef::U}.nc")

varList = [VarDef::U, VarDef::V, VarDef::Temp, VarDef::Ps, VarDef::QH2OVap, VarDef::Height]
gp_U, gp_V, gp_Temp, gp_Ps, gp_QVap, gp_Height = GPhysUtil.get_GPhysObjs(varList)


ofile = NetCDF::create(OutputNCName)
GPhys::IO.each_along_dims_write( 
  [gp_U, gp_V, gp_Temp, gp_Ps, gp_QVap, gp_Height], ofile, AxisDef::Time){
  |u, v, temp, ps, qvap, height|

  time = u.axis("time").pos
  puts "time=#{time.val[0]} [#{time.units}] .."

  msf = @dcpamUtil.calc_MSF(v, ps)

  [ \
    GPhysUtil.redef_GPhysObj( msf/1e9,                  \
                 "MSF", "mass stream function(note: We define 1Sv=10^9 [kg/s].", "Sv")
  ]
}

ofile.close
