#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------

# Configuration **********************************

OutputNCName = 'EngyFlx.nc'
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
gp_U, gp_V, gp_Temp, gp_Ps, gp_QVap, gp_Height = GPhysUtil.get_GPhysObjs(varList, CurrentDir, nil, "_rank*")

p gp_Temp.shape
p gp_Ps.shape
ofile = NetCDF::create(OutputNCName)
GPhys::IO.each_along_dims_write( 
  [gp_U, gp_V, gp_Temp, gp_Ps, gp_QVap, gp_Height], ofile, AxisDef::Time){
  |u, v, temp, ps, qvap, height|

  time = u.axis("time").pos
  puts "time=#{time.val[0]} [#{time.units}] .."
  press = @dcpamUtil.calc_Pressure(ps)
  rho = @dcpamUtil.calc_Density(press, temp)

  dryStatEn = rho*(Atm::CpDry*temp + Grav*height)
  moistStatEn = rho*(LatentHeatV*qvap)


  [ \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalIntLonSig(dryStatEn*v/rho / 1e15, ps),                  \
                 "dryStatEnFlxLat",                                                                     \
                 "amount of meridional dry static energy flux integrated in vertical column", "PW") ,   \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalIntLonSig(moistStatEn*v/rho / 1e15, ps),                \
                 "moistStatEnFlxLat",                                                                   \
                 "amount of meridional moist static energy flux integrated in vertical column", "PW") , \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalIntLonSig((dryStatEn + moistStatEn)*v/rho / 1e15, ps),  \
                 "totStatEnFlxLat",                                                                     \
                 "amount of meridional total static energy flux integrated in vertical column", "PW") , \
  ]
}

ofile.close


