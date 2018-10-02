#!/usr/env ruby
# coding: utf-8

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
require "optparse"

opt = OptionParser.new
options = {}
opt.on("-c", "--const_util  <param>",  "name of module for physical constants"){|v| options[:const_util_name] = v}
opt.parse(ARGV)

CONST_UTIL_NAME = (options[:const_util_name] == nil) ? "ConstUtil" : options[:const_util_name]
p "const_util: #{CONST_UTIL_NAME}"

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
p Atm::CpDry
p LatentHeatV


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
  rho = press/(Atm::GasRDry*temp)

  dryStatEn = rho*(Atm::CpDry*temp + Grav*height)
  latentEn = rho*(LatentHeatV*qvap)


  [ \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalIntLonSig(dryStatEn*v/rho / 1e15, ps),                   \
                 "dryStatEnFlxLat",                                                                     \
                 "amount of meridional dry static energy flux integrated in vertical column", "PW") ,   \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalIntLonSig(latentEn*v/rho / 1e15, ps),                    \
                 "latentEnFlxLat",                                                                      \
                 "amount of meridional latent energy flux integrated in vertical column", "PW") ,       \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalIntLonSig((dryStatEn + latentEn)*v/rho / 1e15, ps),      \
                 "moistStatEnFlxLat",                                                                   \
                 "amount of meridional moist static energy flux integrated in vertical column", "PW") , \
  ]
}
ofile.close



