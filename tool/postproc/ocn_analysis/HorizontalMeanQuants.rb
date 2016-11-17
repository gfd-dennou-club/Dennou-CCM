#!/usr/bin/env ruby

#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------

# Configuration **********************************

OutputNCName = 'HorizontalMeanQuants.nc'
CurrentDir=Dir::pwd
PlanetName = 'Earth'

#**************************************************

require "numru/ggraph"
require File.expand_path(File.dirname(__FILE__) + "/../common/ConstUtil.rb")
require File.expand_path(File.dirname(__FILE__) + "/../common/DCModelIOUtil.rb")
require File.expand_path(File.dirname(__FILE__) + "/DennouOGCMUtil.rb")

eval("include ConstUtil::#{PlanetName}")
include DCModelIOUtil
include NumRu

#-----------------------------------------------------------------

@dsogcmUtil = nil
OutputNCName_SfcFlx = OutputNCName.gsub(".nc", "_SfcFlx.nc")
OutputNCName_SIce = OutputNCName.gsub(".nc", "_SIce.nc")
VarDef = DennouOGCMUtil::VarNameDef
AxisDef = DennouOGCMUtil::AxisNameDef

#------------------------------------------------------------------------------------------

puts "CurrentDir=#{CurrentDir} .."
@dsogcmUtil = DennouOGCMUtil.new(PlanetName, "#{CurrentDir}/#{VarDef::U}.nc")

varList = [ VarDef::PTemp, VarDef::Salt, VarDef::H ]
gp_PTemp, gp_Salt, gp_H  = GPhysUtil.get_GPhysObjs( varList )

#-----------------------------------------------------------------------------------------

ofile = NetCDF::create(OutputNCName)
GPhys::IO.each_along_dims_write( \
  [gp_PTemp, gp_Salt, gp_H], ofile, AxisDef::Time){ \
  |ptemp, salt, h|

  time = ptemp.axis("time")
  puts "CommonVar: time=#{time.pos.val[0]} [#{time.pos.units}] .."

  [ \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanLonLat(ptemp),                      \
                              "PTemp",                                                  \
                              "global horizontal mean of surface temperature", "K" ),   \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanLonLat(salt),                       \
                              "Salt",                                                   \
                              "global horizontal mean of salinity", "psu" ),            \
  ]
}
ofile.close
