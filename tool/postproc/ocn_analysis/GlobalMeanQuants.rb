#!/usr/bin/env ruby

#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------

# Configuration **********************************

OutputNCName = 'GlobalMeanQuants.nc'
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
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMean3D(ptemp, h),                      \
                              "PTemp",                                                 \
                              "global mean of surface temperature", "K" ),             \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMean3D(salt, h),                       \
                              "Salt",                                                  \
                              "global mean of salinity", "psu" ),                      \
  ]
}
ofile.close

#-----------------------------------------------------------------------------------------

varSfcFlxList = [ "SfcHFlxO", "FreshWtFlxS" ]
gp_SfcHFlxO, gp_FreshWtFlxS = GPhysUtil.get_GPhysObjs( varSfcFlxList )

ofile = NetCDF::create(OutputNCName_SfcFlx)
GPhys::IO.each_along_dims_write(
  [gp_SfcHFlxO, gp_FreshWtFlxS, gp_H], ofile, AxisDef::Time){
  | sfcHFlxO, freshWtFlxS, h|

  time = sfcHFlxO.axis("time")
  puts "SfcFlx: time=#{time.pos.val[0]} [#{time.pos.units}] .."
  
  [ \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSfc(sfcHFlxO),                    \
                              "SfcHFlxO",                                             \
                              "global mean of surface heat flux", "W.m-2" ),          \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSfc(freshWtFlxS),                 \
                              "FreshWtFlxS",                                          \
                              "global mean of freshwater flux", "m.s-1" ),	      \
#    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSurf(surfHFlxAI),                 \
#                              "SurfHFlxAI",                                           \
#                              "global mean of surface heatr flux (AI)", "W.s-3" ),    \
#    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSurf(surfHFlxAIO),                 \
#                              "SurfHFlxAIO",                                           \
#                              "global mean of surface heatr flux (AIO)", "W.s-3" )
  ]  
}
ofile.close

#-----------------------------------------------------------------------------------------


ofile = NetCDF::create(OutputNCName_SIce)
varList = [ "IceThick", "SnowThick", "SIceCon", "SIceEn", "Wice" ]
gp_IceThick, gp_SnowThick, gp_SIceCon, gp_SIceEn, gp_Wice = GPhysUtil.get_GPhysObjs(varList, Dir::pwd, "history_sice.nc")

#gp_SIceEnSum = gp_SIceEn.sum("sig2")

GPhys::IO.each_along_dims_write( \
  [gp_IceThick, gp_SnowThick, gp_SIceCon, gp_SIceEn, gp_Wice], ofile, AxisDef::Time){ \
  |iceThick, snowThick, siceCon, siceEn, wice|

  time = iceThick.axis("time")
  puts "time=#{time.pos.val[0]} [#{time.pos.units}] .."

  [ \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSfc(iceThick),                    \
                              "IceThick",                                             \
                              "global mean of ice effective thickness", "m" ), 	      \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSfc(snowThick),                   \
                              "SnowThick",                                            \
                              "global mean of snow effective thickness", "m" ),       \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSfc(siceCon),                     \
                              "SIceCon"  ,                                            \
                              "global mean of sea ice thickness", "1"),               \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSfc(siceEn.sum('sig2')),          \
                              "SIceEn"  ,                                            \
                              "global mean of sea ice energy", "J/m2"),               \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSfc(wice),                       \
                              "Wice"  ,                                                   \
                              "global mean of sea ice creation and melting", "kg.m-2.s-1"),   \
  ]
}
ofile.close

#-----------------------------------------------------------------------------------------

=begin
gp_SurfHFlxO, gp_SurfFwFlxO, gp_SurfHFlxAI, gp_SurfHFlxAO, gp_SIceCon \
 = GPhysUtil.get_GPhysObjs(varSurfFlxList.push("SIceCon"))

ofile = NetCDF::create("SurfHFlxAIO.nc")
GPhys::IO.each_along_dims_write(
	[gp_SurfHFlxAI, gp_SurfHFlxAO, gp_SIceCon], ofile, AxisDef::Time){ \
  |surfHFlxAI, surfHFlxAO, siceCon|

  time = surfHFlxAI.axis("time")
  puts "time=#{time.pos.val[0]} [#{time.pos.units}] .."

  surfHFlxAIO = siceCon*surfHFlxAI + (1.0 - siceCon)*surfHFlxAO
  [ \
    GPhysUtil.redef_GPhysObj( surfHFlxAIO,                  \
                              "SurfHFlxAIO",                                            \
                              "global mean of surface heat flux (AIO)", "W.m-2" ) ]
}
				
ofile.close

#-----------------------------------------------------------------------------------------

ofile = NetCDF::create(OutputNCName_SurfFlx)
GPhys::IO.each_along_dims_write( \
  [gp_SurfHFlxO, gp_SurfFwFlxO, gp_SurfHFlxAI, gp_SurfHFlxAO, gp_SIceCon], ofile, AxisDef::Time){ \
  |surfHFlxO, surfFwFlxO, surfHFlxAI, surfHFlxAO, siceCon|

  time = surfHFlxO.axis("time")
  puts "time=#{time.pos.val[0]} [#{time.pos.units}] .."

  surfHFlxAIO = siceCon*surfHFlxAI + (1.0 - siceCon)*surfHFlxAO

  [ \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSurf(surfHFlxO),                  \
                              "SurfHFlxO",                                            \
                              "global mean of surface heat flux", "W.m-2" ),          \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSurf(surfFwFlxO),                 \
                              "SurfFwFlxO",                                           \
                              "global mean of freshwater flux", "m.s-1" ),	      \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSurf(surfHFlxAI),                 \
                              "SurfHFlxAI",                                           \
                              "global mean of surface heatr flux (AI)", "W.s-3" ),    \
    GPhysUtil.redef_GPhysObj( @dsogcmUtil.globalMeanSurf(surfHFlxAIO),                 \
                              "SurfHFlxAIO",                                           \
                              "global mean of surface heatr flux (AIO)", "W.s-3" )
  ]
}
ofile.close


#-----------------------------------------------------------------------------------------

gp_PTemp_glmean = GPhys::NetCDF_IO.open(OutputNCName, "PTemp")
gp_SurfHFlxO_glmean = GPhys::NetCDF_IO.open(OutputNCName_SurfFlx, "SurfHFlxO")

ax_time = gp_SurfHFlxO_glmean.axis("time")
tlen = ax_time.length
gp_time = ax_time.to_gphys

totTime           = UNumeric[(gp_time.val[tlen-1] - gp_time.val[0])*86400.0, "sec"]
tot_surfHeating = - @dsogcmUtil.globalIntLonLat(gp_SurfHFlxO.sum("time"))*(totTime/tlen.to_f)
ocnMass = @dsogcmUtil.globalIntLonLat(Ocn::Dens0*gp_TotDepthBasic)
dtemp_surfHeating = tot_surfHeating/(Ocn::Cp0*ocnMass)

p "DTemp (due to surface heat flux) = #{dtemp_surfHeating.val[0]} [K]"
p "DTemp/Dt (due to surface heat flux) = #{dtemp_surfHeating.val[0]/totTime*86400.0} [K/day]"
p "----------------------------------"

ax_time = gp_PTemp_glmean.axis("time")
tlen = ax_time.length
dtemp_glmean = gp_PTemp_glmean.val[tlen-1] - gp_PTemp_glmean.val[0]
p "Actual change of global mean temperature = #{dtemp_glmean} [K]"
=end
