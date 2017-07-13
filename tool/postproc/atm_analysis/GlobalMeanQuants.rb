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
require File.expand_path(File.dirname(__FILE__) + "/DCPAMUtil.rb")

eval("include ConstUtil::#{PlanetName}")
include DCModelIOUtil
include NumRu

#-----------------------------------------------------------------

@dcpamUtil = nil
VarDef = DCPAMUtil::VarNameDef
AxisDef = DCPAMUtil::AxisNameDef

#-----------------------------------------------------------------

@dcpamUtil = DCPAMUtil.new(PlanetName,"#{CurrentDir}/#{VarDef::U}.nc")

=begin
o2d_varList = [ "SurfTemp", "SurfAlbedo" ]
varList = [ 
	"PRCP", 
	"OLRA", "OSRA",
	"SSRA",
        "SLRA",
        "SensA",
        "EvapU"
]
p varList
=end
#=begin
o2d_varList = [ "o2d_SfcTemp", "o2d_SfcAlbedo" ]
varList = [ 
	"PRCP", 
	"OLRA", "OSRA",
	"SSRA", "SLRA", "SensA", "EvapA"
]
#=end
gp_PRCP, \
gp_OLR, gp_OSR, \
gp_SSR, gp_SLR, gp_Sens, gp_Evap  = GPhysUtil.get_GPhysObjs(varList, CurrentDir, nil, "_rank*")

gp_SfcTemp, gp_SfcAlbedo = GPhysUtil.get_GPhysObjs(o2d_varList, CurrentDir, nil, "_rank*")
time_pos = gp_OLR.axis("time").pos
gp_SfcTemp = gp_SfcTemp.interpolate(time_pos)
gp_SfcAlbedo = gp_SfcAlbedo.interpolate(time_pos)

ofile = NetCDF::create(OutputNCName)
GPhys::IO.each_along_dims_write(  
 [gp_SfcTemp, gp_SfcAlbedo, gp_PRCP, 
  gp_OLR, gp_OSR, 
  gp_SSR, gp_SLR, gp_Sens, gp_Evap], ofile, AxisDef::Time){
  |sfcTemp, sfcAlbedo, prcp, olr, osr, ssr, slr, sens, evap|

  time = sfcTemp.axis("time").pos
  puts "time=#{time.val[0]} [#{time.units}] .."

  radTOA = olr + osr
  netHFlxBOA = slr + ssr + evap + sens

  [ \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(sfcTemp),           \
                              "o2d_SfcTemp",                                    \
                              "global mean of surface temperature", "K" ),  \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(sfcAlbedo),         \
                              "o2d_SfcAlbedo",                                  \
                              "global mean of surface albedo", "1" ),       \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(prcp),             \
                              "PRCP",                                      \
                              "global mean of PRCP", "kg.m.s-1" ),         \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(olr),              \
                              "OLR",                                       \
                              "global mean of OLR", "W.m-2" ),             \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(-osr),             \
                             "mOSR",                                       \
                             "global mean of minus OSR", "W.m-2"),         \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(olr + osr),        \
                              "RadTOA",                                    \
                              "global mean of net radiation at TOA", "W.m-2" ), \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(slr),              \
                              "SLR",                                       \
                              "global mean of SLR", "W.m-2" ),             \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(ssr),              \
                              "SSR",                                       \
                              "global mean of SSR", "W.m-2" ),             \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(sens),              \
                              "Sens",                                       \
                              "global mean of Sens", "W.m-2" ),             \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(evap),              \
                              "Evap",                                       \
                              "global mean of Evap", "W.m-2" ),             \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(netHFlxBOA),        \
                              "HFlxBOA",                                    \
                              "global mean of SLR+SSR+Sens+Evap", "W.m-2" ),   \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(radTOA - netHFlxBOA),  \
                              "HBudgetAtm",                                    \
                              "heat budget", "W.m-2" ),     
  ]
}
ofile.close
