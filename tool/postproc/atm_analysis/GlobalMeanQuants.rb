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

puts "CurrentDir=#{CurrentDir} .."
@dcpamUtil = DCPAMUtil.new(PlanetName, "#{CurrentDir}/#{VarDef::U}.nc")

varList = [ 
	"SurfTempOcn", "PRCP", 
	"OLRA", "OSRA",
	"SSRA", "SLRA", "SensA", "EvapA"
]

gp_SurfTemp, gp_PRCP, \
gp_OLR, gp_OSR, \
gp_SSR, gp_SLR, gp_Sens, gp_Evap \
 = GPhysUtil.get_GPhysObjs(varList)

ofile = NetCDF::create(OutputNCName)
GPhys::IO.each_along_dims_write(  
 [gp_SurfTemp, gp_PRCP, 
  gp_OLR, gp_OSR, 
  gp_SSR, gp_SLR, gp_Sens, gp_Evap], ofile, AxisDef::Time){
  |surfTemp, prcp, olr, osr, ssr, slr, sens, evap|

  time = surfTemp.axis("time").pos
  puts "time=#{time.val[0]} [#{time.units}] .."

  radTOA = olr + osr
  netHFlxBOA = slr + ssr + evap + sens

  [ \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(surfTemp),         \
                              "SurfTemp",                                  \
                              "global mean of surface temperature", "K" ), \
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
