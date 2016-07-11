#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------

# Configuration **********************************

OutputNCName = 'HeatBudget.nc'
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

varNameList = [ "SurfTempOcn", "OLRA", "OSRA", "RadSDWFLXA", "RadSUWFLXA", "RadLDWFLXA", "RadLUWFLXA" ]

gp_SurfTemp, gp_OLR, gp_OSR, gp_RadSDWFlx, gp_RadSUWFlx, gp_RadLDWFlx, gp_RadLUWFlx,  \
 = GPhysUtil.get_GPhysObjs(varNameList)

ofile = NetCDF::create("RadNetTOA.nc")
GPhys::IO.each_along_dims_write( 
  [gp_OLR, gp_OSR],
  ofile, AxisDef::Time){
  |olr, osr|
  [ \
    GPhysUtil.redef_GPhysObj( olr+osr,                                      \
                              "RadNetTOA",                                  \
                              "global mean of OLR", "W.m-2" ),              \
  ]
}
ofile.close

ofile = NetCDF::create(OutputNCName)
GPhys::IO.each_along_dims_write( 
  [gp_SurfTemp, gp_OLR, gp_OSR, gp_RadSDWFlx, gp_RadSUWFlx, gp_RadLDWFlx, gp_RadLUWFlx],
  ofile, AxisDef::Time){
  |surfTemp, olr, osr, radSDWFlx, radSUWFlx, radLDWFlx, radLUWRFlx|

  time = surfTemp.axis("time").pos
  puts "time=#{time.val[0]} [#{time.units}] .."

  [ \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(radSDWFlx.cut('sigm'=>0.0)),  \
                              "SDWRFlxTOA",                                           \
                              "global mean of downward short wave at TOA", "W.m-2" ), \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(radSUWFlx.cut('sigm'=>0.0)),  \
                              "SUWRFlxTOA",                                           \
                              "global mean of upward short wave at TOA", "W.m-2" ),   \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(radSDWFlx.cut('sigm'=>1.0)),  \
                              "SDWRFlxBOA",                                           \
                              "global mean of upward short wave at BOA", "W.m-2" ),   \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(radSUWFlx.cut('sigm'=>1.0)),  \
                              "SUWRFlxBOA",                                           \
                              "global mean of upward short wave at BOA", "W.m-2" ),   \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(olr),              \
                              "OLR",                                       \
                              "global mean of OLR", "W.m-2" ),             \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(-osr),             \
                             "mOSR",                                       \
                             "global mean of minus OSR", "W.m-2"),         \
    GPhysUtil.redef_GPhysObj( @dcpamUtil.globalMeanSurf(olr + osr),        \
                              "RadTOA",                                    \
                              "global mean of net radiation at TOA", "W.m-2" ), \

    GPhysUtil.redef_GPhysObj( (radSUWFlx -radSDWFlx).mean('lon'),  \
                              "SRFlx",                                           \
                              "net short wave radiation", "W.m-2" ),   \
    GPhysUtil.redef_GPhysObj( (radLUWRFlx -radLDWFlx).mean('lon'),  \
                              "LRFlx",                                           \
                              "net long wave radiation", "W.m-2" ),   \
  ]
}
ofile.close
