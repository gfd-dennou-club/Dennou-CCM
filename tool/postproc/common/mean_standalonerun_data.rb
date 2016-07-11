####################

#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------
#
# * Usage
# If you want to average data of cycle A to cycle B (where A and B are the number of cycle)
# along some axises,  execute the following command. 
#
# ./mean_couplerun_data.rb A B ncname varname axisnames (TopDirPath DistDirPath)
#
# Configuration **********************************

require "numru/gphys"
require "numru/gphys/gpcommon"
include NumRu

CurrentDir=Dir::pwd
PlanetName = 'Earth'

CycleBegin = ARGV[0]
CycleEnd = ARGV[1]
NCName = ARGV[2]
VarName = ARGV[3]
AxisNames = ARGV[4]
topDirPath = ARGV[5]
distDirPath = ARGV[6]

#**************************************************

#-----------------------------------------------------------------

#-----------------------------------------------------------------

topDirPath = CurrentDir if topDirPath == nil
distDirPath = CurrentDir if distDirPath == nil

puts "TopSRCDirPath=#{topDirPath} .."


ncpathList = []
for n in CycleBegin..CycleEnd
  ncpathList.push(topDirPath+"/cycle#{n}-standalone/"+NCName)
end
puts "VarName= #{VarName}"
puts "Cycle  = #{CycleBegin}..#{CycleEnd}"
puts "Axis =#{AxisNames}"
#------------------

ofile = NetCDF::create(distDirPath+"/#{VarName}.nc")

gphys = GPhys::IO.open(ncpathList[0], VarName)
ofile.copy_global_att(gphys)

AxisNames.split(",").each{|axis|
  gphys = gphys.mean(axis)
}

for n in 1..ncpathList.length-1
  gphys_ = GPhys::IO.open(ncpathList[n], VarName)
  AxisNames.split(",").each{|axis|
    gphys_ = gphys_.mean(axis)
  }
  gphys = (gphys*n.to_f + gphys_)/(n + 1).to_f
end

# -- Output
GPhys::IO.write( ofile, gphys )
NetCDF_Conventions.add_history(ofile, File.basename($0)+" "+ncpathList[0])
ofile.close
