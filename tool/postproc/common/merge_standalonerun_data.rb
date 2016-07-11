#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------
#
# * Usage
# If you want to merge data of cycle A to cycle B (where A and B are the number of cycle),
# execute the following command. 
#
# ./merge_couplerun_data.rb A B ncname varname (TopDirPath DistDirPath)
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
topDirPath = ARGV[4]
distDirPath = ARGV[5]

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

#------------------

#ofile = NetCDF::create(distDirPath+"/#{VarName}_standalone.nc")

gphys = GPhys::IO.open(ncpathList[0], VarName)
#ofile.copy_global_att(gphys)

gphys = GPhys::IO.open(ncpathList, VarName)

p gphys

# -- Output
#GPhys::IO.write( ofile, gphys )
#NetCDF_Conventions.add_history(ofile, File.basename($0)+" "+ncpathList[0])
#ofile.close


