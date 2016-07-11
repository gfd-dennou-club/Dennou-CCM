require "numru/ggraph"
include NumRu

SRCTOPDIR = ARGV[0] #"/home/ykawai/workspace/DCPCM/run/com_dir/ocn/data_longInteg"
VARNAME = ARGV[1]
AVGPERIOD_CYCLE_BEGIN = ARGV[2].to_i
AVGPERIOD_CYCLE_END = ARGV[3].to_i
FNAME = ARGV[4]

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
  ncpathList.push(topDirPath+"/cycle#{n}-couple/"+NCName)
end
puts "VarName= #{VarName}"
puts "Cycle  = #{CycleBegin}..#{CycleEnd}"
puts "Axis =#{AxisNames}"
#------------------

ofile = NetCDF::create(distDirPath+"/#{VarName}.nc")

gphys = GPhys::IO.open(ncpathList[0], VarName)
ofile.copy_global_att(gphys)

gphys = GPhys::IO.open(ncpathList, VarName)
AxisNames.split(",").each{|axis|
  gphys = gphys.mean(axis)
}

# -- Output
GPhys::IO.write( ofile, gphys )
NetCDF_Conventions.add_history(ofile, File.basename($0)+" "+ncpathList[0])
ofile.close
