#-------------------------------------------------------------
# Copyright (c) 2015-2015 Kawai Yuta. All rights reserved.
#-------------------------------------------------------------

require "numru/ggraph"
include NumRu

###############################################################

module DCModelIOUtil

  module GPhysUtil

    def self.get_GPhysObjs(varNames, currentDir=Dir::pwd, ncname=nil, ncname_suffix="")

      gp_Array = []
      time_len_min = 1e15
      varNames.each{|varName|
        ncfname = []
        if (ncname == nil) then
          ncfname = Dir.glob("#{currentDir}/#{varName}#{ncname_suffix}.nc")
        else
          ncfname = Dir.glob("#{currentDir}/#{ncname.gsub(".nc","")}#{ncname_suffix}.nc")
        end

        gp = GPhys::IO.open(ncfname,varName)
        if gp.axnames.include?("time") then
          time = gp.axis("time")
          time_len_min = time.length if time.length < time_len_min
        end
      }
      varNames.each{|varName|
        ncfname = []
        if (ncname == nil) then
          ncfname = Dir.glob("#{currentDir}/#{varName}#{ncname_suffix}.nc")
        else
          ncfname = Dir.glob("#{currentDir}/#{ncname.gsub(".nc","")}#{ncname_suffix}.nc")
        end

        gp = GPhys::IO.open(ncfname,varName)
        if gp.axnames.include?("time") then
          time = gp.axis("time")
          gp = gp.cut("time"=>time.pos.val[0]..time.pos.val[time_len_min-1])
        end
        gp_Array.push(gp)
      }

      return gp_Array
    end

    def self.redef_GPhysObj(gphysOri, newname, long_name, units=nil)
      gp = gphysOri.rename(newname)
      if units != nil
        gp.units = units
      else
        gp.units = gphysOri.units
      end
      gp.long_name = long_name
      return gp
    end
  end
  
end
