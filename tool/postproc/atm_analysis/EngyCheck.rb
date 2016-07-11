require "numru/gphys"
include NumRu

NRank = 4
NCVarNames = ['IntEngy', 'KinEngy', 'LatEngy', 'PotEngy', 'TotEngy', 'Mass']
IntWeightFile = "OLR.nc"

#######################################

LonIntWeight = GPhys::IO.open(IntWeightFile, "lon_weight")
LatIntWeight = GPhys::IO.open(IntWeightFile, "lat_weight")

def globalMean(gphys)
  tmp = LonIntWeight*gphys
  return (LatIntWeight*tmp).sum("lon","lat")/(4.0*Math::PI)
end

def output_ncfile(fname, *gphys_ary)
  ofile = NetCDF.create(fname)
  gphys_ary.each{|gphys|
    GPhys::IO.write(ofile, gphys)
  }
  ofile.close
end

#######################################
=begin
gphys_OLR = GPhys::IO.open("OLR.nc", "OLR")
gphys_OSR = GPhys::IO.open("OSR.nc", "OSR")
  
intOLR = globalMean(gphys_OLR)
intOLR.rename("OLR")

intOSR = globalMean(gphys_OSR)
intOSR.rename("OSR")

radNetTOA = gphys_OLR + gphys_OSR
radNetTOA.rename("ONetRTOA")
intRadNetTOA = intOLR + intOSR
intRadNetTOA.rename("ONetR")

output_ncfile("ONetRTOA.nc",radNetTOA)
output_ncfile("RadBudgetTOA.nc", intOLR, intOSR, intRadNetTOA)
=end

###############

NCVarNames.each{|varname|
  fname_myrank = []
  for rank in 0..NRank-1
    fname_myrank[rank] = "#{varname}_rank%06d.nc"%rank
  end
#  puts fname_myrank

  gphys = GPhys::IO.open(fname_myrank[0], varname)
  for rank in 1..NRank-1
    gphys = gphys + GPhys::IO.open(fname_myrank[rank], varname)
  end

  output_ncfile("#{varname}.nc", gphys)
}

################

