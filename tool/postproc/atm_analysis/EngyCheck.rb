require "numru/gphys"
include NumRu

NRank = 4
NCVarNames = ['IntEngy', 'KinEngy', 'LatEngy', 'PotEngy', 'TotEngy', 'Mass']
RPlanet = UNumeric[6.371e6, 'm']

#######################################

def output_ncfile(fname, *gphys_ary)
  ofile = NetCDF.create(fname)
  gphys_ary.each{|gphys|
    GPhys::IO.write(ofile, gphys)
  }
  ofile.close
end

def int_latband(fname, varname)
  gp = GPhys::IO.open(fname, varname)
  latband_weight = GPhys::IO.open(fname, "lat_weight").sum("lat")
  lonband_weight = GPhys::IO.open(fname, "lon_weight").sum("lon")
  area = latband_weight * lonband_weight * RPlanet**2
  return  gp * area
end

###############


NCVarNames.each{|varname|
  fname = Dir.glob("#{varname}_rank*.nc").sort
  p "varname=#{varname} nRank=0..#{fname.length-1}"
  gphys = int_latband(fname[0], varname)
  for rank in 1..fname.length-1
    gphys = gphys + int_latband(fname[rank], varname)
  end
  gphys.units = "J"
  output_ncfile("#{varname}.nc", gphys)
}

################

