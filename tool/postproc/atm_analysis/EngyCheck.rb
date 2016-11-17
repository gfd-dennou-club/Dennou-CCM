require "numru/gphys"
include NumRu

NRank = 4
NCVarNames = ['IntEngy', 'KinEngy', 'LatEngy', 'PotEngy', 'TotEngy', 'Mass']

#######################################

def output_ncfile(fname, *gphys_ary)
  ofile = NetCDF.create(fname)
  gphys_ary.each{|gphys|
    GPhys::IO.write(ofile, gphys)
  }
  ofile.close
end

###############

NCVarNames.each{|varname|
  fname = Dir.glob("#{varname}_rank*.nc").sort
  p "varname=#{varname} nRank=0..#{fname.length-1}"
  gphys = GPhys::IO.open(fname[0], varname).copy
  for rank in 1..fname.length-1
    gphys = gphys + GPhys::IO.open(fname[rank], varname)
  end
  output_ncfile("#{varname}.nc", gphys)
}

################

