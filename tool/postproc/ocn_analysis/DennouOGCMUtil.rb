require "numru/ggraph"
require File.expand_path(File.dirname(__FILE__) + "/../common/ConstUtil.rb")
include NumRu
include ConstUtil

class DennouOGCMUtil
  attr_accessor :ax_Lon, :ax_Lat, :ax_Sig
  attr_accessor :gp_lonIntWt, :gp_latIntWt, :gp_sigIntWt
  attr_accessor :planetName
  attr_accessor :lonAxisName, :latAxisName, :sigAxisName
  attr_accessor :const
  attr_accessor :globalVolume, :globalSurfArea

  class VarNameDef
    U  = 'U'
    V  = 'V'
    W  = 'W'
    OMG  = 'OMG'
    Salt  = 'Salt'
    PTemp  = 'PTemp'
    SSH = 'SSH'
    H  = 'H'
    SfcPres = 'SfcPres'    

    BolusU = 'BolusU'
    BolusV = 'BolusV'
  end

  class AxisNameDef
    Lon = 'lon'
    Lat = 'lat'
    Sig = 'sig'
    Time = 'time'
  end
  
  def initialize(planetName, ncpath)

    @planetName = planetName

    @gp_lon = GPhys::IO.open_gturl(ncpath+"@lon")
    @gp_lat = GPhys::IO.open_gturl(ncpath+"@lat")
    @gp_sig = GPhys::IO.open_gturl(ncpath+"@sig")    
    @gp_lonIntWt = GPhys::IO.open_gturl(ncpath+"@lon_weight")
    @gp_latIntWt = GPhys::IO.open_gturl(ncpath+"@lat_weight")
    @gp_sigIntWt = GPhys::IO.open_gturl(ncpath+"@sig_weight")


    @lonAxisName = "lon"
    @latAxisName = "lat"
    @sigAxisName = "sig"
    @sigmAxisName = "sigm"
    @timeAxisName = "time"

    @ax_Lon = @gp_lon.axis(@lonAxisName)
    @ax_Lat = @gp_lat.axis(@latAxisName)
    @ax_Sig = @gp_sig.axis(@sigAxisName)


    nz = @ax_Sig.pos.length
    na_sigm = NArray.float(nz+1)
    na_sigm[0] = 0.0
    for k in 1..nz
      na_sigm[k] = na_sigm[k-1] - @gp_sigIntWt.val[k-1]
    end
    v_sigm = VArray.new(na_sigm, {"long_name"=>"sigma coordinate(mid point)", "units"=>"1"}, "sigm")
    @ax_sigm = Axis.new.set_pos(v_sigm)
    
    @const = eval("ConstUtil::#{@planetName}")

    @globalSurfArea = 4.0*PI*@const::RPlanet**2
    
#    puts "Initialize an object of DCPAMUtil class.."
  end

  def globalIntLonSig(gphys, gp_H)
    tmp = ((@gp_sigIntWt*gphys*gp_H).sum(@sigAxisName))
    return @const::RPlanet * (PI/180.0*@gp_lat).cos * ( @gp_lonIntWt*tmp).sum(@lonAxisName)
  end
  
  def globalIntLonLat(gphys)
    return @const::RPlanet**2 * (@gp_lonIntWt*(@gp_latIntWt*gphys)).sum(@lonAxisName, @latAxisName)
  end
  
  def globalIntLonLatSig(gphys, gp_H)
    tmp = ((@gp_sigIntWt*gphys*gp_H).sum(@sigAxisName))
    return @const::RPlanet**2 * (@gp_lonIntWt*(@gp_latIntWt*tmp)).sum(@lonAxisName, @latAxisName)
  end

  def globalMeanSfc(gphys)
    return globalIntLonLat(gphys)/@globalSurfArea
  end

  def globalMeanLonLat(gphys)
    return globalIntLonLat(gphys)/@globalSurfArea
  end
  
  def globalMean3D(gphys, gp_H)
    gp_one = gp_H.copy
    gp_one[false] = 1.0
    totVol = globalIntLonLatSig(gp_one, gp_H)
    return globalIntLonLatSig(gphys, gp_H)/totVol
  end

  def calc_MSF(gp_V, gp_H)
    gp_MSF = gen_3DGPysObj("MSF", "mass stream function", "kg.s-1", gp_H.axis(@timeAxisName), @ax_sigm)

    cos_phi = (@gp_lat * PI/180.0).cos
    alph = gp_V * cos_phi * @const::RPlanet * PI * 2.0 * @const::Ocn::Dens0

    kmax = @ax_Sig.pos.length
    gp_MSF[false] = 0.0
    for kk in 1..kmax
      k = kmax - kk
      gp_MSF[false,k,true] = gp_MSF[false,k+1,true] + \
                             alph[false,k,true]*gp_H[false,k,true]*@gp_sigIntWt[k]
    end
    return gp_MSF
  end
  
  def gen_3DGPysObj(name, long_name, units, timeAxis=nil, ax_Z=@ax_Sig)
    na = nil; grid = nil
    if timeAxis != nil
#      p @ax_Lon, @ax_Lat, ax_Z, timeAxis
      na = NArray.sfloat(@ax_Lon.length, @ax_Lat.length, ax_Z.length, timeAxis.length)
      grid = Grid.new(@ax_Lon, @ax_Lat, ax_Z, timeAxis)
    else
      na = NArray.sfloat(@ax_Lon.length, @ax_Lat.length, ax_Z.length)
      grid = Grid.new(@ax_Lon, @ax_Lat, ax_Z)
    end

    va = VArray.new(na, {"name"=>name, "long_name"=>long_name, "units"=>units})
    return GPhys.new(grid, va)
  end
  
end
