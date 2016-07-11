require "numru/ggraph"
include NumRu

module ConstUtil

  GasRUniv = UNumeric[8.3144621, "J.K-1.mol-1"]
  PI = UNumeric[Math::acos(-1.0), "1"]
  StB = UNumeric[5.67e-8, "W.m-2.K-4"]

  module Earth
    RPlanet  = UNumeric[6371e3, "m"]
    Grav     = UNumeric[9.81, "m.s-2"]
    Omega    = UNumeric[2.0*Math::acos(-1.0)/(60.0*60.0*23.9345), "s-1"]
    
    LatentHeatV = UNumeric[2.5e6, "J.kg-1"]
    LatentHeatS = UNumeric[2.834e6, "J.kg-1"]
    LatentHeatF = UNumeric[334e3, "J.kg-1"]

    WaterDens = UNumeric[1e3, "kg.m-3"]
    
    module Atm
      MolWtWet = UNumeric[18.0153e-3, "kg.mol-1"]
      MolWtDry = UNumeric[28.9644e-3, "kg.mol-1"]
      GasRDry = GasRUniv / MolWtDry
      GasRWat = GasRUniv / MolWtWet
      EpsV = MolWtWet/MolWtDry
      CpDry = UNumeric[1004e0, "J.kg-1.K-1"]
      CpWet = UNumeric[1810e0, "J.kg-1.K-1"]
    end

    module Ocn
      Dens0 = UNumeric[1.030e3, "kg.m-3"]
      Cp0 = UNumeric[3986e0, "J.kg-1.K-1"]
    end
    
  end
  
end
