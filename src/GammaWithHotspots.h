/*
 * Authors: Gustavo V. Barroso and Julien Y. Dutheil
 * Created: 16/04/2018
 * Last modified: 16/04/2018
 *
 */


#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Numeric/Constraints.h>

class GammaWithHotspots: //Discretized Gamma + Hotspot intensity as an extra parameter
  public bpp::GammaDiscreteDistribution { 

public:
  GammaWithHotspots(unsigned int numGammaCategories,
                    double initialShape,
                    double initialRate,
                    double initialHeat,
                    const bpp::Vdouble& heatConstraint):
    //AbstractParameterAliasable("Gamma."), //Not needed with bio++ 2.4.0, uncomment if older version.
    GammaDiscreteDistribution(numGammaCategories, initialShape, initialRate)
    { 
      //hotspot intensity
      addParameter_(new bpp::Parameter("Gamma.heat", initialHeat, new bpp::IntervalConstraint(heatConstraint[0], heatConstraint[1], true, true)));    
    }

};
