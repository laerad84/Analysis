#include <list>
#include <vector>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"

#include "klong/Klong.h"

#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

enum{ ENERGY_CUT=0, CSI_FIDUCIAL_CUT, DISTANCE_CUT, 
      KL_VERTEX_CUT, KL_PT_CUT,KL_CHI2_CUT,KL_2NDBEST_CHI2_CUT,
      PI0_MASS_CUT, PI0_Z_CUT, SHAPE_CHI2_CUT};

double getCsiThreshold(double distance);
int kinematicCut( std::vector<Klong> const &klVec);
int shapeCut(std::vector<Klong> const &klVec);
int vetoCut(E14GNAnaDataContainer const &data);
int csiVetoCut(E14GNAnaDataContainer const &data,std::vector<Klong> const &klVec);
bool user_rec(std::list<Gamma> const &glist, std::vector<Klong> &klVec);
void user_cut(E14GNAnaDataContainer &data, std::vector<Klong> const &klVec);

