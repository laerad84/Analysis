#include <iostream>
#include <list>
#include <vector>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"

#include "klong/Klong.h"
#include "klong/RecKlong.h"
#include "pi0/Pi0.h"
#include "cluster/ClusterFinder.h"
#include "gamma/GammaFinder.h"
#include "rec2g/Rec2g.h"

#include "csimap/CsiMap.h"

#include "TMath.h"

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
double rec_mass2g( std::list<Gamma> const & glist, double recPosition);
void recVtxWithConstM( const Gamma& g1, const Gamma& g2, double Mass, double* recZm ,double* recZsig2 );
std::list<Pi0> recPi0withConstM( std::list<Gamma> glist, double mass );
bool user_rec( std::list<Gamma> const &glist, std::list<Pi0>& piList, double &mass ,double recPosition);
void correctTime( Gamma& g1);
