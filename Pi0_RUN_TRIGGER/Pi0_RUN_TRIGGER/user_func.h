#include "TFile.h"
#include "TChain.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/GammaFinder.h"
#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "TApplication.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <list>
#include "pi0/Pi0.h"
#include "gamma/Gamma.h"

double rec_mass2g(std::list<Gamma> const &glist, double recPosition);
void recVtxWithConstM( const Gamma& g1, const Gamma& g2, double Mass,
		       double* recZ, double* recZsig2 );
std::list<Pi0> recPi0withConstM( std::list<Gamma> glist, double mass );
bool user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList,double &mass , double recPosition);
