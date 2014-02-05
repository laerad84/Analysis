#ifndef  _LOCAL_FUNCTION_H_
#define  _LOCAL_FUNCTION_H_

#include "klong/Klong.h"
#include "gamma/Gamma.h"
#include "TMath.h"

#include <list>
double const SpeedOfLight = 299.792458;

void SetGammaTime(Gamma &g);
void SetGammaTime(std::list<Gamma> glist);
double GetWeight(Gamma g);
double GetTiming(Gamma g);
double GetClusterTSigma(Gamma g);
double CalGammaTOF( Klong kl, Gamma g );

#endif //_LOCAL_FUNCTION_H_
