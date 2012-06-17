#ifndef __KL_CALIBRATION__H__
#define __KL_CALIBRATION__H__
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <list>
#include <vector>

#include "TH1.h"
#include "TROOT.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "gamma/Gamma.h"
#include "pi0/Pi0.h"
#include "klong/Klong.h"
#include "klong/RecKlong.h"


#define MAX_ITERACTIONS 30
#define MASS_KL 	497.648
#define MASS_PI0	134.9766

const  int N_CSI  = 2716;
const  int N_EDGE_CSI = 148;
static int edge_csi[N_EDGE_CSI] ={
  2240 , 2241 , 2242 , 2243 , 2244 ,
  2245 , 2246 , 2247 , 2248 , 2249 ,
  2250 , 2251 , 2252 , 2253 , 2254 , 
  2255 , 2264 , 2265 , 2266 , 2267 , 
  2268 , 2269 , 2270 , 2285 , 2286 ,
  2287 , 2288 , 2289 , 2308 , 2309 , 
  2310 , 2311 , 2332 , 2333 , 2334 , 
  2335 , 2358 , 2359 , 2360 , 2361 ,
  2362 , 2363 , 2364 , 2365 , 2368 ,
  2369 , 2370 , 2371 , 2376 , 2377 , 
  2378 , 2385 , 2386 , 2387 , 2394 , 
  2395 , 2396 , 2405 , 2406 , 2407 , 
  2416 , 2417 , 2418 , 2419 , 2428 , 
  2429 , 2430 , 2441 , 2442 , 2453 , 
  2454 , 2465 , 2466 , 2477 , 2478 , 
  2489 , 2490 , 2501 , 2502 , 2513 , 
  2514 , 2525 , 2526 , 2527 , 2536 , 
  2537 , 2538 , 2539 , 2548 , 2549 , 
  2550 , 2559 , 2560 , 2561 , 2568 ,
  2569 , 2570 , 2577 , 2578 , 2579 ,
  2584 , 2585 , 2586 , 2587 , 2590 ,
  2591 , 2592 , 2593 , 2594 , 2595 , 
  2596 , 2597 , 2620 , 2621 , 2622 , 
  2623 , 2644 , 2645 , 2646 , 2647 , 
  2666 , 2667 , 2668 , 2669 , 2670 , 
  2685 , 2686 , 2687 , 2688 , 2689 , 
  2690 , 2691 , 2700 , 2701 , 2702 , 
  2703 , 2704 , 2705 , 2706 , 2707 , 
  2708 , 2709 , 2710 , 2711 , 2712 ,
  2713 , 2714 , 2715}; 

int        csi_cal_factor_count[N_CSI];
double     csi_cal_factor[N_CSI];
double     csi_cal_factor_sum[N_CSI];
double     csi_cal_factor_weight[N_CSI];
double     csi_energy[N_CSI];
double     csi_time[N_CSI];
double     csi_dxy[N_CSI];

CLHEP::Hep3Vector csi_pos[N_CSI];

//double Refit_KL_full(Klong &KL);
void   user_ana_recg6_init(TH1D**,TH1D**);
double Refit_KL_freelast(Klong &KL);
int    CalEnergy_idv(Klong &KL_prefit,TH1D** ,TH1D**);

#endif
