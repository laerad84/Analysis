//
//   Example of a program to fit non-equidistant data points
//   =======================================================
//
//   The fitting function fcn is a simple chisquare function
//   The data consists of 5 data points (arrays x,y,z) + the errors in errorsz
//   More details on the various functions or parameters for these functions
//   can be obtained in an interactive ROOT session with:
//    Root > TMinuit *minuit = new TMinuit(10);
//    Root > minuit->mnhelp("*")  to see the list of possible keywords
//    Root > minuit->mnhelp("SET") explains most parameters
//Author: Rene Brun
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cmath>
#include <cassert>


#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2C.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TMarker.h"
#include "TPad.h"
#include "TMath.h"
#include "TMinuit.h"
//#include "IDHandler.cpp"


//const Int_t N = 2716;
const Int_t N = 20;
static Double_t *offsets;
Double_t (*means)[N];
Double_t (*errors)[N];

Double_t chisq = 0.;
Double_t oldChi2 = 1E10;
Int_t cycle=0;

//IDHandler* handler = new IDHandler( "../crystal.txt", "../FADC.txt");


/*
bool idcut( int region, int id )
{
	
	double x=0, y=0;
	handler->GetMetricPosition( id, x, y );
	
	if( region==0 ){  // All region
		return true;
	}
	
	if( region < 7 ){ // Small bottom region
		if( -600+200*(region-1)<x && x<-400+200*(region-1) && y<-100 ) return true;
		return false;
	}
	
	if( region < 13 ){ // Small top region
		if( -600+200*(region-7)<x && x<-400+200*(region-7) && -100<y ) return true;
		return false;
	}
	
	if( region==13 ){  // Large left region
		if( x<-600 ) return true;
		return false;
	}
	
	if( region==14 ){  // Large right region
		if( 600<x ) return true;
		return false;
	}
		
		
	return false;
}
*/




void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  std::cout<< "fcn " << std::endl;
	//calculate chisquare
	
	Double_t delta = 0.;
	chisq = 0.0;
	for ( int index=0; index<N; index++) {
		for( int jndex=index+1; jndex<N; jndex++ ){
			if( errors[index][jndex]!=0 ){
				delta  = ( means[index][jndex] - par[index] + par[jndex] ) / errors[index][jndex];
				chisq += delta*delta;
			}
		}
	}
						   
	f = chisq;

	
	if( cycle%100 == 0 ){
		cerr << "Iteration:" << cycle << "\tchisquare: " << f << "\t";
		for( int index=0; index<N; index++ ){
			if( index%500==0 ){
				cerr << par[index] << " ";
			}
		}
		cerr << endl;
		oldChi2 = f;
	}
	
	cycle++;
	std::cout<< "end fcn" << std::endl;
}




/*
//______________________________________________________________________________
void cosmicAnalysis3( const char* kind, const int region )
{
	
	
	// Initialize
	offsets = new Double_t [N];
	means   = new Double_t [N][N];
	errors  = new Double_t [N][N];
	
	bool idFlag[N] = {false};
	
	for( int index=0; index<N; index++ ){
		offsets[index] = 0.0;
		for( int jndex=0; jndex<N; jndex++ ){
			means[index][jndex]  = 0.0;
			errors[index][jndex] = 0.0;
		}
	}
	
	
	// Get data points
	TFile* tf = new TFile("timing.root");
	TTree* tr = (TTree*)tf->Get("tr");
	const long numberOfEntries = tr->GetEntries();
	
	Int_t firstID=0, secondID=0, numberOfEvents=0, ndf=0;
	Double_t mean=0.0, meanError=0.0, chi2=0.0;
	{
		tr->SetBranchAddress( "firstID", &firstID );
		tr->SetBranchAddress( "secondID", &secondID );
		tr->SetBranchAddress( Form( "%sMean", kind), &mean );
		tr->SetBranchAddress( Form( "%sMeanErr", kind), &meanError );
		tr->SetBranchAddress( Form( "%sChi2", kind), &chi2 );
		tr->SetBranchAddress( Form( "%sNdf", kind), &ndf );
		tr->SetBranchAddress( "numberOfEvents", &numberOfEvents );
	}
	for( long index=0; index<numberOfEntries; index++ ){
		tr->GetEntry( index );
		if( numberOfEvents>0 && meanError!=1.41421 && ndf>10 && abs(chi2/ndf-1.0)<0.5 && meanError<5. && idcut(region, firstID) && idcut(region, secondID) ){
			means[firstID][secondID]  = mean;
			errors[firstID][secondID] = meanError;
			idFlag[firstID]  = true;
			idFlag[secondID] = true;
		}
	}
	tf->Close();
	delete tf;
	
	
	// Get initial values for iteration
	if( region==0 ){
		for( int index=0; index<N; index++ ){
			if( idFlag[index] ){
				if( index < 2240 ){
					offsets[index] += 14./1196.*139.;
				}else{
					offsets[index] -= 14./1196.*1057.;
				}
			}
		}
	}
	
	
	
	
	//getchar();
	
	TMinuit* gMinuit = new TMinuit(N);
	gMinuit->SetFCN( fcn );
	gMinuit->SetPrintLevel( -1 );
	
	
	
	// Set starting values and step sizes for parameters
	gMinuit->Command( "SET ERR" );

	double minValue=-10.0, maxValue=10.0;
	for( int index=0; index<N; index++ ){
		if( idFlag[index] ){
			if( region ){
				if( gMinuit->DefineParameter( index, Form("p%d",index), offsets[index], 0.1, minValue, maxValue ) ){
					cerr << "DEBUG : index=" << index << endl;
					return;
				}
			}else{
				if( index < 2240 ){
					if( gMinuit->DefineParameter( index, Form("p%d",index), offsets[index], 0.1, minValue, maxValue ) ){
						cerr << "DEBUG : index=" << index << endl;
						return;
					}
				}else{
					if( gMinuit->DefineParameter( index, Form("p%d",index), offsets[index], 0.1, minValue-14., maxValue-14. ) ){
						cerr << "DEBUG : index=" << index << endl;
						return;
					}
				}
			}
		}else{
			if( gMinuit->DefineParameter( index, Form("p%d",index), 0, 1., -2., 2. ) ){
				cerr << "DEBUG : index=" << index << endl;
				return;
			}
			gMinuit->FixParameter( index );
		}
	}

	
	gMinuit->SetMaxIterations(0x7FFFFFFF);
	gMinuit->Migrad();
	
	cerr << "Iteration:" << cycle << "\tchisquare: " << chisq << endl;
	cycle = 0;
	
	ofstream output( Form("%sOffset.txt", kind) );
	double par=0.0, err=0.0;
	for( int index=0; index<N; index++ ){
		gMinuit->GetParameter( index, par, err );
		output << index << "\t" << par << "\t" << err << endl;
	}
	
	

	
	
	delete   gMinuit;
	delete[] offsets;
	delete[] means;
	delete[] errors;

	return;

}
*/



 /*
int connection( void )
{
	
	const int numberOfCrystals = 2716;
	int crystal[numberOfCrystals] = {0};
	double buf[2] = {0.0};
	double offset[numberOfCrystals] = {0.0};
	double error[numberOfCrystals] = {0.0};
	
	
	ofstream output("timeOffset0.txt");
	
	for( int region=1; region<15; region++ ){
		ifstream input( Form( "Region%d.txt", region ) );
		while( input >> crystal[0] >> buf[0] >> buf[1] ){
			if( buf[1] != 0 ){
				offset[crystal[0]] = buf[0];
				error[crystal[0]]  = buf[1];
			}
		}
		input.close();
	}
	
	for( int index=0; index<numberOfCrystals; index++ ){
		crystal[index] = index;
		output << crystal[index] << "\t" << offset[index] << "\t" << error[index] << endl;
	}
	
	output.close();
	
	return 0;
}
*/	
	
	
void testMinuit(){
  std::cout<< __LINE__ << std::endl;
  double Testeans[3][3]  = {{0}};
  double TestErrors[3][3] = {{0}};
  double TestOffsets[3]   = {0};
  std::cout<< __LINE__ << std::endl;
  offsets = new Double_t [N];
  means   = new Double_t [N][N];
  errors  = new Double_t [N][N];
  for( int iIndex  =0; iIndex < N; iIndex++){
    offsets[iIndex] = 0; 
    for( int jIndex = 0; jIndex < N; jIndex++){
      means[iIndex][jIndex] = 0;
      errors[iIndex][jIndex] = 0;
    }
  }
  std::ifstream ifs("Data.txt");
  double mean;
  double error; 
  double entries;
  int iID;
  int jID;
  while( ifs >> iID >> jID >> entries >> mean >> error ){
    means[iID][jID]= mean;
    errors[iID][jID] = error;
  }
  
  std::cout<< __LINE__ << std::endl;
  
  //  N = 3;
  TMinuit *Minuit = new TMinuit(N);
  gMinuit->SetFCN( fcn );
  gMinuit->SetPrintLevel( -1 );
  gMinuit->Command( "SET ERR"  );
  double minValue= -50;
  double maxValue=  50;
  std::cout<< __LINE__ << std::endl;

  for ( int index = 0; index < N ;index ++){
  std::cout<< __LINE__ << std::endl;
    gMinuit->DefineParameter( index ,Form("p%d", index ), offsets[index], 0.1, minValue, maxValue );
  }
  std::cout<< __LINE__ << std::endl;

  gMinuit->SetMaxIterations( 0x7FFFFFFF );
  gMinuit->Migrad();

  std::cerr << "Iteration :" << cycle << "\tChisquare : " << chisq << std::endl;
  cycle = 0; 
  ofstream output ( Form("%sOffset.txt","test"));
  double par=0.0, err= 0.0;
  for( int index = 0; index < N; index++){
    gMinuit->GetParameter(index, par, err );
    output << index << "\t" << par << "\t" << err << std::endl; 
  }


  delete gMinuit; 
  output.close();
  delete[] offsets;
  delete[] means;
  delete[] errors;
}
	



