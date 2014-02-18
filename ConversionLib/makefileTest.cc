//============================================================================
// Name        : makefileTest.cpp
// Author      : Lee Jong-won
// Version     :
// Copyright   : LAERAD
// Description : Hello World in C, Ansi-style
//============================================================================

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TGraph.h>

#include <GsimData/GsimDetectorEventData.h>
#include <GsimData/GsimDigiData.h>
#include "include/ANAconvsimCsI.h"
#include "include/ConvertDigiData.h"
#include "include/ConvertData.h"

int main( int argc, char** argv){
	if( argc != 3){
		puts("./main <input ROOT file> <output ROOT file>");
		return -1;
	}

	puts("!!!Hello World!!!");
	ConvertData* data = new ConvertData(argv[1],argv[2]);
	
	data->Add("CSI");
	data->Add("CC03");
	data->Add("CV");
	data->Add("OEV");
	data->Add("Laser");
	data->Add("Cosmic");
	data->Add("Sci");

	if(!data->AddTrack()){
	  std::cout<< "Add Track is Aborted" << std::endl;
	}	

	std::cout << data->GetNDetector() << std::endl;
	std::cout << "ConvertAll" << std::endl;
	Long_t nentries = data->ConvertAll();
	std::cout << nentries << std::endl;
	data->~ConvertData();

}
