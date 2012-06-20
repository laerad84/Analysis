#include "E14ConvWriter.h"

E14ConvWriter::E14ConvWriter(char* SumupFile,char* outFile){  
  m_mapFilename = SumupFile;
  m_outFilename = outFile;
  m_tf = new TFile( m_outFilename.c_str(), "recreate");
  InitMap();
}

E14Convwriter::~E14ConvWriter(){
  ;
}

bool E14ConvWriter::InitMap(){
  map = new E14MapReader(m_mapFilename.c_str());
  map->Add("Csi");
  map->Add("CC03");
  map->Add("OEV");
  map->Add("Cosmic");
  map->Add("Laser");
  map->Add("CV");
  map->Add("Etc");
  map->SetMap();
  for( int iMod =0 ; iMod< map->GetNModule(); iMod++){
    map->CopyMap( iMod, ModMap[iMod]);
  }

  std::cout<< "Setup Map is Done" << std::endl;
  std::cout<< map->GetNmodule() << std::endl;			 
  return true;
}

bool E14ConvWriter::Branch(char* ModID){
  m_tr->Branch(Form("%s.nFit", &
}
