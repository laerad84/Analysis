/*
 * ANAconvsimCsI.cc
 *
 *  Created on: Apr 1, 2011
 *      Author: jwlee
 */
#ifndef ANACONVSIMCSI_H_
#include <ANAconvsimCsI.h>
#endif

ANAconvSim::ANAconvSim(char* ifileName,char* ofileName,char* ofileoption = "recreate")
{
	this->m_itfile    = new TFile(ifileName);
	this->m_eventTree = (TTree*)this->m_itfile->Get("eventTree00");
	this->m_otfile    = new TFile(ofileName,ofileoption);
	this->m_outTree   = new TTree("Tree","");
	this->m_detectorData = new GsimDetectorData();
	this->m_eventData    = new GsimDetectorEventData();
}

ANAconvSim::~ANAconvSim(){
}
Bool_t ANAconvSim::SetInputFile(char* ifileName){
	this->m_itfile    = new TFile(ifileName);
	this->m_eventTree = (TTree*)this->m_itfile->Get("eventTree00");
	return true;
}
Bool_t ANAconvSim::SetOutputFile(char* ofileName, char* ofileoption = "recreate"){
	this->m_otfile    = new TFile(ofileName,ofileoption);
	this->m_outTree   = new TTree("Tree","");
	return true;
}

Bool_t ANAconvSim::Branch(){
	m_outTree->Branch("eventNum",&eventNum,"eventNum/I");
	m_outTree->Branch("nDigi",&nDigi,"nDigi/I");
	m_outTree->Branch("totalE",&totalE,"totalE/D");
	m_outTree->Branch("ID",ID,"ID[nDigi]/I");//nDigi
	m_outTree->Branch("depE",depE,"depE[nDigi]/D");//nDigi

	return true;
}

Bool_t ANAconvSim::SetBranchAddress(){
	std::cout << m_itfile << std::endl;
	std::cout << m_otfile << std::endl;
	std::cout << m_eventTree << std::endl;
	std::cout << m_outTree << std::endl;


	this->m_eventTree->SetBranchStatus("*",0);
	this->m_eventTree->SetBranchStatus("CSI.nDigi");
	this->m_eventTree->SetBranchStatus("CSI.digi.modID");
	this->m_eventTree->SetBranchStatus("CSI.digi.energy");
	this->m_eventTree->SetBranchAddress("CSI.",&(this->m_eventData));
	return true;
}

Long_t ANAconvSim::Convert(){

	this->SetBranchAddress();
	this->Branch();

	Long_t nEvent = this->m_eventTree->GetEntries();
	TClonesArray* csiarr;
	for( Long_t ievent = 0; ievent<nEvent; ievent++){
		this->m_eventTree->GetEntry(ievent);
		this->eventNum = ievent;
		this->nDigi    = this->m_eventData->nDigi;
		this->totalE   = 0;
		csiarr         = this->m_eventData->digi;
		for( Int_t idigi = 0; idigi < this->nDigi; idigi++){
			GsimDigiData* digi = (GsimDigiData*)csiarr->UncheckedAt(idigi);
			this->ID[idigi]    = digi->modID;
			this->depE[idigi]  = digi->energy;
			this->totalE      += this->depE[idigi];
			//std::cout   << digi->energy      << " : "
			//			<< this->depE[idigi] << " : "
			//			<< this->totalE      << std::endl;
		}
		//std::cout << "TOTAL " << this->totalE << std::endl;
		this->m_outTree->Fill();
	}
	this->m_outTree->Write();
	this->m_otfile->Close();
	return nEvent;
}
