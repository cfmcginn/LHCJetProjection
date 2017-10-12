#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"

#include "include/getLogBins.h"
#include "include/ncollFunctions_5TeV.h"
#include "include/taaVals.h"

int jetProjection(const std::string inFileName)
{
  const Double_t rAAFactor = 0.5;

  const Int_t nPthat = 3;
  Float_t pthats[nPthat+1] = {50., 200., 500., 100000.};
  Float_t weights[nPthat] = {1., .00632675, .00632675*.00349246};

  const Int_t nCentBins = 5;
  const Double_t centBinsLow[nCentBins] = {0,5,10,30,50};
  const Double_t centBinsHigh[nCentBins] = {5,10,30,50,70};
  Double_t centNCollFrac[nCentBins];
  Double_t centJetEta2 = 0;
  Double_t centJetEtaAll = 0;
  
  double totFrac = 0;
  for(int cIter = 0; cIter < nCentBins; ++cIter){
    centNCollFrac[cIter] = findNCollFrac_Cent(centBinsLow[cIter], centBinsHigh[cIter]);
    totFrac += centNCollFrac[cIter];
  }

  double remainingFrac = 0;
  if(centBinsHigh[nCentBins-1] != 100) remainingFrac = findNCollFrac_Cent(centBinsHigh[nCentBins-1], 100);

  std::cout << "Tot frac, remaining: " << totFrac << ", " << remainingFrac << std::endl;
  
  const Double_t pbANucleons = 208;
  const Double_t xsection100Jets = 413.; //http://www.hep.ph.ic.ac.uk/~wstirlin/plots/crosssections2013.jpg
  const Double_t xsectionTotal = 7700000000.;
  //const Double_t xsection100Jets = 100.; //http://www.hep.ph.ic.ac.uk/~wstirlin/plots/crosssections2013.jpg
  const Double_t etaCut = 2.8;

  Double_t pbpbCentXSect[nCentBins];
  for(int cIter = 0; cIter < nCentBins; ++cIter){
    pbpbCentXSect[cIter] = pbANucleons*pbANucleons*xsection100Jets*centNCollFrac[cIter];
    std::cout << "XSect: " << pbpbCentXSect[cIter] << ", " << centNCollFrac[cIter] << std::endl;
  }

  TFile* outFile_p = new TFile("output/proj.root", "RECREATE");
  const Float_t jtPtLow = 100;
  const Float_t jtPtHigh = 1400;
  const Int_t nBins = 14;
  Double_t bins[nBins+1];
  getLogBins(jtPtLow, jtPtHigh, nBins, bins);

  TH1F* proj_NoWeight_p = new TH1F("proj_NoWeight_h", ";Jet p_{T};#frac{1}{N_{evt}} #frac{d^{2}N_{jet}}{dp_{T}d#eta}", nBins, bins);
  TH1F* proj_Weight_p = new TH1F("proj_Weight_h", ";Jet p_{T};#frac{1}{N_{evt}} #frac{d^{2}N_{jet}}{dp_{T}d#eta}", nBins, bins);
  proj_NoWeight_p->Sumw2();
  proj_Weight_p->Sumw2();

  TH1F* proj_Weight_Cent_p[nCentBins];
  TH1F* proj_Weight_Cent_TAA_p[nCentBins];
  
  for(Int_t iter = 0; iter < nCentBins; ++iter){
    const std::string name = "proj_Weight_Cent" + std::to_string((int)centBinsLow[iter]) + "to" + std::to_string((int)centBinsHigh[iter]) + "_h";
    const std::string nameTAA = "proj_Weight_Cent" + std::to_string((int)centBinsLow[iter]) + "to" + std::to_string((int)centBinsHigh[iter]) + "_TAA_h";
    
    proj_Weight_Cent_p[iter] = new TH1F(name.c_str(), ";Jet p_{T} [GeV/c];#frac{1}{N_{evt}} #frac{d^{2}N_{jet}}{dp_{T}d#eta} [1/(GeV/c)]", nBins, bins);
    proj_Weight_Cent_TAA_p[iter] = new TH1F(nameTAA.c_str(), ";Jet p_{T} [GeV/c];#frac{1}{#LTTAA#GT} #frac{1}{N_{evt}} #frac{d^{2}N_{jet}}{dp_{T}d#eta} [nb/(GeV/c)]", nBins, bins);

    proj_Weight_Cent_p[iter]->Sumw2();
    proj_Weight_Cent_TAA_p[iter]->Sumw2();
  }
  std::vector<std::string> fileList;
  if(inFileName.substr(inFileName.size()-5, 5).find(".root") != std::string::npos) fileList.push_back(inFileName);
  else if(inFileName.substr(inFileName.size()-4, 4).find(".txt") != std::string::npos){
    std::ifstream file(inFileName.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      if(tempStr.size() == 0) continue;
      fileList.push_back(tempStr);
    }

    file.close();
  }

  const Int_t nMaxJets = 500;
  Float_t pthat_;
  Int_t ngen_;
  Float_t genpt_[nMaxJets];
  Float_t geneta_[nMaxJets];


  for(unsigned int fileIter = 0; fileIter < fileList.size(); ++fileIter){
    TFile* inFile_p = new TFile(fileList.at(fileIter).c_str(), "READ");
    TTree* jetTree_p = (TTree*)inFile_p->Get("ak4GenJetTree");
    
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("pthat", 1);
    jetTree_p->SetBranchStatus("ngen", 1);
    jetTree_p->SetBranchStatus("genpt", 1);
    jetTree_p->SetBranchStatus("geneta", 1);
    
    jetTree_p->SetBranchAddress("pthat", &pthat_);
    jetTree_p->SetBranchAddress("ngen", &ngen_);
    jetTree_p->SetBranchAddress("genpt", genpt_);
    jetTree_p->SetBranchAddress("geneta", geneta_);
    
    const Int_t nEntries = jetTree_p->GetEntries();
    
    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(entry%1000000 == 0) std::cout << "Entry: " << entry << "/" << nEntries << std::endl;
      
      jetTree_p->GetEntry(entry);
      
      for(Int_t gI = 0; gI < ngen_; ++gI){
	if(genpt_[gI] < jtPtLow) continue;
	if(genpt_[gI] >= jtPtHigh) continue;

	++centJetEtaAll;

	if(TMath::Abs(geneta_[gI]) > etaCut) continue;
	
	++centJetEta2;

	Int_t pthatPos = -1;
	for(int pthatI = 0; pthatI < nPthat; ++pthatI){
	  if(pthats[pthatI] <= pthat_ && pthats[pthatI+1] > pthat_){
	    pthatPos = pthatI;
	    break;
	  }
	}
	if(pthatPos == -1){std::cout << "ERROR - RETURN 1" << std::endl; return 1;}
       
	proj_NoWeight_p->Fill(genpt_[gI]);
	proj_Weight_p->Fill(genpt_[gI], weights[pthatPos]);

	for(int cIter = 0; cIter < nCentBins; ++cIter){
	  proj_Weight_Cent_p[cIter]->Fill(genpt_[gI], weights[pthatPos]);
	}
      } 
    }

    inFile_p->Close();
    delete inFile_p;
  }

  std::cout << centJetEta2/centJetEtaAll << std::endl;

  for(int cIter = 0; cIter < nCentBins; ++cIter){
    pbpbCentXSect[cIter] *= centJetEta2/centJetEtaAll;
  }



  for(int cIter = 0; cIter < nCentBins; ++cIter){
    proj_Weight_Cent_p[cIter]->Scale(rAAFactor*pbpbCentXSect[cIter]*10./proj_Weight_Cent_p[cIter]->Integral());

    for(int binIter = 0; binIter < proj_Weight_Cent_p[cIter]->GetNbinsX(); ++binIter){
      if(proj_Weight_Cent_p[cIter]->GetBinContent(binIter+1) < 1){
	proj_Weight_Cent_p[cIter]->SetBinContent(binIter+1, 0);
	proj_Weight_Cent_p[cIter]->SetBinError(binIter+1, 0);
      }
      else{
	proj_Weight_Cent_p[cIter]->SetBinError(binIter+1, TMath::Sqrt(proj_Weight_Cent_p[cIter]->GetBinContent(binIter+1)));
      }
    }

    proj_Weight_Cent_p[cIter]->Scale(1./(2.*etaCut*xsectionTotal*10.*(centBinsHigh[cIter]-centBinsLow[cIter])/100));
    for(int binIter = 0; binIter < proj_Weight_Cent_p[cIter]->GetNbinsX(); ++binIter){
      proj_Weight_Cent_p[cIter]->SetBinContent(binIter+1, proj_Weight_Cent_p[cIter]->GetBinContent(binIter+1)/proj_Weight_Cent_p[cIter]->GetBinWidth(binIter+1));
      proj_Weight_Cent_p[cIter]->SetBinError(binIter+1, proj_Weight_Cent_p[cIter]->GetBinError(binIter+1)/proj_Weight_Cent_p[cIter]->GetBinWidth(binIter+1));
    }
  }

  proj_NoWeight_p->Scale(1./(2.*etaCut*xsectionTotal*10.));
  proj_Weight_p->Scale(1./(2.*etaCut*xsectionTotal));

  for(Int_t bI = 0; bI < proj_NoWeight_p->GetNbinsX(); ++bI){
    proj_NoWeight_p->SetBinContent(bI + 1, proj_NoWeight_p->GetBinContent(bI + 1)/proj_NoWeight_p->GetBinWidth(bI + 1));
    proj_Weight_p->SetBinContent(bI + 1, proj_Weight_p->GetBinContent(bI + 1)/proj_Weight_p->GetBinWidth(bI + 1));

    proj_NoWeight_p->SetBinError(bI + 1, proj_NoWeight_p->GetBinError(bI + 1)/proj_NoWeight_p->GetBinWidth(bI + 1));
    proj_Weight_p->SetBinError(bI + 1, proj_Weight_p->GetBinError(bI + 1)/proj_Weight_p->GetBinWidth(bI + 1));
  }

  outFile_p->cd();

  proj_NoWeight_p->Write("", TObject::kOverwrite);
  delete proj_NoWeight_p;

  proj_Weight_p->Write("", TObject::kOverwrite);
  delete proj_Weight_p;

  for(int cIter = 0; cIter < nCentBins; ++cIter){
    for(int binIter = 0; binIter < proj_Weight_Cent_TAA_p[cIter]->GetNbinsX(); ++binIter){
      proj_Weight_Cent_TAA_p[cIter]->SetBinContent(binIter+1, proj_Weight_Cent_p[cIter]->GetBinContent(binIter+1));
      proj_Weight_Cent_TAA_p[cIter]->SetBinError(binIter+1, proj_Weight_Cent_p[cIter]->GetBinError(binIter+1));
    }

    proj_Weight_Cent_TAA_p[cIter]->Scale(1000000./getTAA(centBinsLow[cIter], centBinsHigh[cIter]));
  }


  for(int cIter = 0; cIter < nCentBins; ++cIter){
    proj_Weight_Cent_p[cIter]->Write("", TObject::kOverwrite);
    delete  proj_Weight_Cent_p[cIter];

    proj_Weight_Cent_TAA_p[cIter]->Write("", TObject::kOverwrite);
    delete  proj_Weight_Cent_TAA_p[cIter];
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./jetProjection.exe <inFileName>" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += jetProjection(argv[1]);
  return retVal;
}
