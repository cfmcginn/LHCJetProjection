#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TNamed.h"

#include "include/doGlobalDebug.h"

int extractWeightAndErr(const std::string inFileName1, const std::string inFileName2, const int thresh)
{
  const Double_t binLow = thresh-.00001;
  const Double_t binHigh = thresh+100;
  const Int_t nBins = (binHigh - binLow)/5;

  TFile* outFile_p = new TFile("extractFile.root", "RECREATE");
  TH1F* noWeightSpectra_h = new TH1F("noWeightSpectra_h", ";p_{T} Hat;Events", nBins, binLow, binHigh);
  TH1F* weightSpectra_h = new TH1F("weightSpectra_h", ";p_{T} Hat;Events", nBins, binLow, binHigh);
  TH1F* dividedSpectra_h = new TH1F("dividedSpectra_h", ";p_{T} Hat;Events", nBins, binLow, binHigh);

  std::vector<std::string> fileList;
  if(inFileName1.substr(inFileName1.size()-5, 5).find(".root") != std::string::npos) fileList.push_back(inFileName1);
  else if(inFileName1.substr(inFileName1.size()-4, 4).find(".txt") != std::string::npos){
    std::ifstream file(inFileName1.c_str());
    std::string tempStr;

    while(std::getline(file, tempStr)){
      if(tempStr.size() == 0) continue;
      fileList.push_back(tempStr);
    }

    file.close();
  }

  Float_t pthat1_;
  Float_t pthat2_;

  if(doGlobalDebug) std::cout << __LINE__ << ", " << __FILE__ << std::endl;

  for(unsigned int fileIter = 0; fileIter < fileList.size(); ++fileIter){
    TFile* inFile1_p = new TFile(fileList.at(fileIter).c_str(), "READ");
    TTree* inTree1_p = (TTree*)inFile1_p->Get("ak3GenJetTree");
    inTree1_p->SetBranchStatus("*", 0);
    inTree1_p->SetBranchStatus("pthat", 1);

    inTree1_p->SetBranchAddress("pthat", &pthat1_);

    const int nEntries1 = inTree1_p->GetEntries();

    for(int entry = 0; entry < nEntries1; ++entry){
      inTree1_p->GetEntry(entry);
      
      if(pthat1_ < binLow) continue;
      if(pthat1_ >= binHigh) continue;
      
      noWeightSpectra_h->Fill(pthat1_);
      weightSpectra_h->Fill(pthat1_);
      dividedSpectra_h->Fill(pthat1_);
    }
    inFile1_p->Close();
    delete inFile1_p;
  }

  if(doGlobalDebug) std::cout << __LINE__ << ", " << __FILE__ << std::endl;

  TFile* inFile2_p = new TFile(inFileName2.c_str(), "READ");
  TTree* inTree2_p = (TTree*)inFile2_p->Get("ak3GenJetTree");
  inTree2_p->SetBranchStatus("*", 0);
  inTree2_p->SetBranchStatus("pthat", 1);

  inTree2_p->SetBranchAddress("pthat", &pthat2_);

  const int nEntries2 = inTree2_p->GetEntries();

  for(int entry = 0; entry < nEntries2; ++entry){
    inTree2_p->GetEntry(entry);

    if(pthat2_ < binLow) continue;
    if(pthat2_ >= binHigh) continue;

    weightSpectra_h->Fill(pthat2_);
  }

  inFile2_p->Close();
  delete inFile2_p;

  if(doGlobalDebug) std::cout << __LINE__ << ", " << __FILE__ << std::endl;

  noWeightSpectra_h->Sumw2();
  weightSpectra_h->Sumw2();
  dividedSpectra_h->Sumw2();
  
  if(doGlobalDebug) std::cout << __LINE__ << ", " << __FILE__ << std::endl;

  outFile_p->cd();

  noWeightSpectra_h->Write("", TObject::kOverwrite);
  weightSpectra_h->Write("", TObject::kOverwrite);

  if(doGlobalDebug) std::cout << __LINE__ << ", " << __FILE__ << std::endl;


  dividedSpectra_h->Divide(weightSpectra_h);
  
  TF1* f1_p = new TF1("f1_p", "[0]", thresh, thresh+100);

  dividedSpectra_h->Fit("f1_p", "M N", "", thresh, thresh+100);
  dividedSpectra_h->Write("", TObject::kOverwrite);

  TNamed fitPar0("fitPar0", std::to_string(f1_p->GetParameter(0)));
  TNamed fitPar0Err("fitPar0Err", std::to_string(f1_p->GetParError(0)));
  
  fitPar0.Write("", TObject::kOverwrite);
  fitPar0Err.Write("", TObject::kOverwrite);

  if(doGlobalDebug) std::cout << __LINE__ << ", " << __FILE__ << std::endl;

  delete f1_p;

  delete noWeightSpectra_h;
  delete weightSpectra_h;
  delete dividedSpectra_h;


  outFile_p->Close();
  delete outFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 4){
    std::cout << "./Usage: ./extractWeightAndErr.exe <inFileName1> <inFileName2> <thresh>" << std::endl;
    return 1;
  }


  int retVal = 0;
  retVal += extractWeightAndErr(argv[1], argv[2], std::stoi(argv[3]));
  return retVal;
}
