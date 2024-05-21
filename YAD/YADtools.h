#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

#include <TFile.h>
#include <TTree.h>

// using Vec3D = ROOT::Math::XYZVector;
// using Vec4D = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;

namespace yad {

class Truth {
private:

    TFile* file;
    TTree* truth;

public:

    //Events
    size_t Event, Run, SubRun;

    //MCTruth
    size_t NPrt;

    //MCParticle
    vector<int> *PrtPdg=nullptr;
    vector<double>  *PrtStX=nullptr,
                    *PrtStY=nullptr,
                    *PrtStZ=nullptr,
                    *PrtStT=nullptr,
                    *PrtStPx=nullptr,
                    *PrtStPy=nullptr,
                    *PrtStPz=nullptr,
                    *PrtStP=nullptr,
                    *PrtStE=nullptr;

    //SimEnergy*Deposit
    size_t NDep;
    vector<int>  *DepPdg=nullptr;
    vector<double>  *DepX=nullptr,
                    *DepY=nullptr,
                    *DepZ=nullptr,
                    *DepT=nullptr,
                    *DepE=nullptr;

    Truth(const char* filename) : file{new TFile(filename)} {
        truth = file->Get<TTree>("YAD/Truth");

        //Events
        truth->SetBranchAddress("fEvent",   &Event);
        truth->SetBranchAddress("fRun",     &Run);
        truth->SetBranchAddress("fSubrun",  &SubRun);

        //MCTruth
        truth->SetBranchAddress("fNPrt",    &NPrt);

        //MCParticle
        truth->SetBranchAddress("fPrtPdg",  &PrtPdg); 
        truth->SetBranchAddress("fPrtStX",  &PrtStX); 
        truth->SetBranchAddress("fPrtStY",  &PrtStY); 
        truth->SetBranchAddress("fPrtStZ",  &PrtStZ); 
        truth->SetBranchAddress("fPrtStT",  &PrtStT); 
        truth->SetBranchAddress("fPrtStPx", &PrtStPx); 
        truth->SetBranchAddress("fPrtStPy", &PrtStPy); 
        truth->SetBranchAddress("fPrtStPz", &PrtStPz); 
        truth->SetBranchAddress("fPrtStP",  &PrtStP); 
        truth->SetBranchAddress("fPrtStE",  &PrtStE); 

        //SimEnergyDeposit
        truth->SetBranchAddress("fNDep",    &NDep);
        truth->SetBranchAddress("fDepPdg",  &DepPdg); 
        truth->SetBranchAddress("fDepX",    &DepX); 
        truth->SetBranchAddress("fDepY",    &DepY); 
        truth->SetBranchAddress("fDepZ",    &DepZ); 
        truth->SetBranchAddress("fDepT",    &DepT); 
        truth->SetBranchAddress("fDepE",    &DepE); 

    }
    ~Truth() { file-> Close(); }
    
    size_t GetEntries() { return truth->GetEntries(); }
    void GetEntry(size_t i) { truth->GetEntry(i); }
};

class Reco {
private:

    TFile* file;
    TTree* reco;

public:

    //Events
    size_t Event, Run, SubRun;

    //PFParticle
    size_t NPfp;
    vector<int> *PfpID=nullptr,
                *PfpTrkID=nullptr,
                *PfpShwID=nullptr,
                *PfpPdg=nullptr;

    //Track
    size_t NTrk;
    vector<int> *TrkID=nullptr,
                *TrkNPt=nullptr;
    vector<double>  *TrkLength=nullptr;
    vector<vector<double>>  *TrkPtX=nullptr,
                            *TrkPtY=nullptr,
                            *TrkPtZ=nullptr,
                            *TrkDirX=nullptr,
                            *TrkDirY=nullptr,
                            *TrkDirZ=nullptr;

    //Calorimetry
    vector<double>  *TrkCalRange=nullptr;

    //Shower
    size_t NShw;
    vector<int> *ShwID=nullptr;

    //Cluster
    vector<size_t> *PfpNClu=nullptr;
    vector<vector<double>>  *CluNHit=nullptr,
                            *CluPlane=nullptr,
                            *CluIntFit=nullptr,
                            *CluSumADC=nullptr,
                            *CluWidth=nullptr;

    //SpacePoint
    vector<size_t> *PfpNSpt=nullptr;
    vector<vector<double>>  *SptX=nullptr,
                            *SptY=nullptr,
                            *SptZ=nullptr;

    Reco(const char* filename) : file{new TFile(filename)} {
        reco = file->Get<TTree>("YAD/Reco");

        //Events
        reco->SetBranchAddress("fEvent",        &Event);
        reco->SetBranchAddress("fRun",          &Run);
        reco->SetBranchAddress("fSubrun",       &SubRun);

        //PFParticle
        reco->SetBranchAddress("fNPfp",         &NPfp);
        reco->SetBranchAddress("fPfpID",        &PfpID);
        reco->SetBranchAddress("fPfpTrkID",     &PfpTrkID);
        reco->SetBranchAddress("fPfpShwID",     &PfpShwID);
        reco->SetBranchAddress("fPfpPdg",       &PfpPdg);

        //Track
        reco->SetBranchAddress("fNTrk",         &NTrk);
        reco->SetBranchAddress("fTrkID",        &TrkID);
        reco->SetBranchAddress("fTrkLength",    &TrkLength);
        reco->SetBranchAddress("fTrkNPt",       &TrkNPt);
        reco->SetBranchAddress("fTrkPtX",       &TrkPtX);
        reco->SetBranchAddress("fTrkPtY",       &TrkPtY);
        reco->SetBranchAddress("fTrkPtZ",       &TrkPtZ);
        reco->SetBranchAddress("fTrkDirX",      &TrkDirX);
        reco->SetBranchAddress("fTrkDirY",      &TrkDirY);
        reco->SetBranchAddress("fTrkDirZ",      &TrkDirZ);

        //Calorimetry
        reco->SetBranchAddress("fTrkCalRange",  &TrkCalRange);

        //Shower
        reco->SetBranchAddress("fNShw",         &NShw);
        reco->SetBranchAddress("fShwID",        &ShwID);

        //Cluster
        reco->SetBranchAddress("fPfpNClu",      &PfpNClu);
        // reco->SetBranchAddress("fCluNHit",      &CluNHit);
        // reco->SetBranchAddress("fCluPlane",     &CluPlane);
        // reco->SetBranchAddress("fCluIntFit",    &CluIntFit);
        // reco->SetBranchAddress("fCluSumADC",    &CluSumADC);
        // reco->SetBranchAddress("fCluWidth",     &CluWidth);

        //SpacePoint
        reco->SetBranchAddress("fPfpNSpt",      &PfpNSpt);
        reco->SetBranchAddress("fSptX",         &SptX);
        reco->SetBranchAddress("fSptY",         &SptY);
        reco->SetBranchAddress("fSptZ",         &SptZ);

    }
    ~Reco() { file-> Close(); }
    
    size_t GetEntries() { return reco->GetEntries(); }
    void GetEntry(size_t i) { reco->GetEntry(i); }
    
};

vector<string> ReadFileList(size_t n_file, string listname)
{
  vector<string> filelist;
  string filename;
  ifstream file(listname.c_str());

  if(file){
    while(file.good()){
      file >> filename;
      filelist.push_back(filename);
      if(n_file>0 && filelist.size() == n_file) break;
    }
  }
  else{
    cerr << "\e[3mNo file named " << listname << "\e[0m" << endl;
    exit(1);
  }
  file.close();

  return filelist;
}

}