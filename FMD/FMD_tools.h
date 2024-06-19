/*
 *  classes to store and access variables from the output ROOT file of FMD_module
 */


#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <ctime>

using namespace std;


#include <TFile.h>
#include <TTree.h>

// using Vec3D = ROOT::Math::XYZVector;
// using Vec4D = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;

namespace fmd {

class Truth {
private:

    TFile* file;
    TTree* tree;

public:

    //Events
    size_t Event, Run, SubRun;

    //MCTruth
    size_t NPrt;

    //MCParticle
    vector<int> *PrtPdg=nullptr, //[# particle]
                *PrtMomID=nullptr;
    vector<size_t>  *PrtNPt=nullptr,
                    *PrtNDau=nullptr;
    vector<vector<double>>  *PrtDauID=nullptr, //[# particle][# daughter]
                            *PrtX=nullptr, //[# particle][# point]
                            *PrtY=nullptr,
                            *PrtZ=nullptr,
                            *PrtT=nullptr,
                            *PrtPx=nullptr,
                            *PrtPy=nullptr,
                            *PrtPz=nullptr,
                            *PrtP=nullptr,
                            *PrtE=nullptr;

    //SimEnergy*Deposit
    vector<size_t>  *PrtNDep=nullptr; //[# particle]
    vector<vector<double>>  *DepPdg=nullptr, //[# particle][# deposit]
                            *DepX=nullptr,
                            *DepY=nullptr,
                            *DepZ=nullptr,
                            *DepT=nullptr,
                            *DepE=nullptr;

    Truth(const char* filename) {
        file = new TFile(filename);
        tree = file->Get<TTree>("FMD/Truth");

        //Events
        tree->SetBranchAddress("fEvent",    &Event);
        tree->SetBranchAddress("fRun",      &Run);
        tree->SetBranchAddress("fSubrun",   &SubRun);

        //MCTruth

        //MCParticle
        tree->SetBranchAddress("fNPrt",     &NPrt);
        tree->SetBranchAddress("fPrtPdg",   &PrtPdg); 
        tree->SetBranchAddress("fPrtMomID", &PrtMomID); 
        tree->SetBranchAddress("fPrtNDau",  &PrtNDau); 
        tree->SetBranchAddress("fPrtDauID", &PrtDauID); 
        tree->SetBranchAddress("fPrtNPt",   &PrtNPt); 
        tree->SetBranchAddress("fPrtX",     &PrtX); 
        tree->SetBranchAddress("fPrtY",     &PrtY); 
        tree->SetBranchAddress("fPrtZ",     &PrtZ); 
        tree->SetBranchAddress("fPrtT",     &PrtT); 
        tree->SetBranchAddress("fPrtPx",    &PrtPx); 
        tree->SetBranchAddress("fPrtPy",    &PrtPy); 
        tree->SetBranchAddress("fPrtPz",    &PrtPz); 
        tree->SetBranchAddress("fPrtP",     &PrtP); 
        tree->SetBranchAddress("fPrtE",     &PrtE); 

        //SimEnergyDeposit
        tree->SetBranchAddress("fPrtNDep",  &PrtNDep);
        tree->SetBranchAddress("fDepPdg",   &DepPdg); 
        tree->SetBranchAddress("fDepX",     &DepX); 
        tree->SetBranchAddress("fDepY",     &DepY); 
        tree->SetBranchAddress("fDepZ",     &DepZ); 
        tree->SetBranchAddress("fDepT",     &DepT); 
        tree->SetBranchAddress("fDepE",     &DepE); 
    }
    ~Truth() { file-> Close(); }
    
    size_t GetEntries() { return tree->GetEntries(); }
    void GetEntry(size_t i) { tree->GetEntry(i); }
}; //end of fmd::Truth

class Reco {
private:

    TFile* file;
    TTree* tree;

public:

    //Events
    size_t Event, Run, SubRun;

    //PFParticle
    size_t NPfp;
    vector<int> *PfpTrkID=nullptr, //[# pfparticle]
                *PfpShwID=nullptr,
                *PfpPdg=nullptr;

    //Track
    size_t NTrk;
    vector<double>  *TrkLength=nullptr; //[# track]
    vector<size_t>  *TrkNPt=nullptr; //[# track]
    vector<vector<double>>  *TrkPtX=nullptr, //[# track][# point]
                            *TrkPtY=nullptr,
                            *TrkPtZ=nullptr,
                            *TrkDirX=nullptr,
                            *TrkDirY=nullptr,
                            *TrkDirZ=nullptr;

    //Calorimetry
    vector<size_t>  *TrkCalPlane=nullptr, //[# track]
                    *TrkCalNPt=nullptr;
    vector<float>   *TrkCalRange=nullptr,
                    *TrkCalKineticE=nullptr;
    vector<vector<double>>  *TrkCaldEdx=nullptr, //[# track][# calo point]
                            *TrkCaldQdx=nullptr,
                            *TrkCalResRange=nullptr;

    //Shower
    size_t NShw;
    vector<int> *ShwID=nullptr; //[# shower]

    //Cluster
    vector<size_t> *PfpNClu=nullptr; //[# pfparticle]
    vector<vector<double>>  *CluNHit=nullptr, //[# pfparticle][# cluster]
                            *CluPlane=nullptr,
                            *CluIntegral=nullptr,
                            *CluSumADC=nullptr,
                            *CluWidth=nullptr;

    //Hit
    vector<size_t> *PfpNHit=nullptr; //[# pfparticle]
    vector<vector<double>>  *HitCluID=nullptr, //[# pfparticle][# hit]
                            *HitPlane=nullptr, 
                            *HitSumADC=nullptr,
                            *HitIntegral=nullptr,
                            *HitNSpt=nullptr;

    //SpacePoint
    vector<size_t> *PfpNSpt=nullptr; //[# pfparticle]
    vector<vector<double>>  *SptHitID=nullptr, //[# pfparticle][# spacepoint]
                            *SptX=nullptr, 
                            *SptY=nullptr,
                            *SptZ=nullptr;

    Reco(const char* filename) {
        file = new TFile(filename);
        tree = file->Get<TTree>("FMD/Reco");

        //Events
        tree->SetBranchAddress("fEvent",        &Event);
        tree->SetBranchAddress("fRun",          &Run);
        tree->SetBranchAddress("fSubrun",       &SubRun);

        //PFParticle
        tree->SetBranchAddress("fNPfp",         &NPfp);
        tree->SetBranchAddress("fPfpTrkID",     &PfpTrkID);
        tree->SetBranchAddress("fPfpShwID",     &PfpShwID);
        tree->SetBranchAddress("fPfpPdg",       &PfpPdg);

        //Track
        tree->SetBranchAddress("fNTrk",         &NTrk);
        tree->SetBranchAddress("fTrkLength",    &TrkLength);
        tree->SetBranchAddress("fTrkNPt",       &TrkNPt);
        tree->SetBranchAddress("fTrkPtX",       &TrkPtX);
        tree->SetBranchAddress("fTrkPtY",       &TrkPtY);
        tree->SetBranchAddress("fTrkPtZ",       &TrkPtZ);
        tree->SetBranchAddress("fTrkDirX",      &TrkDirX);
        tree->SetBranchAddress("fTrkDirY",      &TrkDirY);
        tree->SetBranchAddress("fTrkDirZ",      &TrkDirZ);

        //Calorimetry
        tree->SetBranchAddress("fTrkCalPlane",   &TrkCalPlane);
        tree->SetBranchAddress("fTrkCalRange",   &TrkCalRange);
        tree->SetBranchAddress("fTrkCalKineticE",&TrkCalKineticE);
        tree->SetBranchAddress("fTrkCalNPt",     &TrkCalNPt);
        tree->SetBranchAddress("fTrkCaldEdx",    &TrkCaldEdx);
        tree->SetBranchAddress("fTrkCaldQdx",    &TrkCaldQdx);
        tree->SetBranchAddress("fTrkCalResRange",&TrkCalResRange);

        //Shower
        tree->SetBranchAddress("fNShw",         &NShw);
        tree->SetBranchAddress("fShwID",        &ShwID);

        //Cluster
        tree->SetBranchAddress("fPfpNClu",      &PfpNClu);
        tree->SetBranchAddress("fCluNHit",      &CluNHit);
        tree->SetBranchAddress("fCluPlane",     &CluPlane);
        tree->SetBranchAddress("fCluIntegral",  &CluIntegral);
        tree->SetBranchAddress("fCluSumADC",    &CluSumADC);
        tree->SetBranchAddress("fCluWidth",     &CluWidth);

        //Hit
        tree->SetBranchAddress("fPfpNHit",      &PfpNHit);
        tree->SetBranchAddress("fHitCluID",     &HitCluID);
        tree->SetBranchAddress("fHitPlane",     &HitPlane);
        tree->SetBranchAddress("fHitSumADC",    &HitSumADC);
        tree->SetBranchAddress("fHitIntegral",  &HitIntegral);
        tree->SetBranchAddress("fHitNSpt",      &HitNSpt);

        //SpacePoint
        tree->SetBranchAddress("fPfpNSpt",      &PfpNSpt);
        tree->SetBranchAddress("fSptHitID",     &SptHitID);
        tree->SetBranchAddress("fSptX",         &SptX);
        tree->SetBranchAddress("fSptY",         &SptY);
        tree->SetBranchAddress("fSptZ",         &SptZ);
    }
    ~Reco() { file-> Close(); }
    
    size_t GetEntries() { return tree->GetEntries(); }
    void GetEntry(size_t i) { tree->GetEntry(i); }
}; //end of fmd::Reco

vector<string> readFileList(size_t n_file, string listname) {
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
} //end of yad::readFileList

bool isInside (vector<double> Xs,vector<double> Ys,vector<double> Zs, double Xmin=-350, double Xmax=350, double Ymin=-350, double Ymax=350, double Zmin=0, double Zmax=200) {
    bool is_in=true;
    vector<size_t> indexes = {0,Xs.size()-1};
    for (size_t i : indexes) {
        is_in = Xmin < Xs[i] && Xs[i] < Xmax &&
                Ymin < Ys[i] && Ys[i] < Ymax &&
                Zmin < Zs[i] && Zs[i] < Zmax;
        if (!is_in) break;
    }
    return is_in;
} //end of fmd::isInside

double distance(double X1, double Y1, double Z1, double X2, double Y2, double Z2) {
    return TMath::Sqrt(
        TMath::Power(X1-X2,2)
        +TMath::Power(Y1-Y2,2)
        +TMath::Power(Z1-Z2,2)
    );
} //end of fmd::distance

} //end of namespace fmd