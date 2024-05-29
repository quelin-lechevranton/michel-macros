/*
 *  classes to store and access variables from the output ROOT file of YAD_module
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

namespace yad {

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
    vector<int> *PrtPdg=nullptr,
                *PrtMomID=nullptr;
    vector<size_t>  *PrtNPt=nullptr,
                    *PrtNDau=nullptr;
    vector<vector<double>>  *PrtDauID=nullptr,
                            *PrtX=nullptr,
                            *PrtY=nullptr,
                            *PrtZ=nullptr,
                            *PrtT=nullptr,
                            *PrtPx=nullptr,
                            *PrtPy=nullptr,
                            *PrtPz=nullptr,
                            *PrtP=nullptr,
                            *PrtE=nullptr;

    //SimEnergy*Deposit
    vector<size_t>  *PrtNDep=nullptr;
    vector<vector<double>>  *DepPdg=nullptr,
                            *DepX=nullptr,
                            *DepY=nullptr,
                            *DepZ=nullptr,
                            *DepT=nullptr,
                            *DepE=nullptr;

    Truth(const char* filename) {
        file = new TFile(filename);
        tree = file->Get<TTree>("YAD/Truth");

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
}; //end of yad::Truth

class Reco {
private:

    TFile* file;
    TTree* tree;

public:

    //Events
    size_t Event, Run, SubRun;

    //PFParticle
    size_t NPfp;
    vector<int> *PfpTrkID=nullptr,
                *PfpShwID=nullptr,
                *PfpPdg=nullptr;

    //Track
    size_t NTrk;
    vector<double>  *TrkLength=nullptr;
    vector<size_t>  *TrkNPt=nullptr;
    vector<vector<double>>  *TrkPtX=nullptr,
                            *TrkPtY=nullptr,
                            *TrkPtZ=nullptr,
                            *TrkDirX=nullptr,
                            *TrkDirY=nullptr,
                            *TrkDirZ=nullptr;

    //Calorimetry
    vector<size_t>  *TrkCalPlane=nullptr,
                    *TrkCalNPt=nullptr;
    vector<float>   *TrkCalRange=nullptr,
                    *TrkCalKineticE=nullptr;
    vector<vector<double>>  *TrkCaldEdx=nullptr,
                            *TrkCaldQdx=nullptr,
                            *TrkCalResRange=nullptr;

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

    Reco(const char* filename) {
        file = new TFile(filename);
        tree = file->Get<TTree>("YAD/Reco");

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
        tree->SetBranchAddress("fCluIntFit",    &CluIntFit);
        tree->SetBranchAddress("fCluSumADC",    &CluSumADC);
        tree->SetBranchAddress("fCluWidth",     &CluWidth);

        //SpacePoint
        tree->SetBranchAddress("fPfpNSpt",      &PfpNSpt);
        tree->SetBranchAddress("fSptX",         &SptX);
        tree->SetBranchAddress("fSptY",         &SptY);
        tree->SetBranchAddress("fSptZ",         &SptZ);

    }
    ~Reco() { file-> Close(); }
    
    size_t GetEntries() { return tree->GetEntries(); }
    void GetEntry(size_t i) { tree->GetEntry(i); }
}; //end of yad::Reco

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

bool isInside (vector<double> Xs,vector<double> Ys,vector<double> Zs, double Xmin=-350, double Xmax=350, double Ymin=-350, double Ymax=350, double Zmin=0, double Zmax=300) {
    bool is_in=true;
    vector<size_t> indexes = {0,Xs.size()-1};
    for (size_t i : indexes) {
        is_in = Xmin < Xs[i] && Xs[i] < Xmax &&
                Ymin < Ys[i] && Ys[i] < Ymax &&
                Zmin < Zs[i] && Zs[i] < Zmax;
        if (!is_in) break;
    }
    return is_in;
} //end of yad::isInside

double distance(double X1, double Y1, double Z1, double X2, double Y2, double Z2) {
    return TMath::Sqrt(
        TMath::Power(X1-X2,2)
        +TMath::Power(Y1-Y2,2)
        +TMath::Power(Z1-Z2,2)
    );
} //end of yad::distance

} //end of namespace yad