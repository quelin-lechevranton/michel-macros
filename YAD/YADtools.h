/*
 *  classes to store and access variables from the output ROOT file of YAD_module
 */


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
    TTree* tree;


    vector<int> *PrtTrkID=nullptr,
                *PrtMomTrkID=nullptr,
                *DepTrkID=nullptr;
    vector<vector<double>>  *PrtDauTrkID=nullptr;
    
    vector<int> PrtMomID_value,
                DepPrtID_value;
    vector<vector<int>> PrtDauID_value;

    //return index i_prt that can be used as: PrtXxx->at(i_prt)
    int GetIDfromTrkID(int TrkID) {
        if (TrkID==0) {
            return -1;
        }
        vector<int>::iterator it = find(PrtTrkID->begin(), PrtTrkID->end(), TrkID);
        if (it==PrtTrkID->end()) {
            cerr << "i_prt not found" << endl;
            return 0;
        }
        return distance(PrtTrkID->begin(),it);
    }

public:

    //Events
    size_t Event, Run, SubRun;

    //MCTruth
    size_t NPrt;

    //MCParticle
    vector<int> *PrtPdg=nullptr,
                *PrtMomID=&PrtMomID_value;
    vector<size_t>  *PrtNPt=nullptr,
                    *PrtNDau=nullptr;
    vector<vector<int>> *PrtDauID=&PrtDauID_value;
    vector<vector<double>>  *PrtX=nullptr,
                            *PrtY=nullptr,
                            *PrtZ=nullptr,
                            *PrtT=nullptr,
                            *PrtPx=nullptr,
                            *PrtPy=nullptr,
                            *PrtPz=nullptr,
                            *PrtP=nullptr,
                            *PrtE=nullptr;

    //SimEnergy*Deposit
    size_t NDep;
    vector<int> *DepPdg=nullptr,
                *DepPrtID=&DepPrtID_value;
    vector<double>  *DepX=nullptr,
                    *DepY=nullptr,
                    *DepZ=nullptr,
                    *DepT=nullptr,
                    *DepE=nullptr;

    Truth(const char* filename) : file{new TFile(filename)} {
        tree = file->Get<TTree>("YAD/Truth");

        //Events
        tree->SetBranchAddress("fEvent",    &Event);
        tree->SetBranchAddress("fRun",      &Run);
        tree->SetBranchAddress("fSubrun",   &SubRun);

        //MCTruth
        tree->SetBranchAddress("fNPrt",     &NPrt);

        //MCParticle
        tree->SetBranchAddress("fPrtPdg",   &PrtPdg); 
        tree->SetBranchAddress("fPrtTrkID", &PrtTrkID); 
        tree->SetBranchAddress("fPrtMomTrkID",&PrtMomTrkID); 
        tree->SetBranchAddress("fPrtNDau",  &PrtNDau); 
        tree->SetBranchAddress("fPrtDauTrkID",&PrtDauTrkID); 
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
        tree->SetBranchAddress("fNDep",     &NDep);
        tree->SetBranchAddress("fDepPdg",   &DepPdg); 
        tree->SetBranchAddress("fDepTrkID", &DepTrkID); 
        tree->SetBranchAddress("fDepX",     &DepX); 
        tree->SetBranchAddress("fDepY",     &DepY); 
        tree->SetBranchAddress("fDepZ",     &DepZ); 
        tree->SetBranchAddress("fDepT",     &DepT); 
        tree->SetBranchAddress("fDepE",     &DepE); 
    }
    ~Truth() { file-> Close(); }
    
    size_t GetEntries() { return tree->GetEntries(); }
    void GetEntry(size_t i) { 
        tree->GetEntry(i); 

        for (size_t i_prt=0; i_prt < NPrt; i_prt++) {

            PrtMomID_value.push_back(GetIDfromTrkID(PrtMomTrkID->at(i_prt)));

            vector<int> tpPrtDauID_value;
            for (size_t i_dau=0; i_dau < PrtNDau->at(i_prt); i_dau++) {

                tpPrtDauID_value.push_back(GetIDfromTrkID((*PrtDauTrkID)[i_prt][i_dau]));
            }
            PrtDauID_value.push_back(tpPrtDauID_value);
        }

        for (size_t i_dep=0; i_dep < NDep; i_dep++) {
            DepPrtID_value.push_back(GetIDfromTrkID(DepTrkID->at(i_dep)));
        }
    }
};

class Reco {
private:

    TFile* file;
    TTree* tree;

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

    Reco(const char* filename) : file{new TFile(filename)} {
        tree = file->Get<TTree>("YAD/Reco");

        //Events
        tree->SetBranchAddress("fEvent",        &Event);
        tree->SetBranchAddress("fRun",          &Run);
        tree->SetBranchAddress("fSubrun",       &SubRun);

        //PFParticle
        tree->SetBranchAddress("fNPfp",         &NPfp);
        tree->SetBranchAddress("fPfpID",        &PfpID);
        tree->SetBranchAddress("fPfpTrkID",     &PfpTrkID);
        tree->SetBranchAddress("fPfpShwID",     &PfpShwID);
        tree->SetBranchAddress("fPfpPdg",       &PfpPdg);

        //Track
        tree->SetBranchAddress("fNTrk",         &NTrk);
        tree->SetBranchAddress("fTrkID",        &TrkID);
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