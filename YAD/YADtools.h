#include <vector>
#include <string>
#include <iostream>

using namespace std;

#include <TFile.h>
#include <TTree.h>
// #include <Vector3D.h>
// #include "Math/GenVector/LorentzVector.h"

using Vec3D = ROOT::Math::XYZVector;
using Vec4D = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;

namespace yad {

class Truth {
private:

    TFile* file;
    TTree* truth;


    vector<double>  *PrtStX=nullptr,
                    *PrtStY=nullptr,
                    *PrtStZ=nullptr,
                    *PrtStT=nullptr,
                    *PrtStPx=nullptr,
                    *PrtStPy=nullptr,
                    *PrtStPz=nullptr,
                    *PrtStE=nullptr;

    vector<double>  *DepX=nullptr,
                    *DepY=nullptr,
                    *DepZ=nullptr,
                    *DepT=nullptr;

public:

    //Events
    int Event, Run, SubRun;

    //MCTruth
    int NPrt;

    //MCParticle
    vector<int> *PrtPdg=nullptr;
    vector<Vec4D>   PrtStPt,
                    PrtStP;

    //SimEnergy*Deposit
    int NDep;
    vector<int>  *DepPdg=nullptr;
    vector<Vec4D>   Dep;
    vector<double>  *DepE=nullptr;

    Truth(const char* filename) : file{new TFile(filename)} {
        truth = file->Get<TTree>("YAD/Truth");

        //Events
        truth->SetBranchAddress("fEvent",          &Event);
        truth->SetBranchAddress("fRun",            &Run);
        truth->SetBranchAddress("fSubrun",         &SubRun);

        //MCTruth
        truth->SetBranchAddress("fNPrt",   &NPrt);

        //MCParticle
        truth->SetBranchAddress("fPrtPdg",      &PrtPdg); 
        truth->SetBranchAddress("fPrtStX",   &PrtStX); 
        truth->SetBranchAddress("fPrtStY",   &PrtStY); 
        truth->SetBranchAddress("fPrtStZ",   &PrtStZ); 
        truth->SetBranchAddress("fPrtStT",   &PrtStT); 
        truth->SetBranchAddress("fPrtStPx",  &PrtStPx); 
        truth->SetBranchAddress("fPrtStPy",  &PrtStPy); 
        truth->SetBranchAddress("fPrtStPz",  &PrtStPz); 
        truth->SetBranchAddress("fPrtStE",   &PrtStE); 

        //SimEnergyDeposit
        truth->SetBranchAddress("fNDep",    &NDep);
        truth->SetBranchAddress("fDepPdg",       &DepPdg); 
        truth->SetBranchAddress("fDepX",         &DepX); 
        truth->SetBranchAddress("fDepY",         &DepY); 
        truth->SetBranchAddress("fDepZ",         &DepZ); 
        truth->SetBranchAddress("fDepT",         &DepT); 
        truth->SetBranchAddress("fDepE",         &DepE); 

    }
    ~Truth() { file-> Close(); }
    
    int GetEntries() { return truth->GetEntries(); }
    void GetEntry(int i) { 
        
        truth->GetEntry(i); 

        for(int i_prt=0; i_prt<NPrt; i_prt++) {

            PrtStPt.emplace_back( 
                PrtStX->at(i_prt),
                PrtStY->at(i_prt),
                PrtStZ->at(i_prt),
                PrtStT->at(i_prt)
            );
            PrtStP.emplace_back(
                PrtStPx->at(i_prt),
                PrtStPy->at(i_prt),
                PrtStPz->at(i_prt),
                PrtStE->at(i_prt)
            );
        }

        for(int i_dep=0; i_dep<NDep; i_dep++) {
            Dep.emplace_back(
                DepX->at(i_dep),
                DepY->at(i_dep),
                DepZ->at(i_dep),
                DepT->at(i_dep)
            );
        }
    }
};

class Reco {
private:

    TFile* file;
    TTree* reco;

    vector<vector<double>>  *TrkPtX=nullptr,
                            *TrkPtY=nullptr,
                            *TrkPtZ=nullptr,
                            *TrkDirX=nullptr,
                            *TrkDirY=nullptr,
                            *TrkDirZ=nullptr;

    vector<vector<double>>  *SptX=nullptr,
                            *SptY=nullptr,
                            *SptZ=nullptr;

public:

    //Events
    int Event, Run, SubRun;

    //PFParticle
    int NPfp;
    vector<int> *PfpTrkID=nullptr,
                *PfpShwID=nullptr,
                *PfpID=nullptr,
                *PfpPdg=nullptr;

    //Track
    int NTrk;
    vector<int> *TrkID=nullptr,
                *TrkNPt=nullptr;
    vector<double>  *TrkLength=nullptr;
    vector<vector<Vec3D>>   TrkPt,
                            TrkDir;

    //Calorimetry
    vector<double>  *TrkCalRange=nullptr;

    //Shower
    int NShw;
    vector<int> *ShwID=nullptr;

    //Cluster
    vector<int> *PfpNClu=nullptr;
    vector<vector<int>>     *CluNHit=nullptr,
                            *CluPlane=nullptr;
    vector<vector<double>>  *CluIntFit=nullptr,
                            *CluSumADC=nullptr,
                            *CluWidth=nullptr;

    //SpacePoint
    vector<int> *PfpNSpt=nullptr;
    vector<vector<Vec3D>>   Spt;

    Reco(const char* filename) : file{new TFile(filename)} {
        reco = file->Get<TTree>("YAD/Reco");

        //Events
        reco->SetBranchAddress("fEvent",          &Event);
        reco->SetBranchAddress("fRun",            &Run);
        reco->SetBranchAddress("fSubrun",         &SubRun);

        //PFParticle
        reco->SetBranchAddress("fNPfp",  &NPfp);
        reco->SetBranchAddress("fPfpTrkID",    &PfpTrkID);
        reco->SetBranchAddress("fPfpShwID",     &PfpShwID);
        reco->SetBranchAddress("fPfpID",         &PfpID);
        reco->SetBranchAddress("fPfpPdg",        &PfpPdg);

        //Track
        reco->SetBranchAddress("fNTrk",       &NTrk);
        reco->SetBranchAddress("fTrkID",         &TrkID);
        reco->SetBranchAddress("fTrkLength",     &TrkLength);
        reco->SetBranchAddress("fTrkNPt",    &TrkNPt);
        reco->SetBranchAddress("fTrkPtX",        &TrkPtX);
        reco->SetBranchAddress("fTrkPtY",        &TrkPtY);
        reco->SetBranchAddress("fTrkPtZ",        &TrkPtZ);
        reco->SetBranchAddress("fTrkDirX",       &TrkDirX);
        reco->SetBranchAddress("fTrkDirY",       &TrkDirY);
        reco->SetBranchAddress("fTrkDirZ",       &TrkDirZ);

        //Calorimetry
        reco->SetBranchAddress("fTrkCalRange",     &TrkCalRange);

        //Shower
        reco->SetBranchAddress("fNShw",      &NShw);
        reco->SetBranchAddress("fShwID",        &ShwID);

        //Cluster
        reco->SetBranchAddress("fPfpNClu",     &PfpNClu);
        // reco->SetBranchAddress("fCluNHit",      &CluNHit);
        reco->SetBranchAddress("fCluIntFit",     &CluIntFit);
        reco->SetBranchAddress("fCluSumADC",     &CluSumADC);
        reco->SetBranchAddress("fCluWidth",      &CluWidth);
        // reco->SetBranchAddress("fCluPlane",    &CluPlane);

        //SpacePoint
        reco->SetBranchAddress("fPfpNSpt",       &PfpNSpt);
        reco->SetBranchAddress("fSptX",        &SptX);
        reco->SetBranchAddress("fSptY",        &SptY);
        reco->SetBranchAddress("fSptZ",        &SptZ);

    }
    ~Reco() { file-> Close(); }
    
    int GetEntries() { return reco->GetEntries(); }
    void GetEntry(int i) { 

        reco->GetEntry(i); 

        for (int i_trk=0; i_trk<NTrk; i_trk++) {

            vector<Vec3D> tpTrkPt, tpTrkDir;
            for (size_t i_tpt=0; i_tpt<TrkPtX->at(i_trk).size(); i_tpt++) {

                tpTrkPt.emplace_back(
                    TrkPtX->at(i_trk)[i_tpt],
                    TrkPtY->at(i_trk)[i_tpt],
                    TrkPtZ->at(i_trk)[i_tpt]
                );

                tpTrkDir.emplace_back(
                    TrkDirX->at(i_trk)[i_tpt],
                    TrkDirY->at(i_trk)[i_tpt],
                    TrkDirZ->at(i_trk)[i_tpt]
                );
            }

            TrkPt.push_back(tpTrkPt);
            TrkDir.push_back(tpTrkDir);

        }
        
        for (int i_pfp=0; i_pfp<NPfp; i_pfp++) {

            vector<Vec3D> tpSpt;
            for (int i_spt=0; i_spt<PfpNSpt->at(i_pfp); i_spt++) {
                
                tpSpt.emplace_back(
                    SptX->at(i_pfp)[i_spt],
                    SptY->at(i_pfp)[i_spt],
                    SptZ->at(i_pfp)[i_spt]
                );
            }

            Spt.push_back(tpSpt);
        }

    }
    
};

vector<string> ReadFileList(int nfiles, string FileName)
{
  vector<string> vec;
  string file;
  string filename = FileName;
  ifstream File(filename.c_str());

  if(File){
    cout << "    ReadFileList - " << filename << " opened"<< endl;
    while(File.good()){
      File >> file;
      vec.push_back(file);
      if(nfiles>0 && vec.size() == nfiles) break;
    }
  }
  else{
    cerr << "    ReadFileList - " << filename << " not found. Exit." << endl;
    exit(1);
  }

  File.close();

  if(nfiles == 0) vec.pop_back();
  cout << "    ReadFileList - " << vec.size() << " files found in " << FileName << endl;
  cout << "    finished" << endl;

  return vec;
}

}