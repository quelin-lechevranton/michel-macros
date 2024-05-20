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

namespace YAD {

class Reco {
private:

    TFile* file;
    TTree* reco;
public:

    //Events
    int Event, Run, SubRun;

    //MCTruth
    int NPrt;

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
    int NDep;
    vector<double>  *DepPdg=nullptr,
                    *DepX=nullptr,
                    *DepY=nullptr,
                    *DepZ=nullptr,
                    *DepT=nullptr,
                    *DepE=nullptr;

    //PFParticle
    int NPfp;
    vector<int> *PfpTrackID=nullptr,
                *PfpShowerID=nullptr,
                *PfpID=nullptr,
                *PfpPdg=nullptr;

    //Track
    int NTrk;
    vector<int> *TrkID=nullptr,
                *TrkNPt=nullptr;
    vector<double>  *TrkLength=nullptr;
    vector<vector<Vec3D>>   *TrkPt=nullptr,
                            *TrkDir=nullptr;

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
    vector<vector<Vec3D>>   *Spt=nullptr;

    Reco(const char* filename) : file{new TFile(filename)} {
        reco = file->Get<TTree>("YAD/Reco");

        reco->SetBranchAddress("fTrkLength",  &TrkLength);

    }
    ~Reco() { file-> Close(); }
    
    int GetEntries() { return reco->GetEntries(); }
    void GetEntry(int i) { reco->GetEntry(i); }
    
};



}