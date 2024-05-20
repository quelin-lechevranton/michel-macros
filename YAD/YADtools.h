#include <vector>
#include <string>
#include <iostream>

using namespace std;

#include <TFile.h>
#include <TTree.h>
// #include <Vector3D.h>
#include "Math/GenVector/LorentzVector.h"

using Vec3D = ROOT::Math::XYZVector;
using Vec4D = ROOT::Math::Lorentz<ROOT::Math::PxPyPzEVector>;

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
    vector<int> PrtPdg;
    vector<double>  PrtStX,
                    PrtStY,
                    PrtStZ,
                    PrtStT,
                    PrtStPx, 
                    PrtStPy, 
                    PrtStPz, 
                    PrtStP, 
                    PrtStE;

    //SimEnergyDeposit
    int NDep;
    vector<double>  DepPdg,
                    DepX,
                    DepY,
                    DepZ,
                    DepT,
                    DepE;

    //PFParticle
    int NPfp;
    vector<int> PfpTrackID,
                PfpShowerID,
                PfpID,
                PfpPdg;

    //Track
    int NTrk;
    vector<int> TrkID;
    vector<double>  TrkLength;
    vector<int> TrkNPt;
    vector<vector<Vec3D>>   TrkPt,
                            TrkDir;

    //Calorimetry
    vector<double>  TrkCalRange;

    //Shower
    int NShw;
    vector<int> ShwID;

    //Cluster
    vector<int> PfpNClu;
    vector<vector<int>>     CluNHit,
                            CluPlane;
    vector<vector<double>>  CluIntFit,
                            CluSumADC,
                            CluWidth;

    //SpacePoint
    vector<int> PfpNSpt;
    vector<vector<Vec3D>>  Spt;

    Reco(const char* filename) : file{new TFile(filename)} {
        reco = file->Get<TTree>("YAD/Reco");

        reco->SetBranchAddress("fTrkLength",  &TrkLength);

    }
    ~Reco() { file-> Close; }
    
    int GetEntries() { return reco->GetEntries(); }
    void GetEntry(int i) { reco->GetEntry(i); }
    
};



}