/*
 *  classes to store and access variables from the output ROOT file of OmegaDumper_module
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

namespace omega {

class Truth {
private:

    TFile* file;
    TTree* tree;

    class MCParticles {
    public:
        size_t N;
        size_t NDau;
        bool isOrphelin;
        int Pdg;
        size_t NPt;
        vector<double*> X,
                        Y,
                        Z,
                        // T,
                        // Px,
                        // Py,
                        // Pz,
                        P,
                        E;
        
        void reset() {
            N=0;
            NDau=0;
            isOrphelin=true;
            Pdg=0;
            X.clear();
            Y.clear();
            Z.clear();
            // T.clear();
            // Px.clear();
            // Py.clear();
            // Pz.clear();
            P.clear();
            E.clear();
        }
    };

    class Deposits {
    public:
        size_t N;
        // vector<int*> Pdg;
        vector<double*> X,
                        Y,
                        Z,
                        // T,
                        E;
    
        void reset() {
            N=0;
            // Pdg.clear();
            X.clear();
            Y.clear();
            Z.clear();
            // T.clear();
            E.clear();
        }    
    };

public:
    //Events
    size_t Event, Run, SubRun;

private:
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
                            // *PrtT=nullptr,
                            // *PrtPx=nullptr,
                            // *PrtPy=nullptr,
                            // *PrtPz=nullptr,
                            *PrtP=nullptr,
                            *PrtE=nullptr;

    //SimEnergy*Deposit
    vector<size_t>  *PrtNDep=nullptr; //[# particle]
    vector<vector<double>>  *DepPdg=nullptr, //[# particle][# deposit]
                            *DepX=nullptr,
                            *DepY=nullptr,
                            *DepZ=nullptr,
                            // *DepT=nullptr,
                            *DepE=nullptr;

public:

    size_t N;
    MCParticles Prt;
    MCParticles Mom;
    MCParticles Dau;
    Deposits Dep;

    Truth(const char* filename) {
        file = new TFile(filename);
        tree = file->Get<TTree>("OmegaDumper/Truth");
        N=tree->GetEntries(); 

        //Events
        tree->SetBranchAddress("fEvent",    &Event);
        tree->SetBranchAddress("fRun",      &Run);
        tree->SetBranchAddress("fSubrun",   &SubRun);

        //MCTruth

        //MCParticlR.Spt.X[i_spt]e
        tree->SetBranchAddress("fNPrt",     &NPrt);
        tree->SetBranchAddress("fPrtPdg",   &PrtPdg); 
        tree->SetBranchAddress("fPrtMomID", &PrtMomID); 
        tree->SetBranchAddress("fPrtNDau",  &PrtNDau); 
        tree->SetBranchAddress("fPrtDauID", &PrtDauID); 
        tree->SetBranchAddress("fPrtNPt",   &PrtNPt); 
        tree->SetBranchAddress("fPrtX",     &PrtX); 
        tree->SetBranchAddress("fPrtY",     &PrtY); 
        tree->SetBranchAddress("fPrtZ",     &PrtZ); 
        // tree->SetBranchAddress("fPrtT",     &PrtT); 
        // tree->SetBranchAddress("fPrtPx",    &PrtPx); 
        // tree->SetBranchAddress("fPrtPy",    &PrtPy); 
        // tree->SetBranchAddress("fPrtPz",    &PrtPz); 
        tree->SetBranchAddress("fPrtP",     &PrtP); 
        tree->SetBranchAddress("fPrtE",     &PrtE); 

        //SimEnergyDeposit
        tree->SetBranchAddress("fPrtNDep",  &PrtNDep);
        tree->SetBranchAddress("fDepPdg",   &DepPdg); 
        tree->SetBranchAddress("fDepX",     &DepX); 
        tree->SetBranchAddress("fDepY",     &DepY); 
        tree->SetBranchAddress("fDepZ",     &DepZ); 
        // tree->SetBranchAddress("fDepT",     &DepT); 
        tree->SetBranchAddress("fDepE",     &DepE);
    }
    ~Truth() { file-> Close(); }
    
    void GetEvt(size_t i_evt) { 
        tree->GetEntry(i_evt);
        Prt.reset();
        Prt.N = NPrt;
    }

    void GetPrt(size_t i_prt) {
        Prt.reset();
        Prt.N = NPrt;
        Prt.NDau = PrtNDau->at(i_prt);
        Prt.isOrphelin = PrtMomID->at(i_prt) == -1;
        Prt.Pdg = PrtPdg->at(i_prt);
        Prt.NPt = PrtNPt->at(i_prt);
        for (size_t i_ppt=0; i_ppt < PrtNPt->at(i_prt); i_ppt++) {

            // Prt.index
            Prt.X.push_back(&(*PrtX)[i_prt][i_ppt]);
            Prt.Y.push_back(&(*PrtY)[i_prt][i_ppt]);
            Prt.Z.push_back(&(*PrtZ)[i_prt][i_ppt]);
            // Prt.T.push_back(&(*PrtT)[i_prt][i_ppt]);
            // Prt.Px.push_back(&(*PrtPx)[i_prt][i_ppt]);
            // Prt.Py.push_back(&(*PrtPy)[i_prt][i_ppt]);
            // Prt.Pz.push_back(&(*PrtPz)[i_prt][i_ppt]);
            Prt.P.push_back(&(*PrtP)[i_prt][i_ppt]);
            Prt.E.push_back(&(*PrtE)[i_prt][i_ppt]);
        }
    }
    void GetPrtMom(size_t i_prt) {
        Mom.reset();
        size_t i_mom = PrtMomID->at(i_prt);
        Mom.NDau = PrtNDau->at(i_mom);
        Mom.isOrphelin = PrtMomID->at(i_mom) == -1;
        Mom.NPt = PrtNPt->at(i_mom);
        for (size_t i_ppt=0; i_ppt < PrtNPt->at(i_mom); i_ppt++) {

            Mom.X.push_back(&(*PrtX)[i_mom][i_ppt]);
            Mom.Y.push_back(&(*PrtY)[i_mom][i_ppt]);
            Mom.Z.push_back(&(*PrtZ)[i_mom][i_ppt]);
            // Mom.T.push_back(&(*PrtT)[i_mom][i_ppt]);
            // Mom.Px.push_back(&(*PrtPx)[i_mom][i_ppt]);
            // Mom.Py.push_back(&(*PrtPy)[i_mom][i_ppt]);
            // Mom.Pz.push_back(&(*PrtPz)[i_mom][i_ppt]);
            Mom.P.push_back(&(*PrtP)[i_mom][i_ppt]);
            Mom.E.push_back(&(*PrtE)[i_mom][i_ppt]);
        }
    }
    void GetPrtDau(size_t i_prt, size_t j_dau) {
        Dau.reset();
        size_t i_dau = (*PrtDauID)[i_prt][j_dau];
        Dau.NDau = PrtNDau->at(i_dau);
        Dau.isOrphelin = PrtMomID->at(i_dau) == -1;
        Dau.NPt = PrtNPt->at(i_dau);
        for (size_t i_ppt=0; i_ppt < PrtNPt->at(i_dau); i_ppt++) {

            Dau.X.push_back(&(*PrtX)[i_dau][i_ppt]);
            Dau.Y.push_back(&(*PrtY)[i_dau][i_ppt]);
            Dau.Z.push_back(&(*PrtZ)[i_dau][i_ppt]);
            // Dau.T.push_back(&(*PrtT)[i_dau][i_ppt]);
            // Dau.Px.push_back(&(*PrtPx)[i_dau][i_ppt]);
            // Dau.Py.push_back(&(*PrtPy)[i_dau][i_ppt]);
            // Dau.Pz.push_back(&(*PrtPz)[i_dau][i_ppt]);
            Dau.P.push_back(&(*PrtP)[i_dau][i_ppt]);
            Dau.E.push_back(&(*PrtE)[i_dau][i_ppt]);
        }
    }
    void GetPrtDep(size_t i_prt) {
        Dep.reset();
        Dep.N = PrtNDep->at(i_prt);
        for (size_t i_dep=0; i_dep < PrtNDep->at(i_prt); i_dep++) {

            Dep.X.push_back(&(*DepX)[i_prt][i_dep]);
            Dep.Y.push_back(&(*DepY)[i_prt][i_dep]);
            Dep.Z.push_back(&(*DepZ)[i_prt][i_dep]);
            // Dep.T.push_back(&(*DepT)[i_prt][i_dep]);
            Dep.E.push_back(&(*DepE)[i_prt][i_dep]);
        }
    }
}; //end of omega::Truth


class Reco {
private:

    TFile* file;
    TTree* tree;

    class PFParticles {
    public:
        size_t N;
        // vector<size_t>  index;
        vector<bool>    isTrk;
                        // isShw;
        // vector<size_t*> NClu,
        vector<size_t*> NSpt,
                        NHit;

        void reset() {
            N=0;
            // index.clear();
            isTrk.clear();
            // isShw.clear();
            // NClu.clear();
            NSpt.clear();
            NHit.clear();
        }
    };
    class Tracks {
    public:
        size_t  N,
                NHit,
                NPt;
        double Length;
        vector<double*> X,
                        Y,
                        Z;
                        // DirX,
                        // DirY,
                        // DirZ;

        void reset() {
            N=0;
            NHit=0;
            NPt=0;
            Length=0;
            X.clear();
            Y.clear();
            Z.clear();
            // DirX.clear();
            // DirY.clear();
            // DirZ.clear();
        }
    };
    class Calorimetry {
    public:
        size_t Plane,
               NPt;
        float Range;
            //   KineticE;
        vector<double*> dEdx,
                        dQdx,
                        ResRange;

        void reset() {
            Plane=0;
            NPt=0;
            Range=0;
            // KineticE=0;
            dEdx.clear();
            dQdx.clear();
            ResRange.clear();
        }
    };
    // class Clusters {
    // public:
    //     size_t N;
    //     vector<size_t>  index;
    //     vector<size_t*> NHit,
    //                     Plane;
    //     vector<float*>  SumADC,
    //                     Integral,
    //                     Width;
        
    //     void reset() {
    //         N=0;
    //         index.clear();
    //         NHit.clear();
    //         Plane.clear();
    //         SumADC.clear();
    //         Integral.clear();
    //         Width.clear();
    //     }
    // };
    class SpacePoints {
    public: 
        size_t N;
        vector<size_t>  index;
        vector<size_t*> NHit;
        vector<double*> X,
                        Y,
                        Z;

        void reset() {
            N=0;
            index.clear();
            NHit.clear();
            X.clear();
            Y.clear();
            Z.clear();
        }
    };
    class Hits {
    public: 
        size_t N;
        vector<size_t>  index;
        vector<size_t*> NTrk,
                        NSpt,
                        // NClu,
                        Plane;
        vector<float*>  SumADC,
                        Integral;

        void reset() {
            N=0;
            index.clear();
            NTrk.clear();
            NSpt.clear();
            // NClu.clear();
            Plane.clear();
            SumADC.clear();
            Integral.clear();
        }
    };

public:
    //Events
    size_t Event, Run, SubRun;

private:
    //PFParticle
    size_t NPfp;
    vector<int> *PfpTrkID=nullptr; //[# pfparticle]
                // *PfpShwID=nullptr;
    // vector<size_t>  *PfpNClu=nullptr,
    vector<size_t>  *PfpNSpt=nullptr,
                    *PfpNHit=nullptr;
    // vector<vector<double>>  *PfpCluKey=nullptr, //should be ?
    vector<vector<double>>  *PfpSptKey=nullptr;

    //Track
    size_t NTrk;
    vector<size_t>  *TrkKey=nullptr, //[# track]
                    *TrkNPt=nullptr, 
                    *TrkNHit=nullptr;
    vector<double>  *TrkLength=nullptr; //[# track]
    vector<vector<double>>  *TrkPtX=nullptr, //[# track][# point]
                            *TrkPtY=nullptr,
                            *TrkPtZ=nullptr,
                            // *TrkDirX=nullptr,
                            // *TrkDirY=nullptr,
                            // *TrkDirZ=nullptr,
                            *TrkHitKey=nullptr;

    //Calorimetry
    vector<size_t>  *TrkCalPlane=nullptr, //[# track]
                    *TrkCalNPt=nullptr;
    vector<float>   *TrkCalRange=nullptr;
                    // *TrkCalKineticE=nullptr;
    vector<vector<double>>  *TrkCaldEdx=nullptr, //[# track][# calo point]
                            *TrkCaldQdx=nullptr,
                            *TrkCalResRange=nullptr;

    //Shower
    // size_t NShw;
    // vector<int> *ShwID=nullptr; //[# shower]

    //Cluster
    // vector<size_t>  *CluKey=nullptr,
    //                 *CluNHit=nullptr,
    //                 *CluPlane=nullptr; //should be unsigned int
    // vector<float>   *CluSumADC=nullptr,
    //                 *CluIntegral=nullptr,
    //                 *CluWidth=nullptr;
    // vector<vector<double>>  *CluHitKey=nullptr; //should be ?

    //SpacePoint
    vector<size_t>  *SptKey=nullptr,
                    *SptNHit=nullptr;
    vector<double>  *SptX=nullptr, //should be Double32_t
                    *SptY=nullptr,
                    *SptZ=nullptr;
    vector<vector<double>>  *SptHitKey=nullptr; //should be ?

    //Hit
    size_t NHit;
    vector<size_t>  *HitKey=nullptr,
                    *HitNTrk=nullptr,
                    *HitNSpt=nullptr,
                    // *HitNClu=nullptr,  
                    *HitPlane=nullptr; //should be unsigned int
    vector<float>   *HitSumADC=nullptr, //should be float
                    *HitIntegral=nullptr; //should be float
    vector<vector<double>>  *HitTrkKey=nullptr, //should be ?
                            *HitSptKey=nullptr;
                            // *HitCluKey=nullptr;

public:

    size_t N;
    PFParticles Pfp;
    Tracks Trk;
    Calorimetry Cal;
    // Clusters Clu;
    SpacePoints Spt;
    Hits Hit;

    Reco(const char* filename) {
        file = new TFile(filename);
        tree = file->Get<TTree>("OmegaDumper/Reco");
        N=tree->GetEntries();

        //Events
        tree->SetBranchAddress("fEvent",        &Event);
        tree->SetBranchAddress("fRun",          &Run);
        tree->SetBranchAddress("fSubrun",       &SubRun);

        //PFParticle
        tree->SetBranchAddress("fNPfp",         &NPfp);
        tree->SetBranchAddress("fPfpTrkID",     &PfpTrkID);
        // tree->SetBranchAddress("fPfpShwID",     &PfpShwID);
        // tree->SetBranchAddress("fPfpNClu",      &PfpNClu);
        tree->SetBranchAddress("fPfpNSpt",      &PfpNSpt);
        tree->SetBranchAddress("fPfpNHit",      &PfpNHit);
        // tree->SetBranchAddress("fPfpCluKey",    &PfpCluKey);
        tree->SetBranchAddress("fPfpSptKey",    &PfpSptKey);

        //Track
        tree->SetBranchAddress("fNTrk",         &NTrk);
        tree->SetBranchAddress("fTrkKey",       &TrkKey);
        tree->SetBranchAddress("fTrkNHit",      &TrkNHit);
        tree->SetBranchAddress("fTrkNPt",       &TrkNPt);
        tree->SetBranchAddress("fTrkLength",    &TrkLength);
        tree->SetBranchAddress("fTrkPtX",       &TrkPtX);
        tree->SetBranchAddress("fTrkPtY",       &TrkPtY);
        tree->SetBranchAddress("fTrkPtZ",       &TrkPtZ);
        // tree->SetBranchAddress("fTrkDirX",      &TrkDirX);
        // tree->SetBranchAddress("fTrkDirY",      &TrkDirY);
        // tree->SetBranchAddress("fTrkDirZ",      &TrkDirZ);
        tree->SetBranchAddress("fTrkHitKey",    &TrkHitKey);

        //Calorimetry
        tree->SetBranchAddress("fTrkCalPlane",   &TrkCalPlane);
        tree->SetBranchAddress("fTrkCalRange",   &TrkCalRange);
        // tree->SetBranchAddress("fTrkCalKineticE",&TrkCalKineticE);
        tree->SetBranchAddress("fTrkCalNPt",     &TrkCalNPt);
        tree->SetBranchAddress("fTrkCaldEdx",    &TrkCaldEdx);
        tree->SetBranchAddress("fTrkCaldQdx",    &TrkCaldQdx);
        tree->SetBranchAddress("fTrkCalResRange",&TrkCalResRange);

        //Shower
        // tree->SetBranchAddress("fNShw",         &NShw);
        // tree->SetBranchAddress("fShwID",        &ShwID);

        //Cluster
        // tree->SetBranchAddress("fCluKey",       &CluKey);
        // tree->SetBranchAddress("fCluNHit",      &CluNHit);
        // tree->SetBranchAddress("fCluPlane",     &CluPlane);
        // tree->SetBranchAddress("fCluIntegral",  &CluIntegral);
        // tree->SetBranchAddress("fCluSumADC",    &CluSumADC);
        // tree->SetBranchAddress("fCluWidth",     &CluWidth);
        // tree->SetBranchAddress("fCluHitKey",    &CluHitKey);

        //SpacePoint
        tree->SetBranchAddress("fSptKey",       &SptKey);
        tree->SetBranchAddress("fSptNHit",      &SptNHit);
        tree->SetBranchAddress("fSptX",         &SptX);
        tree->SetBranchAddress("fSptY",         &SptY);
        tree->SetBranchAddress("fSptZ",         &SptZ);
        tree->SetBranchAddress("fSptHitKey",    &SptHitKey);

        //Hit
        tree->SetBranchAddress("fNHit",         &NHit);
        tree->SetBranchAddress("fHitKey",       &HitKey);
        tree->SetBranchAddress("fHitNTrk",      &HitNTrk);
        tree->SetBranchAddress("fHitNSpt",      &HitNSpt);
        // tree->SetBranchAddress("fHitNClu",      &HitNClu);
        tree->SetBranchAddress("fHitPlane",     &HitPlane);
        tree->SetBranchAddress("fHitSumADC",    &HitSumADC);
        tree->SetBranchAddress("fHitIntegral",  &HitIntegral);
        tree->SetBranchAddress("fHitTrkKey",    &HitTrkKey);
        tree->SetBranchAddress("fHitSptKey",    &HitSptKey);
        // tree->SetBranchAddress("fHitSptKey",    &HitCluKey);

    }
    ~Reco() { file-> Close(); }
    
    void GetEvtPfp(size_t i_evt) { 
        tree->GetEntry(i_evt); 
        Pfp.reset();
        Pfp.N = NPfp;
        Trk.N = NTrk;
        for (size_t i_pfp=0; i_pfp < NPfp; i_pfp++) {

            // Pfp.index.push_back(i_pfp);
            Pfp.isTrk.push_back(PfpTrkID->at(i_pfp) >= 0);
            // Pfp.isShw.push_back(PfpShwID->at(i_pfp) >= 0);
            // Pfp.NClu.push_back(&PfpNClu->at(i_pfp));
            Pfp.NSpt.push_back(&PfpNSpt->at(i_pfp));
            Pfp.NHit.push_back(&PfpNHit->at(i_pfp));
        }
    }
    void GetEvtHit(size_t i_evt) { 
        tree->GetEntry(i_evt); 
        Hit.reset();
        Hit.N = NHit;
        for (size_t i_hit=0; i_hit < NHit; i_hit++) {

            Hit.index.push_back(i_hit);
            Hit.NTrk.push_back(&HitNTrk->at(i_hit));
            Hit.NSpt.push_back(&HitNSpt->at(i_hit));
            // Hit.NClu.push_back(&HitNClu->at(i_hit));
            Hit.Plane.push_back(&HitPlane->at(i_hit));
            Hit.SumADC.push_back(&HitSumADC->at(i_hit));
            Hit.Integral.push_back(&HitIntegral->at(i_hit));
        }
    }
    void GetPfpTrk(size_t i_pfp) {
        Trk.reset();
        size_t i_trk = PfpTrkID->at(i_pfp);
        Trk.Length = TrkLength->at(i_trk);
        Trk.NPt = TrkNPt->at(i_trk);
        for (size_t i_tpt=0; i_tpt < TrkNPt->at(i_trk); i_tpt++) {
            Trk.X.push_back(&(*TrkPtX)[i_trk][i_tpt]);
            Trk.Y.push_back(&(*TrkPtY)[i_trk][i_tpt]);
            Trk.Z.push_back(&(*TrkPtZ)[i_trk][i_tpt]);
            // Trk.DirX.push_back(&(*TrkDirX)[i_trk][i_tpt]);
            // Trk.DirY.push_back(&(*TrkDirY)[i_trk][i_tpt]);
            // Trk.DirZ.push_back(&(*TrkDirZ)[i_trk][i_tpt]);
        }
        Cal.reset();
        Cal.Plane = TrkCalPlane->at(i_trk);
        Cal.NPt = TrkCalNPt->at(i_trk);
        Cal.Range = TrkCalRange->at(i_trk);
        // Cal.KineticE = TrkCalKineticE->at(i_trk);
        for (size_t i_cal=0; i_cal < TrkCalNPt->at(i_trk); i_cal++) {
            Cal.dEdx.push_back(&(*TrkCaldEdx)[i_trk][i_cal]);
            Cal.dQdx.push_back(&(*TrkCaldQdx)[i_trk][i_cal]);
            Cal.ResRange.push_back(&(*TrkCalResRange)[i_trk][i_cal]);
        }
    }
    // void GetPfpClu(size_t i_pfp) {
    //     Clu.reset();
    //     Clu.N = PfpNClu->at(i_pfp);
    //     for (size_t key : PfpCluKey->at(i_pfp)) {

    //         size_t i_clu = std::distance(CluKey->begin(), std::find(CluKey->begin(), CluKey->end(), key));

    //         Clu.index.push_back(i_clu);
    //         Clu.NHit.push_back(&CluNHit->at(i_clu));
    //         Clu.Plane.push_back(&CluPlane->at(i_clu));
    //         Clu.SumADC.push_back(&CluSumADC->at(i_clu));
    //         Clu.Integral.push_back(&CluIntegral->at(i_clu));
    //         Clu.Width.push_back(&CluWidth->at(i_clu));
    //     }
    // }
    void GetPfpSpt(size_t i_pfp) {
        Spt.reset();
        Spt.N = PfpNSpt->at(i_pfp);
        for (size_t key : PfpSptKey->at(i_pfp)) {

            size_t i_spt = std::distance(SptKey->begin(), std::find(SptKey->begin(), SptKey->end(), key));

            Spt.index.push_back(i_spt);
            Spt.NHit.push_back(&SptNHit->at(i_spt));
            Spt.X.push_back(&SptX->at(i_spt));
            Spt.Y.push_back(&SptY->at(i_spt));
            Spt.Z.push_back(&SptZ->at(i_spt));
        }
    }
    void GetTrkHit(size_t i_trk) {
        Hit.reset();
        Hit.N = TrkNHit->at(i_trk);
        for (size_t key : TrkHitKey->at(i_trk)) {

            size_t i_hit = std::distance(HitKey->begin(), std::find(HitKey->begin(), HitKey->end(), key));

            Hit.index.push_back(i_hit);
            Hit.NTrk.push_back(&HitNTrk->at(i_hit));
            Hit.NSpt.push_back(&HitNSpt->at(i_hit));
            // Hit.NClu.push_back(&HitNClu->at(i_hit));
            Hit.Plane.push_back(&HitPlane->at(i_hit));
            Hit.SumADC.push_back(&HitSumADC->at(i_hit));
            Hit.Integral.push_back(&HitIntegral->at(i_hit));
        }
    }
    // void GetCluHit(size_t i_clu) {
    //     Hit.reset();
    //     Hit.N = CluNHit->at(i_clu);
    //     for (size_t key : CluHitKey->at(i_clu)) {

    //         size_t i_hit = std::distance(HitKey->begin(), std::find(HitKey->begin(), HitKey->end(), key));

    //         Hit.index.push_back(i_hit);
    //         Hit.NTrk.push_back(&HitNTrk->at(i_hit));
    //         Hit.NSpt.push_back(&HitNSpt->at(i_hit));
    //         Hit.NClu.push_back(&HitNClu->at(i_hit));
    //         Hit.Plane.push_back(&HitPlane->at(i_hit));
    //         Hit.SumADC.push_back(&HitSumADC->at(i_hit));
    //         Hit.Integral.push_back(&HitIntegral->at(i_hit));
    //     }
    // }
    void GetSptHit(size_t i_spt) {
        Hit.reset();
        Hit.N = SptNHit->at(i_spt);
        for (size_t key : SptHitKey->at(i_spt)) {

            size_t i_hit = std::distance(HitKey->begin(), std::find(HitKey->begin(), HitKey->end(), key));

            Hit.index.push_back(i_hit);
            Hit.NSpt.push_back(&HitNSpt->at(i_hit));
            // Hit.NClu.push_back(&HitNClu->at(i_hit));
            Hit.Plane.push_back(&HitPlane->at(i_hit));
            Hit.SumADC.push_back(&HitSumADC->at(i_hit));
            Hit.Integral.push_back(&HitIntegral->at(i_hit));
        }
    }
    void GetHitTrk(size_t i_hit) {
        Trk.reset();
        size_t key = (*HitTrkKey)[i_hit][0];
        size_t i_trk = std::distance(TrkKey->begin(), std::find(TrkKey->begin(), TrkKey->end(), key));
        Trk.Length = TrkLength->at(i_trk);
        Trk.NPt = TrkNPt->at(i_trk);
        for (size_t i_tpt=0; i_tpt < TrkNPt->at(i_trk); i_tpt++) {
            Trk.X.push_back(&(*TrkPtX)[i_trk][i_tpt]);
            Trk.Y.push_back(&(*TrkPtY)[i_trk][i_tpt]);
            Trk.Z.push_back(&(*TrkPtZ)[i_trk][i_tpt]);
            // Trk.DirX.push_back(&(*TrkDirX)[i_trk][i_tpt]);
            // Trk.DirY.push_back(&(*TrkDirY)[i_trk][i_tpt]);
            // Trk.DirZ.push_back(&(*TrkDirZ)[i_trk][i_tpt]);
        }
        Cal.reset();
        Cal.Plane = TrkCalPlane->at(i_trk);
        Cal.NPt = TrkCalNPt->at(i_trk);
        Cal.Range = TrkCalRange->at(i_trk);
        // Cal.KineticE = TrkCalKineticE->at(i_trk);
        for (size_t i_cal=0; i_cal < TrkCalNPt->at(i_trk); i_cal++) {
            Cal.dEdx.push_back(&(*TrkCaldEdx)[i_trk][i_cal]);
            Cal.dQdx.push_back(&(*TrkCaldQdx)[i_trk][i_cal]);
            Cal.ResRange.push_back(&(*TrkCalResRange)[i_trk][i_cal]);
        }
    }
    // void GetHitClu(size_t i_hit) {
    //     Clu.reset();
    //     Clu.N = HitNClu->at(i_hit);
    //     for (size_t key : HitCluKey->at(i_hit)) {

    //         size_t i_clu = std::distance(CluKey->begin(), std::find(CluKey->begin(), CluKey->end(), key));

    //         Clu.index.push_back(i_clu);
    //         Clu.NHit.push_back(&CluNHit->at(i_clu));
    //         Clu.Plane.push_back(&CluPlane->at(i_clu));
    //         Clu.SumADC.push_back(&CluSumADC->at(i_clu));
    //         Clu.Integral.push_back(&CluIntegral->at(i_clu));
    //         Clu.Width.push_back(&CluWidth->at(i_clu));
    //     }
    // }
    void GetHitSpt(size_t i_hit) {
        Spt.reset();
        Spt.N = HitNSpt->at(i_hit);
        for (size_t key : HitSptKey->at(i_hit)) {

            size_t i_spt = std::distance(SptKey->begin(), std::find(SptKey->begin(), SptKey->end(), key));

            Spt.index.push_back(i_spt);
            Spt.NHit.push_back(&SptNHit->at(i_spt));
            Spt.X.push_back(&SptX->at(i_spt));
            Spt.Y.push_back(&SptY->at(i_spt));
            Spt.Z.push_back(&SptZ->at(i_spt));
        }
    }

}; //end of omega::Reco

vector<string> ReadFileList(size_t n_file, string listname) {
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
} //end of yad::ReadFileList

typedef struct {
    size_t  n;
    double  min,
            max;
} Binning;

typedef struct {
    double  Xmin, //cm
            Xmax,
            Ymin,
            Ymax,
            Zmin,
            Zmax;
} Limits;
const Limits cryostat = {-350,350,-350,350,0,300};
const Limits cryostat2 = {-345,345,-345,345,5,295};
const Limits fiducial = {-320,350,-320,320,20,280};

bool IsInside (vector<double*> Xs,vector<double*> Ys,vector<double*> Zs, Limits lim = cryostat) {
    bool is_in=true;
    vector<size_t> indexes = {0,Zs.size()-1};
    for (size_t i : indexes) {
        is_in = lim.Xmin < *Xs[i] && *Xs[i] < lim.Xmax &&
                lim.Ymin < *Ys[i] && *Ys[i] < lim.Ymax &&
                lim.Zmin < *Zs[i] && *Zs[i] < lim.Zmax;
        if (!is_in) break;
    }
    return is_in;
} //end of omega::IsInside

double Distance(double *X1, double *Y1, double *Z1, double *X2, double *Y2, double *Z2) {
    return TMath::Sqrt(
        TMath::Power(*X1-*X2,2)
        +TMath::Power(*Y1-*Y2,2)
        +TMath::Power(*Z1-*Z2,2)
    );
} //end of omega::Distance

} //end of namespace omega