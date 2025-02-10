/*
 * Based on Dumper module by Laura Perez Molina
 */


// Art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
// #include "larsim/MCCheater/BackTrackerService.h"
// #include "larsim/MCCheater/PhotonBackTrackerService.h"
// #include "larsim/MCCheater/ParticleInventoryService.h"
// #include "larsim/MCCheater/BackTracker.h"
// #include "larsim/MCCheater/BackTrackerService.h"
#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"
#include "protoduneana/Utilities/ProtoDUNEPFParticleUtils.h"
#include "protoduneana/Utilities/ProtoDUNETrackUtils.h"

// DUNE includes
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"

// ROOT includes
#include <TTree.h>
// #include "Math/GenVector/LorentzVector.h"
// #include "Math/GenVector/PositionVector3D.h"
// #include <TVector3.h>

// std includes
#include <vector>
#include <iterator> 
#include <string>

using namespace std;

namespace ana { class YAD; }

class ana::YAD : public art::EDAnalyzer {
public:

    explicit YAD(fhicl::ParameterSet const& fcl);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    YAD(YAD const &) = delete;
    YAD(YAD &&) = delete;
    YAD & operator = (YAD const &) = delete;
    YAD & operator = (YAD &&) = delete;
    
    // Required functions.
    void analyze(art::Event const & evt) override; 
    void beginJob() override;
    // void endJob()   override;
    void reset();

private:

    protoana::ProtoDUNEPFParticleUtils util_pfp;
    protoana::ProtoDUNETrackUtils util_trk;

    vector<string> fTrees;
    vector<vector<string>> fProducts;

    TTree *fTruth;
    TTree *fReco;

    art::InputTag   tag_tru,
                    tag_prt,
                    tag_dep,
                    tag_pfp,
                    tag_clu,
                    tag_trk,
                    tag_shw,
                    tag_spt,
                    tag_cal;


    //Events
    size_t fEvent, fRun, fSubRun;

    //MCTruth
    size_t fNPrt;

    //MCParticle
    vector<int> fPrtPdg;
    vector<int> fPrtMomID;
    vector<size_t>  fPrtNPt,
                    fPrtNDau;
    vector<vector<double>>  fPrtDauID, //should be int
                            fPrtX,
                            fPrtY,
                            fPrtZ,
                            fPrtT,
                            fPrtPx, 
                            fPrtPy, 
                            fPrtPz, 
                            fPrtP, 
                            fPrtE;

    //SimEnergyDeposit
    vector<size_t>  fPrtNDep;
    vector<vector<double>>  fDepPdg, //should be int
                            fDepX,
                            fDepY,
                            fDepZ,
                            fDepT,
                            fDepE;

    //PFParticle
    size_t fNPfp;
    vector<int> fPfpTrkID,
                fPfpShwID,
                fPfpPdg; //useless ?

    //Track
    size_t fNTrk;
    vector<double>  fTrkLength;
    vector<size_t> fTrkNPt;
    vector<vector<double>>  fTrkPtX, //should be Double_t
                            fTrkPtY,
                            fTrkPtZ,
                            fTrkDirX,
                            fTrkDirY,
                            fTrkDirZ;

    //Calorimetry
    vector<size_t>  fTrkCalPlane,   //should be unsigned int
                                    //TrackUtils::GetCalorimetry only gives calorimetry on the collection: Plane = geo::kW = 2
                    fTrkCalNPt;
    vector<float>   fTrkCalRange,
                    fTrkCalKineticE;
    vector<vector<double>>  fTrkCaldEdx, //should be float
                            fTrkCaldQdx, //should be float
                            fTrkCalResRange; //should be float

    //Shower
    size_t fNShw;
    vector<int> fShwID; //useless?

    //Cluster
    vector<size_t> fPfpNClu;
    vector<vector<double>>  fCluNHit, //should be size_t
                            fCluPlane; //should be unsigned int
    vector<vector<double>>  fCluIntFit, //should be float
                            fCluSumADC, //should be float
                            fCluWidth; //should be float

    //SpacePoint
    vector<size_t> fPfpNSpt;
    vector<vector<double>>  fSptX, //should be Double32_t
                            fSptY,
                            fSptZ;

};

ana::YAD::YAD(fhicl::ParameterSet const & fcl) :
    EDAnalyzer{fcl} 
{
    art::ServiceHandle<art::TFileService> tfs;
    fTrees =    fcl.get<vector<string>>("Trees");
    fProducts = fcl.get<vector<vector<string>>>("Products");

    for (vector<string> prod : fProducts) {

        string  label=prod[0],
                instance=prod[1],
                object=prod[2],
                process=prod[3];

        if (object == "simb::MCTruth")              tag_tru= art::InputTag(label,instance);
        else if (object == "simb::MCParticle")      tag_prt= art::InputTag(label,instance);
        else if (object == "sim::SimEnergyDeposit") tag_dep= art::InputTag(label,instance);
        else if (object == "recob::PFParticle")     tag_pfp= art::InputTag(label,instance);
        else if (object == "recob::Cluster")        tag_clu= art::InputTag(label,instance);
        else if (object == "recob::Track")          tag_trk= art::InputTag(label,instance);
        else if (object == "anab::Calorimetry")     tag_cal= art::InputTag(label,instance);
        else if (object == "recob::Shower")         tag_shw= art::InputTag(label,instance);
        else if (object == "recob::SpacePoint")     tag_spt= art::InputTag(label,instance);
    } //end prod loop

    if ( find(fTrees.begin(), fTrees.end(), "Truth") != fTrees.end() ) {

        fTruth = tfs->make<TTree>("Truth","Truth");

        //Events
        fTruth->Branch("fEvent",        &fEvent);
        fTruth->Branch("fRun",          &fRun);
        fTruth->Branch("fSubrun",       &fSubRun);

        //MCTruth
        // if (!tag_tru.empty()) {

        // }

        //MCParticle
        if (!tag_prt.empty()) {
        fTruth->Branch("fNPrt",         &fNPrt);
        fTruth->Branch("fPrtPdg",       &fPrtPdg); 
        fTruth->Branch("fPrtMomID",     &fPrtMomID); 
        fTruth->Branch("fPrtNDau",      &fPrtNDau); 
        fTruth->Branch("fPrtDauID",     &fPrtDauID); 
        fTruth->Branch("fPrtNPt",       &fPrtNPt); 
        fTruth->Branch("fPrtX",         &fPrtX); 
        fTruth->Branch("fPrtY",         &fPrtY); 
        fTruth->Branch("fPrtZ",         &fPrtZ); 
        fTruth->Branch("fPrtT",         &fPrtT); 
        fTruth->Branch("fPrtPx",        &fPrtPx); 
        fTruth->Branch("fPrtPy",        &fPrtPy); 
        fTruth->Branch("fPrtPz",        &fPrtPz); 
        fTruth->Branch("fPrtP",         &fPrtP); 
        fTruth->Branch("fPrtE",         &fPrtE); 
        }

        //SimEnergyDeposit
        if (!tag_dep.empty()) {
        fTruth->Branch("fPrtNDep",      &fPrtNDep);
        fTruth->Branch("fDepPdg",       &fDepPdg); 
        fTruth->Branch("fDepX",         &fDepX); 
        fTruth->Branch("fDepY",         &fDepY); 
        fTruth->Branch("fDepZ",         &fDepZ); 
        fTruth->Branch("fDepT",         &fDepT); 
        fTruth->Branch("fDepE",         &fDepE); 
        }

    } //end Truth tree

    if ( find(fTrees.begin(), fTrees.end(), "Reco") != fTrees.end() ) {

        fReco = tfs->make<TTree>("Reco","Reco");

        //Events
        fReco->Branch("fEvent",         &fEvent);
        fReco->Branch("fRun",           &fRun);
        fReco->Branch("fSubrun",        &fSubRun);

        //PFParticle
        if (!tag_pfp.empty()) {
        fReco->Branch("fNPfp",          &fNPfp);
        fReco->Branch("fPfpTrkID",      &fPfpTrkID);
        fReco->Branch("fPfpShwID",      &fPfpShwID);
        fReco->Branch("fPfpPdg",        &fPfpPdg);
        }

        //Track
        if (!tag_trk.empty()) {
        fReco->Branch("fNTrk",          &fNTrk);
        fReco->Branch("fTrkLength",     &fTrkLength);
        fReco->Branch("fTrkNPt",        &fTrkNPt);
        fReco->Branch("fTrkPtX",        &fTrkPtX);
        fReco->Branch("fTrkPtY",        &fTrkPtY);
        fReco->Branch("fTrkPtZ",        &fTrkPtZ);
        fReco->Branch("fTrkDirX",       &fTrkDirX);
        fReco->Branch("fTrkDirY",       &fTrkDirY);
        fReco->Branch("fTrkDirZ",       &fTrkDirZ);
        }

        //Calorimetry
        if (!tag_cal.empty()) {
        fReco->Branch("fTrkCalPlane",   &fTrkCalPlane);
        fReco->Branch("fTrkCalRange",   &fTrkCalRange);
        fReco->Branch("fTrkCalKineticE",&fTrkCalKineticE);
        fReco->Branch("fTrkCalNPt",     &fTrkCalNPt);
        fReco->Branch("fTrkCaldEdx",    &fTrkCaldEdx);
        fReco->Branch("fTrkCaldQdx",    &fTrkCaldQdx);
        fReco->Branch("fTrkCalResRange",&fTrkCalResRange);
        }

        //Shower
        if (!tag_shw.empty()) {
        fReco->Branch("fNShw",          &fNShw);
        fReco->Branch("fShwID",         &fShwID);
        }

        //Cluster
        if (!tag_clu.empty()) {
        fReco->Branch("fPfpNClu",       &fPfpNClu);
        fReco->Branch("fCluNHit",       &fCluNHit);
        fReco->Branch("fCluIntFit",     &fCluIntFit);
        fReco->Branch("fCluSumADC",     &fCluSumADC);
        fReco->Branch("fCluWidth",      &fCluWidth);
        fReco->Branch("fCluPlane",      &fCluPlane);
        }

        //SpacePoint
        if (!tag_spt.empty()) {
        fReco->Branch("fPfpNSpt",       &fPfpNSpt);
        fReco->Branch("fSptX",          &fSptX);
        fReco->Branch("fSptY",          &fSptY);
        fReco->Branch("fSptZ",          &fSptZ);
        }

    } //end Reco tree

} //end YAD()

void ana::YAD::beginJob() {} //end beginJob()

void ana::YAD::analyze(const art::Event & evt) {
    reset();
    fEvent  = evt.id().event(); 
    fRun    = evt.id().run();
    fSubRun = evt.id().subRun();

    //Truth Tree

    //simb::MCTruth
    // if (!tag_tru.empty()) {
    // auto vh_tru = evt.getValidHandle<vector<simb::MCTruth>>(tag_tru);

    // for (simb::MCTruth const & tru : *vh_tru) {

    // } //end vh_tru loop
    // } //end simb::MCTruth

    //simb::MCParticle
    if (!tag_prt.empty()) {
    auto vh_prt = evt.getValidHandle<vector<simb::MCParticle>>(tag_prt);
    
    fNPrt=vh_prt->size();

    for (simb::MCParticle const & prt : *vh_prt) {

        fPrtPdg.        push_back(prt.PdgCode());
        fPrtMomID.      push_back(prt.Mother()-1); //IDs start at 1, we want 0 for vector indexing
        fPrtNDau.       push_back(prt.NumberDaughters());

        vector<double> tpPrtDauID;
        for (int i_dau=0; i_dau < prt.NumberDaughters(); i_dau++) {
            tpPrtDauID.push_back(prt.Daughter(i_dau)-1); //IDs start at 1, we want 0 for vector indexing
        }
        fPrtDauID.   push_back(tpPrtDauID);

        fPrtNPt.        push_back(prt.NumberTrajectoryPoints());

        vector<double> tpPrtX,tpPrtY,tpPrtZ,tpPrtT,tpPrtPx,tpPrtPy,tpPrtPz,tpPrtP,tpPrtE;
        for (size_t i_ppt=0; i_ppt < prt.NumberTrajectoryPoints(); i_ppt++) {
            tpPrtX.     push_back(prt.Vx(i_ppt));
            tpPrtY.     push_back(prt.Vy(i_ppt));
            tpPrtZ.     push_back(prt.Vz(i_ppt));
            tpPrtT.     push_back(prt.T(i_ppt));
            tpPrtPx.    push_back(prt.Px(i_ppt));
            tpPrtPy.    push_back(prt.Py(i_ppt));
            tpPrtPz.    push_back(prt.Pz(i_ppt));
            tpPrtP.     push_back(prt.P(i_ppt));
            tpPrtE.     push_back(prt.E(i_ppt));
        }
        fPrtX.          push_back(tpPrtX);
        fPrtY.          push_back(tpPrtY);
        fPrtZ.          push_back(tpPrtZ);
        fPrtT.          push_back(tpPrtT);
        fPrtPx.         push_back(tpPrtPx);
        fPrtPy.         push_back(tpPrtPy);
        fPrtPz.         push_back(tpPrtPz);
        fPrtP.          push_back(tpPrtP);
        fPrtE.          push_back(tpPrtE);

        //sim::SimEnergyDeposit
        if (!tag_dep.empty()) {
        auto vh_dep = evt.getValidHandle<vector<sim::SimEnergyDeposit>>(tag_dep);

        vector<double> tpDepPdg,tpDepX,tpDepY,tpDepZ,tpDepT,tpDepE;
        for (sim::SimEnergyDeposit const & dep : *vh_dep ) {

            if (dep.TrackID()!=prt.TrackId()) continue;

            tpDepPdg.   push_back(dep.PdgCode());
            tpDepX.     push_back(dep.X());
            tpDepY.     push_back(dep.Y());
            tpDepZ.     push_back(dep.Z());
            tpDepT.     push_back(dep.T());
            tpDepE.     push_back(dep.E());

        } //end vh_dep loop
        fDepPdg.        push_back(tpDepPdg);
        fDepX.          push_back(tpDepX);
        fDepY.          push_back(tpDepY);
        fDepZ.          push_back(tpDepZ);
        fDepT.          push_back(tpDepT);
        fDepE.          push_back(tpDepE);

        fPrtNDep.       push_back(tpDepPdg.size());

        } //sim::SimEnergyDeposit
    } //end vh_prt loop
    } //end simb::MCParticle

    fTruth->Fill();

    //Reco Tree

    //recob::PFParticle
    if (!tag_pfp.empty()) { 
    art::ValidHandle<vector<recob::PFParticle>> const vh_pfp = evt.getValidHandle<vector<recob::PFParticle>>(tag_pfp);
    
    fNPfp = vh_pfp->size();

    for (recob::PFParticle const & pfp : *vh_pfp) {

        fPfpPdg.            push_back(pfp.PdgCode());

        //recob::Track
        if (!tag_trk.empty()) {
        if (!util_pfp.IsPFParticleTracklike(pfp, evt, tag_pfp.label(), tag_trk.label())) {

            fPfpTrkID.      push_back(-1);

        }
        else {
            const recob::Track* p_trk = util_pfp.GetPFParticleTrack(pfp,evt,tag_pfp.label(),tag_trk.label());

            fPfpTrkID.      push_back(p_trk->ID());

            fTrkLength.     push_back(p_trk->Length());

            vector<double>  tpTrkPtX,tpTrkPtY,tpTrkPtZ,
                            tpTrkDirX,tpTrkDirY,tpTrkDirZ;
            for (size_t i_tpt = p_trk->FirstPoint(); i_tpt < p_trk->LastPoint(); i_tpt++) {
                if (!p_trk->HasValidPoint(i_tpt)) {continue;}

                tpTrkPtX.   push_back(p_trk->LocationAtPoint(i_tpt).X());
                tpTrkPtY.   push_back(p_trk->LocationAtPoint(i_tpt).Y());
                tpTrkPtZ.   push_back(p_trk->LocationAtPoint(i_tpt).Z());
                tpTrkDirX.  push_back(p_trk->MomentumVectorAtPoint(i_tpt).X());
                tpTrkDirY.  push_back(p_trk->MomentumVectorAtPoint(i_tpt).Y());
                tpTrkDirZ.  push_back(p_trk->MomentumVectorAtPoint(i_tpt).Z());
            }
            fTrkPtX.        push_back(tpTrkPtX);
            fTrkPtY.        push_back(tpTrkPtY);
            fTrkPtZ.        push_back(tpTrkPtZ);
            fTrkDirX.       push_back(tpTrkDirX);
            fTrkDirY.       push_back(tpTrkDirY);
            fTrkDirZ.       push_back(tpTrkDirZ);

            fTrkNPt.        push_back(tpTrkPtX.size());

            //anab::Calorimetry
            if (!tag_cal.empty()) {
            vector<anab::Calorimetry> v_cal = util_trk.GetRecoTrackCalorimetry(*p_trk,evt,tag_trk.label(),tag_cal.label());

            for (anab::Calorimetry const & cal : v_cal) {
                if (cal.PlaneID().Plane != geo::kW) continue; 

                fTrkCalPlane.   push_back(cal.PlaneID().Plane);
                fTrkCalRange.   push_back(cal.Range());
                fTrkCalKineticE.push_back(cal.KineticEnergy());
                fTrkCalNPt.     push_back(cal.dEdx().size());
                // fTrkCaldEdx.    push_back(cal.dEdx());
                // fTrkCaldQdx.    push_back(cal.dQdx());
                // fTrkCalResRange.push_back(cal.ResidualRange());

                vector<double> tpTrkCaldEdx,tpTrkCaldQdx,tpTrkCalResRange;
                for (double dEdx : cal.dEdx()) {tpTrkCaldEdx.push_back(dEdx);}
                for (double dQdx : cal.dQdx()) {tpTrkCaldQdx.push_back(dQdx);}
                for (double ResR : cal.ResidualRange()) {tpTrkCalResRange.push_back(ResR);}
                fTrkCaldEdx.    push_back(tpTrkCaldEdx);
                fTrkCaldQdx.    push_back(tpTrkCaldQdx);
                fTrkCalResRange.push_back(tpTrkCalResRange);

                break; //only one Calorimetry per Track, may be unecessary but I didn't check...
            } //end v_cal loop
            } //end anab::Calorimetry
        } //end IsTrack condition
        } //end recob::Track

        //recob::Shower
        if (!tag_shw.empty()) {
        if (!util_pfp.IsPFParticleShowerlike(pfp,evt,tag_pfp.label(),tag_shw.label())) {
            
            fPfpShwID.      push_back(-1);

        }
        else {
            const recob::Shower *p_shw = util_pfp.GetPFParticleShower(pfp,evt,tag_pfp.label(),tag_shw.label());

            fPfpShwID.      push_back(p_shw->ID());

            fShwID.         push_back(p_shw->ID());

        } //end IsShower condition
        } //end recob::Shower

        //recob::Cluster
        if (!tag_clu.empty()) {
        const vector<const recob::Cluster*> vp_clu = util_pfp.GetPFParticleClusters(pfp,evt,tag_pfp.label());

        fPfpNClu.           push_back(vp_clu.size());

        vector<double> tpCluNHits,tpCluPlane,tpCluIntFit,tpCluSumADC,tpCluWidth;
        for (const recob::Cluster *p_clu : vp_clu) {

            tpCluNHits.     push_back(p_clu->NHits());
            tpCluPlane.     push_back(p_clu->Plane().Plane);
            tpCluIntFit.    push_back(p_clu->Integral());
            tpCluSumADC.    push_back(p_clu->SummedADC());
            tpCluWidth.     push_back(p_clu->Width());

        } //end vp_clu loop
        fCluNHit.           push_back(tpCluNHits);
        fCluPlane.          push_back(tpCluPlane);
        fCluIntFit.         push_back(tpCluIntFit);
        fCluSumADC.         push_back(tpCluSumADC);
        fCluWidth.          push_back(tpCluWidth);
        } //end recob::Cluster

        //recob::SpacePoint
        if (!tag_spt.empty()) {
        const vector<const recob::SpacePoint*> vp_spt = util_pfp.GetPFParticleSpacePoints(pfp,evt,tag_pfp.label());

        fPfpNSpt.           push_back(vp_spt.size());

        vector<double> tpSptX,tpSptY,tpSptZ;
        for (const recob::SpacePoint *p_spt : vp_spt) {

            tpSptX.         push_back(p_spt->position().X());
            tpSptY.         push_back(p_spt->position().Y());
            tpSptZ.         push_back(p_spt->position().Z());

        } //end vp_spt loop
        fSptX.              push_back(tpSptX);
        fSptY.              push_back(tpSptY);
        fSptZ.              push_back(tpSptZ);
        } //end recob::SpacePoint

    } //end vp_pfp loop

    if (!tag_trk.empty()) fNTrk=fTrkLength.size();
    if (!tag_shw.empty()) fNShw=fShwID.size();
    } //end recob::PFParticle
    
    fReco->Fill();

} //end analyze()

void ana::YAD::reset() {
    
    //Events
    fEvent=0;
    fRun=0;
    fSubRun=0;

    //MCTruth
    fNPrt=0;

    //MCParticle
    fPrtPdg.clear();
    fPrtMomID.clear();
    fPrtNDau.clear();
    fPrtDauID.clear();
    fPrtNPt.clear();
    fPrtX.clear();
    fPrtY.clear();
    fPrtZ.clear();
    fPrtT.clear();
    fPrtPx.clear();
    fPrtPy.clear();
    fPrtPz.clear();
    fPrtP.clear();
    fPrtE.clear();

    //SimEnergyDeposit
    fPrtNDep.clear();
    fDepPdg.clear();
    fDepX.clear();
    fDepY.clear();
    fDepZ.clear();
    fDepT.clear();
    fDepE.clear();

    //PFParticle
    fNPfp=0;
    fPfpTrkID.clear();
    fPfpShwID.clear();
    fPfpPdg.clear();

    //Track
    fNTrk=0;
    fTrkLength.clear();

    fTrkNPt.clear();
    fTrkPtX.clear();
    fTrkPtY.clear();
    fTrkPtZ.clear();

    //Calorimetry
    fTrkCalPlane.clear();
    fTrkCalRange.clear();
    fTrkCalKineticE.clear();
    fTrkCalNPt.clear();
    fTrkCaldEdx.clear();
    fTrkCaldQdx.clear();
    fTrkCalResRange.clear();

    //Shower
    fNShw=0;
    fShwID.clear();

    //Cluster
    fPfpNClu.clear();
    fCluNHit.clear();
    fCluIntFit.clear();
    fCluSumADC.clear();
    fCluWidth.clear();
    fCluPlane.clear();

    //SpacePoint
    fPfpNSpt.clear();
    fSptX.clear();
    fSptY.clear();
    fSptZ.clear();

} //end reset()

DEFINE_ART_MODULE(ana::YAD)
