#include "art/Framework/Core/EDAnalyzer.h"

/*
 * Based on Dumper module by Laura Perez Molina
 */


namespace ana {

class JDumper : public art::EDAnalyzer {
public:

    explicit JDumper(fhicl::ParameterSet const& fcl);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    JDumper(JDumper const &) = delete;
    JDumper(JDumper &&) = delete;
    JDumper & operator = (JDumper const &) = delete;
    JDumper & operator = (JDumper &&) = delete;
    
    // Required functions.
    void analyze(art::Event const & evt) override; 
    void beginJob() override;
    // void endJob()   override;
    void reset();

private:

    vector<string> fTrees;
    vector<vector<string>> fProducts;

    bool fRollUpUnsavedIDs;
    const geo::Geometry* fGeom;

    TTree *fTruth;
    TTree *fReco;
    TTree *fDetSim;

    unsigned int fEvent, fRun, fSubRun;
}

JDumper::JDumper(fhicl::ParameterSet const & fcl) :
    EDAnalyzer{p} 
{
    art::ServiceHandle<art::TFileService> tfs;
    fTrees =    fcl.get<vector<string>>("Trees");
    fProducts = fcl.get<vector<vector<string>>>("Products");



    if ( find(fTrees.begin(), fTrees.end(), "Truth") != fTress.end() ) {

        fTruth = tfs->make<TTree>("Truth","Truth");


    // BAD PRODUCT LOOP ###############################
        for (vector<string> prod : fProducts) {
            string  label=prod[0],
                    instance=prod[1],
                    object=prod[2],
                    process=prod[3];
            
            if (object == "simb::MCParticle") {

                /* ->Branch */

            } //end MCParticle
            else if (object == "sim::SimEnergyDeposit") {

                /* ->Branch */

            } //end SimEnergyDeposit
        } //end prod loop
    } //end Truth tree

    if ( find(fTrees.begin(), fTrees.end(), "DetSim") != fTress.end() ) {

        fDetSim = tfs->make<TTree>("DetSim","DetSim");

        /* ->Branch */

    } //end DetSim tree

    if ( find(fTrees.begin(), fTrees.end(), "Reco") != fTress.end() ) {

        fReco = tfs->make<TTree>("Reco","Reco");

        /* ->Branch */

    } //end Reco tree

    fRollUpUnsavedIDs = fcl.get<bool>("RollUpUnsavedIDs");
    fGeom             = &*art::ServiceHandle<geo::Geometry>();
} //end JDumper()

void JDumper::beginJob() {} //end beginJob()

void JDumper::analyze(const art::Event & evt) {
    reset();
    fEvent  = evt.id().event(); 
    fRun    = evt.id().run();
    fSubRun = evt.id().subRun();

    bool    is_truth=false,
            is_part=false,
            is_depo=false,
            is_pfp=false,
            is_clu=false,
            is_track=false,
            is_sho=false;

    string  label_truth, inst_truth,
            label_part, inst_part,
            label_depo, inst_depo,
            label_pfp, inst_pfp,
            label_clu, inst_clu,
            label_track, inst_track, 
            label_sho, inst_sho;

    for (vector<string> prod : fProducts) {

        string  label=prod[0],
                instance=prod[1],
                object=prod[2],
                process=prod[3];

        if (object == "simb::MCTruth") {
            is_truth=true;
            label_truth=label;
            inst_truth=instance;
        }
        else if (object == "simb::MCParticle") {
            is_part=true;
            label_part=label;
            inst_part=instance;
        }
        else if (object == "sim::SimEnergyDeposit") {
            is_depo=true;
            label_depo=label;
            inst_depo=instance;
        }
        else if (object == "recob::PFParticle") {
            is_pfp=true;
            label_pfp=label;
            inst_pfp=instance;
        }
        else if (object == "recob::Cluster") {
            is_clu=true;
            label_clu=label;
            inst_clu=instance;
        }
        else if (object == "recob::Track") {
            is_track=true;
            label_track=label;
            inst_track=instance;
        } 
        else if (object == "recob::Shower") {
            is_sho=true;
            label_sho=label;
            inst_sho=instance;
        }
        
    } //end prod loop

    if ( find(fTrees.begin(), fTrees.end(), "Truth") != fTress.end() ) {

        /* reinitialize vectors ??????????????????????????????? */

        if (is_truth) {
            // if(evt.isRealData()) #####################
            art::InputTag tag_truth(label_truth,inst_truth);
            auto vh_truth = evt.getValidHandle<vector<simb::MCTruth>>(tag_truth);

            /* push_backs */

            for (int i_gen=0; i_gen< vh_truth->size(); i_gen++) {

                // #################################################
                simb::MCTruth const & truth = vh_truth->at(i_gen);
                // #################################################
                simb::MCParticle const & part = truth.GetParticle(i_gen)
                // #################################################

                int i_mcpart=0;
                for (vector<string> p : fProducts) {
                    if(p[2]=="simb::MCParticle") break;
                    i_mcpart++;
                } //search for MCParticle label
                string label_mcpart = fProducts[i_mcpart][0];

                // ######################################
                art::FindManyP<simb::MCParticle> fmp_part(vh_truth,evt,label_mcpart);

                auto truth_part = fmp_part.at(i_gen);

                for (auto part_ptr : truth_part) {

                    if(part_ptr->Mother() !=0) continue;

                    /* push_backs */
                    
                } //end ?? ######################
            } //end vh_truth loop
        } //end simb::MCTruth
        
        if(is_depo) {

            art::InputTag tag_depo(label_depo,instance_depo);
            auto vh_depo = evt.getValidHandle<(tag_depo);

            for (auto depo : *vh_depo ) {

                /* push_backs */

            } //end vh_depo loop
        } //end sim::SimEnergyDeposit

        fTruth->Fill();

    } //end Truth tree

    if ( find(fTrees.begin(), fTrees.end(), "DetSim") != fTress.end() ) {

        /*  */

    } //end DetSim tree

    if ( find(fTrees.begin(), fTrees.end(), "Reco") != fTress.end() ) {

        /* .clear() reinitialize vectors */

        /* some declarations */

        if (is_pfp) {
            // const vector<art::Ptr<recob::PFParticle>> ...
            // LauraPDumper#L1066
            // dune_ana::DUNEAnaEventUtils

            // larsim/Utils/TruthMatchUtils ???

            /* ->Branch */
        } //end recob::PFParticle

        if (is_clu) {

        } //end recob::Cluster

        if (is_track) {

        } //end recob::Track

        if (is_sho) {

        } //end recob::Shower
        
        fReco->Fill();

    } //end Reco tree
} //end analyze()

void JDumper::reset() {
    //.clear() all vectors;
} //end reset()

} //end namespace ana

DEFINE_ART_MODULE(ana::JDumper)
// #############################################
