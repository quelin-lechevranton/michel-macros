R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "tools.h"

int test(string a, int b, int c);

int main() {
    string file_name="/eos/user/t/thoudy/pdvd/sims/out/protodunevd_10_muon_reco.root";
    test(file_name,1,10);
    return 0;
}

int test(string file_name, int i_first_event, int i_last_event) {

    int pdg = 13;

    art::InputTag monte_tag("generator::SinglesGen");
    art::InputTag depo_tag("largeant:LArG4DetectorServicevolTPCActive");
    art::InputTag point_tag("pandora");
    
    vector<string> file_list = { file_name };

    for (
    gallery::Event ev(file_list); !ev.atEnd(); ev.next()) {
        
        auto i_event = ev.eventAuxiliary().event();
        if ( i_event < i_first_event ) continue;
        if ( i_event > i_last_event ) break;

        cout << "Event #" << i_event << endl;

        // auto const point_list = ev.getValidHandle<vector<recob::SpacePoints>>(point_tag);
    
        // cout << point_list << endl;

        
        auto const depo_list = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(depo_tag);

        for (size_t i_depo=0; i_depo<depo_list->size(); i_depo++) {

            const sim::SimEnergyDeposit& depo = depo_list->at(i_depo);

            if(depo.PdgCode() != pdg) {continue;}

            // auto Theta = depo.startPos.fCoordinates.Theta();
            //error: /afs/cern.ch/work/t/thoudy/DUNE/analyses/ProtoDUNE/michel/test.C:43:31: error: 'startPos' is a private member of 'sim::SimEnergyDeposit'

            auto len = depo.Length_t;


            cout << "length: " << len << "at event.depo: " << i_depo << i_event << endl; 


        }

        // auto const truth_list = ev.getValidHandle<vector<simb::MCTruths>>(monte_tag);
        // auto const particle_list = ev.getValidHandle<vector<simb::MCParticles>>(monte_tag);

        // for (size_t i_truth=0; i_truth<truth_list->size(); i_truth++) {

        //     const simb::MCTruths& truth = truth_list->at(i_truth);

        //     if(truth.fPartList.fpdgCode != )

        // }



    }

    return 0;
}

// Events->simb::MCTruths_cosmicgenerator__SinglesGen.obj.fPartList.Py()/simb::MCTruths_cosmicgenerator__SinglesGen.obj.fPartList.P()
    