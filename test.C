R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "tools.h"

int test(string file_name, int i_first_event, int i_last_event) {

    int pdg = 13;

    art::InputTag depo_tag("largeant:LArG4DetectorServicevolTPCActive:G4Stage1");
    art::InputTag monte_tag("generator::SinglesGen");
    
    vector<string> file_list = { file_name };
    gallery::Event ev(file_list);

    for ( ; !ev.atEnd() ; ev.next() ) {
        
        auto i_event = ev.eventAuxiliary().event();
        if ( i_event < i_first_event ) continue;
        if ( i_event > i_last_event ) break;

        
        auto const depo_list = ev.getValidHandle<vector<sim::SimEnergyDeposits>>(depo_tag);

        for (size_t i_depo=0; i_depo<depo_list->size(); i_depo++) {

            const sim::SimEnergyDeposit& depo = depo_list->at(i_depo);

            if(depo.PdgCode() != pdg) {continue;}

            auto Theta = depo.startPos.fCoordinates.Theta();

            cout << "Theta :" << Theta << endl; 


        }

        // auto const truth_list = ev.getValidHandle<vector<simb::MCTruths>>(monte_tag);
        // auto const particle_list = ev.getValidHandle<vector<simb::MCParticles>>(monte_tag);

        // for (size_t i_truth=0; i_truth<truth_list->size(); i_truth++) {

        //     const simb::MCTruths& truth = truth_list->at(i_truth);

        //     if(truth.fPartList.fpdgCode != )

        // }



    }
}

// Events->simb::MCTruths_cosmicgenerator__SinglesGen.obj.fPartList.Py()/simb::MCTruths_cosmicgenerator__SinglesGen.obj.fPartList.P()
    