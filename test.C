R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "tools.h"


int test(string file_name, int i_first_event, int i_last_event, int pdg);

int main() {
    // string file_name="/homeijclab/quelin-lechevranton/Documents/out/protodunevd_10_muon_g4_stage1.root";
    // string file_name="/eos/user/t/thoudy/pdvd/sims/out/protodunevd_10_muon_g4_stage1.root";
    string file_name="/eos/user/t/thoudy/pdvd/sims/out/protodunevd_10_muon_reco.root";
    test(file_name,1,1,13);
    return 0;
}

int test(string file_name, int i_first_event, int i_last_event, int pdg) {

    art::InputTag monte_tag("generator::SinglesGen");
    art::InputTag depo_tag("largeant:LArG4DetectorServicevolTPCActive");
    art::InputTag point_tag("pandora");
    
    vector<string> file_list = { file_name };

    // TH3D* TH_depo = new TH3D("TH_depo",  //name
    //     "SimEnergyDeposit",             //title
    //     1000,                           //n_bin_X
    //     -400,                           //X_min
    //     400,                            //X_max
    //     1000,                           //n_bin_Y
    //     -400,                           //Y_min
    //     400,                            //Y_max
    //     1000,                           //n_bin_Z
    //     -400,                           //Z_min
    //     400                             //Z_max
    // );
    
    vector<string> xtitle = {"Z (cm)","Y (cm)","Z (cm)"};
    vector<string> ytitle = {"X (cm)","X (cm)","Y (cm)"};

    vector<TGraph*> TG_depo(3);
    for(int i=0; i<TG_depo.size(); i++) {
        TG_depo[i] = new TGraph();
        TG_depo[i]->SetName("SimEnergyDeposit");
        TG_depo[i]->SetMarkerColorAlpha(kPink,.7);
        TG_depo[i]->GetXaxis()->SetTitle(xtitle[i].c_str());
        TG_depo[i]->GetYaxis()->SetTitle(ytitle[i].c_str());
    }   
    int i_depo_total = 0;

    vector<TGraph*> TG_point(3);
    for(int i=0; i<TG_point.size(); i++) {
        TG_point[i] = new TGraph();
        TG_point[i]->SetName("SpacePoint");
        TG_point[i]->SetMarkerColorAlpha(kBlue,.7);
        TG_point[i]->GetXaxis()->SetTitle(xtitle[i].c_str());
        TG_point[i]->GetYaxis()->SetTitle(ytitle[i].c_str());
    }   
    int i_depo_total = 0;
    

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

            // geo::Length_t len=depo.StepLength(); //equals 0.03 cm until the few last deposits

            
            geo::Point_t depo_point = depo.MidPoint();

            // TH_depo->Fill(depo_point.X(),depo_point.Y(),depo_point.Z());

            TG_depo[0]->SetPoint(i_depo_total,depo_point.Z(),depo_point.X());
            TG_depo[1]->SetPoint(i_depo_total,depo_point.Y(),depo_point.X());
            TG_depo[2]->SetPoint(i_depo_total,depo_point.Z(),depo_point.Y());
            i_depo_total++;

            // cout << "step length: " << len << " at event.depo: " << i_depo << i_event << endl; 




            // auto Theta = depo.startPos.fCoordinates.Theta();
            //error: 'startPos' is a private member of 'sim::SimEnergyDeposit'

            // auto len = depo.Length_t;
            //error: cannot refer to type member 'Length_t' in 'const sim::SimEnergyDeposit' with '.'


        }

        auto const point_list = ev.getValidHandle<vector<recob::SpacePoint>>(point_tag);

        for (size_t i_point=0; i_point<point_list->size(); i_point++) {

            const recob::SpacePoint& point = point_list->at(i_point);

            geo::Point_t point_point = point.position();

            TG_point[0]->SetPoint(i_point_total,point_point.Z(),point_point.X());
            TG_point[1]->SetPoint(i_point_total,point_point.Y(),point_point.X());
            TG_point[2]->SetPoint(i_point_total,point_point.Z(),point_point.Y());
            i_point_total++;



        // auto const truth_list = ev.getValidHandle<vector<simb::MCTruths>>(monte_tag);
        // auto const particle_list = ev.getValidHandle<vector<simb::MCParticles>>(monte_tag);

        // for (size_t i_truth=0; i_truth<truth_list->size(); i_truth++) {

        //     const simb::MCTruths& truth = truth_list->at(i_truth);

        //     if(truth.fPartList.fpdgCode != )

        // }



    }

    TCanvas* canvas = new TCanvas("canvas",   //name
        "SimEnergyDeposit"                      //title
    );
    // canvas->cd();
    // TH_depo->Draw();

    canvas->Divide(2,2);
    for (int i=0; i<TG_depo.size(); i++) {
        canvas->cd(i+1);
        TG_depo[i]->Draw("AP");
        TG_point[i]->Draw("P");
    }
    // canvas->Update(); //Is this useful ??

    return 0;
}

// Events->simb::MCTruths_cosmicgenerator__SinglesGen.obj.fPartList.Py()/simb::MCTruths_cosmicgenerator__SinglesGen.obj.fPartList.P()
    