R__ADD_INCLUDE_PATH("gallery/Event.h")
R__ADD_INCLUDE_PATH("lardataobj/Simulation/SimEnergyDeposit.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_types.h")
R__ADD_INCLUDE_PATH("larcoreobj/SimpleTypesAndConstants/geo_vectors.h")
R__ADD_INCLUDE_PATH("nusimdata/v1_27_01/include/nusimdata/SimulationBase")
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "tools.h"


art::InputTag monte_tag("generator::SinglesGen");
art::InputTag depo_tag("largeant:LArG4DetectorServicevolTPCActive");
art::InputTag point_tag("pandora");
art::InputTag track_tag("pandoraTrack");


int plot(vector<string> file_list, int i_first_event, int i_last_event, int pdg);

int test(int i_file=0, int i_first_event=0, int i_last_event=10, int pdg=13) {
    vector<string> files = ReadFileList(4,"file.list");
    vector<string> file_list = { files[i_file] };
    plot(file_list,i_first_event,i_last_event,pdg);
    return 0;
}


int plot(vector<string> file_list, int i_first_event, int i_last_event, int pdg=13) {


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
    
    vector<int> Xaxis = {2,1,2};
    vector<int> Yaxis = {0,0,1};
    vector<string> axis_title = {"X (cm)","Y (cm)","Z (cm)"};
    const int n_graph=Xaxis.size();

    vector<TGraph*> TG_depo(n_graph);
    for(int i=0; i<TG_depo.size(); i++) {
        TG_depo[i] = new TGraph();
        TG_depo[i]->SetName("SimEnergyDeposit");
        TG_depo[i]->SetMarkerColorAlpha(kPink,.7);
        TG_depo[i]->GetXaxis()->SetTitle(axis_title[Xaxis[i]].c_str());
        TG_depo[i]->GetYaxis()->SetTitle(axis_title[Yaxis[i]].c_str());
    }   
    int i_depo_total = 0;

    vector<TGraph*> TG_point(n_graph);
    for(int i=0; i<TG_point.size(); i++) {
        TG_point[i] = new TGraph();
        TG_point[i]->SetName("SpacePoint");
        TG_point[i]->SetMarkerColorAlpha(kBlue,.7);
        // TG_point[i]->GetXaxis()->SetTitle(axis_title[Xaxis[i]].c_str());
        // TG_point[i]->GetYaxis()->SetTitle(axis_title[Yaxis[i]].c_str());
    }   
    int i_point_total = 0;

    vector<TGraph*> TG_track(n_graph);
    for(int i=0; i<TG_track.size(); i++) {
        TG_track[i] = new TGraph();
        TG_track[i]->SetName("Track");
        TG_track[i]->SetMarkerColorAlpha(kOrange,.7);
        // TG_track[i]->GetXaxis()->SetTitle(axis_title[Xaxis[i]].c_str());
        // TG_track[i]->GetYaxis()->SetTitle(axis_title[Yaxis[i]].c_str());
    }   
    int i_track_total = 0;

    // vector<TGraph*> TG_track_valid(n_graph);
    // for(int i=0; i<TG_track_valid.size(); i++) {
    //     TG_track_valid[i] = new TGraph();
    //     TG_track_valid[i]->SetName("TrackValid");
    //     TG_track_valid[i]->SetMarkerColorAlpha(kSpring,.7);
    //     // TG_track_valid[i]->GetXaxis()->SetTitle(axis_title[Xaxis[i]].c_str());
    //     // TG_track_valid[i]->GetYaxis()->SetTitle(axis_title[Yaxis[i]].c_str());
    // }   
    // int i_track_valid_total = 0;
    

    for (
    gallery::Event ev(file_list); !ev.atEnd(); ev.next()) {
        
        auto i_event = ev.eventAuxiliary().event();
        if ( i_event < i_first_event ) continue;
        if ( i_event > i_last_event ) break;

        cout << "Event #" << i_event << endl;

        // auto const point_list = ev.getValidHandle<vector<recob::SpacePoints>>(point_tag);
    
        // cout << point_list << endl;

        /*SIM ENERGY DEPOSIT*******************/ 
        auto const depo_list = ev.getValidHandle<vector<sim::SimEnergyDeposit>>(depo_tag);

        for (size_t i_depo=0; i_depo<depo_list->size(); i_depo++) {

            const sim::SimEnergyDeposit& depo = depo_list->at(i_depo);

            if(depo.PdgCode() != pdg) {continue;} //keep only muons

            // geo::Length_t len=depo.StepLength(); //equals 0.03 cm until the few last deposits
            
            // geo::Point_t depo_point = depo.MidPoint();
            // double XYZ[3];
            // depo_point.GetCoordinates(XYZ);

            // TH_depo->Fill(depo_point.X(),depo_point.Y(),depo_point.Z());

            double XYZ[3];
            depo.MidPoint().GetCoordinates(XYZ);
            for (int i=0; i<n_graph; i++) {
                TG_depo[i]->SetPoint(i_depo_total,XYZ[Xaxis[i]],XYZ[Yaxis[i]]);
            }
            i_depo_total++;
        }
        /*END SIM ENERGY DEPOSIT************/

        /*SPACE POINTS*********************/
        auto const point_list = ev.getValidHandle<vector<recob::SpacePoint>>(point_tag);

        for (size_t i_point=0; i_point<point_list->size(); i_point++) {

            const recob::SpacePoint& point = point_list->at(i_point);

            // geo::Point_t point_point = point.position();


            double XYZ[3];
            point.position().GetCoordinates(XYZ);
            for (int i=0; i<n_graph; i++) {
                TG_depo[i]->SetPoint(i_point_total,XYZ[Xaxis[i]],XYZ[Yaxis[i]]);
            }
            i_point_total++;
        }
        /*END SPACE POINTS***************/

        /*TRACKS************************/
        auto const track_list = ev.getValidHandle<vector<recob::Track>>(track_tag);

        for (size_t i_track=0; i_track<track_list->size(); i_track++ ) {

            const recob::Track& track = track_list->at(i_track);

            for (size_t i_track_point = track.FirstPoint();
            i_track_point < track.LastPoint();
            i_track_point++ ) {

                // geo::Point_t track_point = track.LocationAtPoint(i_track_point);

                double XYZ[3];
                track.LocationAtPoint(i_track_point).GetCoordinates(XYZ);

                for (int i=0; i<n_graph; i++) {
                    TG_track[i]->SetPoint(i_track_total,XYZ[Xaxis[i]],XYZ[Yaxis[i]]);
                }
                i_track_total++;
            }

        }

        // auto const truth_list = ev.getValidHandle<vector<simb::MCTruths>>(monte_tag);
        // auto const particle_list = ev.getValidHandle<vector<simb::MCParticles>>(monte_tag);

        // for (size_t i_truth=0; i_truth<truth_list->size(); i_truth++) {

        //     const simb::MCTruths& truth = truth_list->at(i_truth);

        //     if(truth.fPartList.fpdgCode != )

        // }



    }

    TCanvas* canvas = new TCanvas("canvas",   //name
        "Bonjour"                      //title
    );
    // canvas->cd();
    // TH_depo->Draw();

    canvas->Divide(2,2);
    for (int i=0; i<3; i++) {
        canvas->cd(i+1);
        TG_depo[i]->Draw("AP");
        TG_point[i]->Draw("P");
        TG_track[i]->Draw("P");
        // TG_track_valid[i]->Draw("P");
    }
    // canvas->Update(); //Is this useful ??

    return 0;
}

// Events->simb::MCTruths_cosmicgenerator__SinglesGen.obj.fPartList.Py()/simb::MCTruths_cosmicgenerator__SinglesGen.obj.fPartList.P()
    