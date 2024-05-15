#include <vector>
#include <string>

#include <TChain.h>

using namespace std;

const vector<string> filelist = {
    "~/Code/out/protodunevd_10_muon_500MeV_Jdumped.root"
    // "/silver/DUNE/quelin-lechevranton/out/protodunevd_100_muon_800MeV_Jdumped.root",
    // "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_2GeV_allangles_Jdumped.root",
    // "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_1500MeV_Jdumped.root",
    // "/silver/DUNE/quelin-lechevranton/out/pdvd_100_muon_1GeV_Jdumped.root"
};

// const vector<string> filelist = {
//     "/silver/DUNE/quelin-lechevranton/out/protodunevd_10_muon_500MeV_Jdumped.root"
    // "/silver/DUNE/quelin-lechevranton/out/protodunevd_100_muon_800MeV_Jdumped.root",
    // "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_2GeV_allangles_Jdumped.root",
    // "/silver/DUNE/quelin-lechevranton/out/pdvd_1k_muon_1500MeV_Jdumped.root",
    // "/silver/DUNE/quelin-lechevranton/out/pdvd_100_muon_1GeV_Jdumped.root"
// };

void TrackEnds();
void TrackLength();

void MacroDumper() {
    // TrackEnds();
    TrackLength();
}

void TrackLength() {

    int nTracks=0;
    
    vector<double>  *TrackLength=nullptr;

    int nGraph = 3;
    vector<TH1D*> h(nGraph);

    int n_bin=100, x_min=0, x_max=700;
    h[0] = new TH1D("hTrackLength",";TrackLength (cm);count",n_bin,x_min,x_max);
    h[0]->SetLineColor(kBlack);
    h[0]->SetLineWidth(2);

    h[1] = new TH1D("hMaxLength",";TrackLength (cm);count",n_bin,x_min,x_max);
    h[1]->SetLineColor(kRed);
    h[1]->SetLineWidth(2);

    h[2] = new TH1D("hMinLength",";TrackLength (cm);count",n_bin,x_min,x_max);
    h[2]->SetLineColor(kBlue);
    h[2]->SetLineWidth(2);


    int nTrackless=0;

    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];

        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        TFile file(filename.c_str());
        TTree *Reco = (TTree*) file.Get("JDumper/Reco");
        // TTree *Truth = (TTree*) file.Get("JDumper/Truth");
        int n_ev=Reco->GetEntries();
        
        Reco->SetBranchAddress("fNTracks",      &nTracks);
        Reco->SetBranchAddress("fTrackLength",  &TrackLength);
        // Truth->SetBranchAddress("", &Y);

        for (int i_ev=0; i_ev < n_ev; i_ev++) {

            Reco->GetEntry(i_ev);
            // Truth->GetEntry(i_ev);

            if (nTracks==0) {
                nTrackless++;
            }
            else if (nTracks==1) {
                h[0]->Fill(TrackLength->at(0));
            }
            else {
                double MaxLength=0;
                double MinLength=0;

                for(int i_trk=0; i_trk < nTracks ; i_trk++) {
                    double TrLen = TrackLength->at(i_trk);

                    MaxLength > TrLen ? MaxLength : TrLen;
                    MinLength > TrLen ? TrLen : MinLength;
                
                }

                h[1]->Fill(MaxLength);
                h[2]->Fill(MinLength);
            }
        
        } //end of event loop

        file.Close();
    
    } //end of file loop

    TCanvas* canvas = new TCanvas("c","");
    canvas->cd();
    h[0]->Draw("AP");
    // h[1]->Draw("P");
    // h[2]->Draw("P");
}




void TrackEnds() {

    vector<TGraph2D*> graph(2);
    graph[0] = new TGraph2D();
    graph[0]->SetName("Track Start Positions");
    graph[0]->SetMarkerColor(kBlue);
    graph[0]->SetMarkerStyle(20);
    // graph[0]->SetMarkerSize(2);

    graph[1] = new TGraph2D();
    graph[1]->SetName("Track End Positions");
    // graph[1]->SetMarkerColor(kRed);
    graph[1]->SetMarkerStyle(20);
    // graph[1]->SetMarkerSize(2);

    graph[0]->GetXaxis()->SetTitle("X (cm)");
    graph[0]->GetYaxis()->SetTitle("Y (cm)");
    graph[0]->GetZaxis()->SetTitle("Z (cm)");

    Int_t j_total=0;

    vector<vector<double>*> TrackStart = {nullptr,nullptr,nullptr};
    vector<vector<double>*> TrackEnd = {nullptr,nullptr,nullptr};


    for (string filename : filelist) {

        cout << "\e[3mOpening file: " << filename << "\e[0m" << endl;

        TFile file(filename.c_str());
        TTree *Reco = (TTree*) file.Get("JDumper/Reco");
        Int_t n_event=Reco->GetEntries();

        Reco->SetBranchAddress("fTrackStartX", &(TrackStart[0]));
        Reco->SetBranchAddress("fTrackStartY", &(TrackStart[1]));
        Reco->SetBranchAddress("fTrackStartZ", &(TrackStart[2]));
        Reco->SetBranchAddress("fTrackEndX",   &(TrackEnd[0]));
        Reco->SetBranchAddress("fTrackEndY",   &(TrackEnd[1]));
        Reco->SetBranchAddress("fTrackEndZ",   &(TrackEnd[2]));

        cout << string(100,'-') << endl; 
        for (Int_t i_event=0; i_event < n_event; i_event++) {

            cout << "Event #" << i_event << endl; 
            Reco->GetEntry(i_event);

            for (size_t j=0; j< TrackStart[0]->size(); j++) {        

                graph[0]->SetPoint(
                    j_total,
                    TrackStart[0]->at(j),
                    TrackStart[1]->at(j),
                    TrackStart[2]->at(j)
                );

                cout << "\tTrack #" << j << endl; 
                cout << "\t\t( " << TrackStart[0]->at(j) << " , " <<
                                TrackStart[1]->at(j) << " , " <<
                                TrackStart[2]->at(j) << " ) --> ( " <<
                                TrackEnd[0]->at(j) << " , " <<
                                TrackEnd[1]->at(j) << " , " <<
                                TrackEnd[2]->at(j) << " )" << endl;

                graph[1]->SetPoint(
                    j_total++,
                    TrackEnd[0]->at(j),
                    TrackEnd[1]->at(j),
                    TrackEnd[2]->at(j)
                );
            }
            cout << string(100,'-') << endl;
        }
        file.Close();
    }

    // cout << j_total << endl;
    // cout << graph[0]->GetN() << endl;
    // cout << graph[1]->GetN() << endl;

    TCanvas* canvas = new TCanvas("c","Track Ends");
    canvas->cd();
    graph[0]->Draw("AP");
    graph[1]->Draw("P");  
}