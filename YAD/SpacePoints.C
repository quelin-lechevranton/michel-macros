#include "YADtools.h"

const size_t n_file=2;
const vector<string> filelist = yad::ReadFileList(n_file,"ijclab.list");

void SpacePoints() {

    vector<TGraph2D*> gTrkEnds(2);
    gTrkEnds[0] = new TGraph2D();
    gTrkEnds[0]->SetName("True Inside Track Start Positions");
    gTrkEnds[0]->SetMarkerColor(kBlue);

    gTrkEnds[1] = new TGraph2D();
    gTrkEnds[1]->SetName("True Inside Track End Positions");
    gTrkEnds[1]->SetMarkerColor(kRed);

    for (size_t i=0; i<gTrkEnds.size(); i++) {
        gTrkEnds[i]->SetMarkerStyle(20);
        gTrkEnds[i]->SetMarkerSize(0.8);
    }

    TGraph2D* gSpt = new TGraph2D;
    gSpt->SetName("Space Points");

    size_t N_evt=0;
    size_t N_trk=0;
    size_t N_trkIn=0;
    size_t iSpt=0;

    // size_t i_evt=73; 

    size_t i_file=0;
    for (string filename : filelist) {

        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Reco R(filename.c_str());
        yad::Truth T(filename.c_str());
        
        size_t n_evt=R.GetEntries();

        N_evt+=n_evt;

        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            R.GetEntry(i_evt);

            // cout << "\tEvent#" << i_evt+1 << endl;

            N_trk+=R.NTrk;

            for(size_t i_pfp=0; i_pfp < R.NPfp ; i_pfp++) {

                int i_trk = R.PfpTrkID->at(i_pfp);
                if (i_trk < 0) {continue;}

                bool is_inside = true;

                for (size_t i_spt=0; i_spt < R.PfpNSpt->at(i_pfp); i_spt++) {              

                    double X = (*R.SptX)[i_pfp][i_spt];     
                    double Y = (*R.SptY)[i_pfp][i_spt];     
                    double Z = (*R.SptZ)[i_pfp][i_spt];     
                
                    bool x_inside = -320 <= X && X <= 350;
                    bool y_inside = -317 <= Y && Y <= 317;
                    bool z_inside = 20 <= Z && Z <= 180;

                    is_inside = is_inside && x_inside && y_inside && z_inside;

                    if (!is_inside) {break;}
                } //end of spt loop 

                if (!is_inside) {continue;}

                double StX =  (*R.TrkPtX)[i_trk][0];
                double StY =  (*R.TrkPtY)[i_trk][0];
                double StZ =  (*R.TrkPtZ)[i_trk][0];
                double EndX = (*R.TrkPtX)[i_trk].back();
                double EndY = (*R.TrkPtY)[i_trk].back();
                double EndZ = (*R.TrkPtZ)[i_trk].back();
                
                for (size_t i_spt=0; i_spt < R.PfpNSpt->at(i_pfp); i_spt++) {              

                    double X = (*R.SptX)[i_pfp][i_spt];     
                    double Y = (*R.SptY)[i_pfp][i_spt];     
                    double Z = (*R.SptZ)[i_pfp][i_spt];     

                    gSpt->SetPoint(iSpt++,Y,Z,X);
                }

                int iSt=0, iEnd=1;
                if (EndX > StX) {iSt=1; iEnd=0;}

                gTrkEnds[iSt]->SetPoint(N_trkIn,
                    StY,
                    StZ,
                    StX
                );

                gTrkEnds[iEnd]->SetPoint(N_trkIn,
                    EndY,
                    EndZ,
                    EndX
                );                

                N_trkIn++;
            }
        
        } //end of event loop
    
    } //end of file loop

    cout << "nEvent=" << N_evt << endl;
    cout << "nTrack=" << N_trk << " (" << (double) N_trk/N_evt << " per event)" << endl;
    cout << "nInsde=" << N_trkIn << " (" << 100.*N_trkIn/N_trk << "%)" << endl;
    cout << "nSpPts=" << iSpt << endl;

    TCanvas* c1 = new TCanvas("c1","Track Ends");
    c1->cd();
    gSpt->Draw("p");
    gTrkEnds[0]->Draw("samep");
    gTrkEnds[1]->Draw("samep");  

    // stringstream ss;
    // ss << i_evt << ".root";
    // c1->SaveAs(ss.str().c_str());
}