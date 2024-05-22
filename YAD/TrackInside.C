#include "YADtools.h"

const size_t n_file=2;
const vector<string> filelist = yad::ReadFileList(n_file,"ijclab.list");

void TrackInside() {

    vector<TGraph2D*> gTrkEnds(8);
    gTrkEnds[0] = new TGraph2D();
    gTrkEnds[0]->SetName("All Track Start Positions");
    gTrkEnds[0]->SetMarkerColor(kBlue);

    gTrkEnds[1] = new TGraph2D();
    gTrkEnds[1]->SetName("All Track End Positions");
    gTrkEnds[1]->SetMarkerColor(kRed);

    gTrkEnds[2] = new TGraph2D();
    gTrkEnds[2]->SetName("Inside Track Start Positions");
    gTrkEnds[2]->SetMarkerColor(kBlue);

    gTrkEnds[3] = new TGraph2D();
    gTrkEnds[3]->SetName("Inside Track End Positions");
    gTrkEnds[3]->SetMarkerColor(kRed);

    gTrkEnds[4] = new TGraph2D();
    gTrkEnds[4]->SetName("True Track Start Positions");
    gTrkEnds[4]->SetMarkerColor(kBlue);

    gTrkEnds[5] = new TGraph2D();
    gTrkEnds[5]->SetName("True Track End Positions");
    gTrkEnds[5]->SetMarkerColor(kRed);

    gTrkEnds[6] = new TGraph2D();
    gTrkEnds[6]->SetName("True Inside Track Start Positions");
    gTrkEnds[6]->SetMarkerColor(kBlue);

    gTrkEnds[7] = new TGraph2D();
    gTrkEnds[7]->SetName("True Inside Track End Positions");
    gTrkEnds[7]->SetMarkerColor(kRed);

    for (size_t i=0; i<gTrkEnds.size(); i++) {
        gTrkEnds[i]->SetMarkerStyle(20);
        gTrkEnds[i]->SetMarkerSize(0.8);
    }

    vector<TGraph2D*> gSpt(2);
    gSpt[0] = new TGraph2D;
    gSpt[0]->SetName("Space Points");
    gSpt[1] = new TGraph2D;
    gSpt[1]->SetName("Space Points");

    TGraph2D* gTrkPt = new TGraph2D();
    gTrkPt->SetName("Space Points");
    gTrkPt->SetMarkerStyle(2);
    gTrkPt->SetMarkerColor(kGreen);


    size_t N_evt=0;
    size_t N_trk=0;
    size_t N_trkIn=0;

    size_t iTrkAll=0;
    size_t iTrkTru=0;
    size_t iTrkIn=0;
    size_t iTrkTruIn=0;
    size_t iSpt=0;
    size_t iTrkPt=0;

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

                    gSpt[0]->SetPoint(iSpt++,Y,Z,X);
                
                    bool x_inside = -320 <= X && X <= 350;
                    bool y_inside = -317 <= Y && Y <= 317;
                    bool z_inside = 20 <= Z && Z <= 180;

                    is_inside = is_inside && x_inside && y_inside && z_inside;

                    if (!is_inside) {break;}
                } //end of spt loop 

                double StX =  (*R.TrkPtX)[i_trk][0];
                double StY =  (*R.TrkPtY)[i_trk][0];
                double StZ =  (*R.TrkPtZ)[i_trk][0];
                double EndX = (*R.TrkPtX)[i_trk].back();
                double EndY = (*R.TrkPtY)[i_trk].back();
                double EndZ = (*R.TrkPtZ)[i_trk].back();

                gTrkEnds[0]->SetPoint(iTrkAll,
                    StY,
                    StZ,
                    StX
                );

                gTrkEnds[1]->SetPoint(iTrkAll++,
                    EndY,
                    EndZ,
                    EndX
                );

                size_t iSt=4, iEnd=5;
                if (EndX > 300) {iSt=5; iEnd=4;}

                gTrkEnds[iSt]->SetPoint(iTrkTru,
                    StY,
                    StZ,
                    StX
                );

                gTrkEnds[iEnd]->SetPoint(iTrkTru++,
                    EndY,
                    EndZ,
                    EndX
                );

                // for (size_t i_tpt=0; i_tpt<R.TrkNPt->at(i_trk); i_tpt++) {
                //     gTrkPt->SetPoint(iTrkPt++,R.TrkPtX->at(i_trk)[i_tpt],R.TrkPtY->at(i_trk)[i_tpt],R.TrkPtZ->at(i_trk)[i_tpt]);
                // }

                if (!is_inside) {continue;}
                
                for (size_t i_spt=0; i_spt < R.PfpNSpt->at(i_pfp); i_spt++) {              

                    double X = (*R.SptX)[i_pfp][i_spt];     
                    double Y = (*R.SptY)[i_pfp][i_spt];     
                    double Z = (*R.SptZ)[i_pfp][i_spt];     

                    gSpt[1]->SetPoint(iSpt++,Y,Z,X);
                }

                gTrkEnds[2]->SetPoint(iTrkIn,
                    StY,
                    StZ,
                    StX
                );

                gTrkEnds[3]->SetPoint(iTrkIn++,
                    EndY,
                    EndZ,
                    EndX
                );
                
                iSt=6, iEnd=7;
                if (EndX > 300) {iSt=7; iEnd=6;}

                gTrkEnds[iSt]->SetPoint(iTrkTruIn,
                    StY,
                    StZ,
                    StX
                );

                gTrkEnds[iEnd]->SetPoint(iTrkTruIn++,
                    EndY,
                    EndZ,
                    EndX
                );                

                N_trkIn++;

            }
        
        } //end of event loop
    
    } //end of file loop

    cout << "nEvent=" << N_evt << endl;
    cout << "nTrack=" << N_trk << "\t\t" << iTrkAll << "\t\t" << gTrkEnds[0]->GetN() << endl;
    cout << "nInsde=" << N_trkIn << endl;
    cout << "nSpPts=" << iSpt << endl;
    cout << "nTrkPt=" << iTrkPt << endl;

    TCanvas* c1 = new TCanvas("c1","Track Inside");
    c1->Divide(2,2);
    c1->cd(1);
    // gTrkEnds[0]->GetXaxis()->SetLimits(-500,500);
    gTrkEnds[0]->Draw("p");
    gTrkEnds[1]->Draw("samep");
    gSpt[0]->Draw("samep");
    // gTrkPt->Draw("samep");
    c1->cd(2);
    gTrkEnds[2]->Draw("p");
    gTrkEnds[3]->Draw("samep");  
    gSpt[1]->Draw("samep");
    c1->cd(3);
    gTrkEnds[4]->Draw("p");
    gTrkEnds[5]->Draw("samep");  
    gSpt[0]->Draw("samep");
    c1->cd(4);
    gTrkEnds[6]->Draw("p");
    gTrkEnds[7]->Draw("samep");  
    gSpt[1]->Draw("samep");
}