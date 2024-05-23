#include "YADtools.h"

const size_t n_file=1;
const vector<string> filelist = yad::ReadFileList(n_file,"jeremy.list");

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

    vector<TGraph2D*> gSpt(2);
    gSpt[0] = new TGraph2D();
    gSpt[0]->SetName("Track Space Points");
    gSpt[0]->SetMarkerColor(kViolet-7);
    gSpt[0]->SetMarkerStyle(20);
    gSpt[0]->SetMarkerSize(0.3);

    gSpt[1] = new TGraph2D();
    gSpt[1]->SetName("non-Track Space Points");

    vector<TGraph2D*> gDep(3);
    gDep[0] = new TGraph2D();
    gDep[0]->SetName("Muon deposits");
    gDep[0]->SetMarkerColor(kOrange+7);
    gDep[0]->SetMarkerStyle(20);
    gDep[0]->SetMarkerSize(0.3);

    gDep[1] = new TGraph2D();
    gDep[1]->SetName("Electron deposits");
    gDep[1]->SetMarkerColor(kBlue-4);
    gDep[1]->SetMarkerStyle(20);
    gDep[1]->SetMarkerSize(0.1);

    gDep[2] = new TGraph2D();
    gDep[2]->SetName("Other deposits");

    size_t N_evt=0;
    size_t N_trk=0;
    size_t N_trkIn=0;
    size_t iSpt=0;
    size_t iSpt1=0;
    size_t iDep0=0;
    size_t iDep1=0;
    size_t iDep2=0;

    // size_t i_evt=73; 

    size_t i_file=0;
    for (string filename : filelist) {

        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Reco R(filename.c_str());
        yad::Truth T(filename.c_str());
        
        size_t n_evt=R.GetEntries();

        N_evt+=n_evt;

        // for (size_t i_evt=0; i_evt < n_evt; i_evt++) {
        size_t i_evt=45; {

            R.GetEntry(i_evt);
            T.GetEntry(i_evt);

            for (size_t i_dep=0; i_dep< T.DepX->size(); i_dep++) {

                double X = (*T.DepX)[i_dep];
                double Y = (*T.DepY)[i_dep];
                double Z = (*T.DepZ)[i_dep];

                if (T.DepPdg->at(i_dep)==13) {
                    gDep[0]->SetPoint(iDep0++,Y,Z,X);
                } else if (T.DepPdg->at(i_dep) == 11) {
                    gDep[1]->SetPoint(iDep1++,Y,Z,X);
                } else {
                    gDep[2]->SetPoint(iDep2++,Y,Z,X);
                }
            }


            N_trk+=R.NTrk;

            for(size_t i_pfp=0; i_pfp < R.NPfp ; i_pfp++) {

                int i_trk = R.PfpTrkID->at(i_pfp);
                
                if (i_trk < 0) {
                    for (size_t i_spt=0; i_spt < R.PfpNSpt->at(i_pfp); i_spt++) {              

                    double X = (*R.SptX)[i_pfp][i_spt];     
                    double Y = (*R.SptY)[i_pfp][i_spt];     
                    double Z = (*R.SptZ)[i_pfp][i_spt];     
                    
                    gSpt[1]->SetPoint(iSpt1++,Y,Z,X);
                }
                    
                }
                else {
                bool is_inside = true;

                // if (R.TrkLength->at(i_trk) > 15) {continue;}
                // if (i_trk != 0) {continue;}

                for (size_t i_spt=0; i_spt < R.PfpNSpt->at(i_pfp); i_spt++) {              

                    double X = (*R.SptX)[i_pfp][i_spt];     
                    double Y = (*R.SptY)[i_pfp][i_spt];     
                    double Z = (*R.SptZ)[i_pfp][i_spt];     
                
                    bool x_inside = -320 <= X && X <= 350;
                    bool y_inside = -317 <= Y && Y <= 317;
                    bool z_inside = 20 <= Z && Z <= 280;

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

                    gSpt[0]->SetPoint(iSpt++,Y,Z,X);
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
            }
        
        } //end of event loop
    
    } //end of file loop

    // cout << "nEvent=" << N_evt << endl;
    // cout << "nTrack=" << N_trk << " (" << (double) N_trk/N_evt << " per event)" << endl;
    // cout << "nInsde=" << N_trkIn << " (" << 100.*N_trkIn/N_trk << "%)" << endl;
    // cout << "nSpPts=" << iSpt << endl;

    TCanvas* c1 = new TCanvas("c1","Track Ends");
    c1->Divide(2,1);
    c1->cd(1);
    gSpt[0]->Draw("p");
    // gSpt[1]->Draw("samep");
    gTrkEnds[0]->Draw("samep");
    gTrkEnds[1]->Draw("samep");  

    c1->cd(2);
    gDep[0]->Draw("p");
    gDep[1]->Draw("samep");
    // gDep[2]->Draw("samep");

    // stringstream ss;
    // ss << i_evt << ".root";
    // c1->SaveAs(ss.str().c_str());
}