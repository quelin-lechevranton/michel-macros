#include "YADtools.h"

using namespace std;

const vector<string> filelist = yad::ReadFileList(1,"jeremy.list");

void TrackPosition() {

    vector<TGraph2D*> gTrkEnds(8);
    gTrkEnds[0] = new TGraph2D();
    gTrkEnds[0]->SetName("All Track Start Positions");
    gTrkEnds[0]->SetMarkerColor(kBlue);
    gTrkEnds[0]->SetMarkerStyle(20);
    gTrkEnds[0]->SetMarkerSize(0.8);

    gTrkEnds[1] = new TGraph2D();
    gTrkEnds[1]->SetName("All Track End Positions");
    gTrkEnds[1]->SetMarkerColor(kRed);
    gTrkEnds[1]->SetMarkerStyle(20);
    gTrkEnds[1]->SetMarkerSize(0.8);

    gTrkEnds[2] = new TGraph2D();
    gTrkEnds[2]->SetName("Inside Track Start Positions");
    gTrkEnds[2]->SetMarkerColor(kBlue);
    gTrkEnds[2]->SetMarkerStyle(20);
    gTrkEnds[2]->SetMarkerSize(0.8);

    gTrkEnds[3] = new TGraph2D();
    gTrkEnds[3]->SetName("Inside Track End Positions");
    gTrkEnds[3]->SetMarkerColor(kRed);
    gTrkEnds[3]->SetMarkerStyle(20);
    gTrkEnds[3]->SetMarkerSize(0.8);

    gTrkEnds[4] = new TGraph2D();
    gTrkEnds[4]->SetName("True Track Start Positions");
    gTrkEnds[4]->SetMarkerColor(kBlue);
    gTrkEnds[4]->SetMarkerStyle(20);
    gTrkEnds[4]->SetMarkerSize(0.8);

    gTrkEnds[5] = new TGraph2D();
    gTrkEnds[5]->SetName("True Track End Positions");
    gTrkEnds[5]->SetMarkerColor(kRed);
    gTrkEnds[5]->SetMarkerStyle(20);
    gTrkEnds[5]->SetMarkerSize(0.8);

    gTrkEnds[6] = new TGraph2D();
    gTrkEnds[6]->SetName("True Inside Track Start Positions");
    gTrkEnds[6]->SetMarkerColor(kBlue);
    gTrkEnds[6]->SetMarkerStyle(20);
    gTrkEnds[6]->SetMarkerSize(0.8);

    gTrkEnds[7] = new TGraph2D();
    gTrkEnds[7]->SetName("True Inside Track End Positions");
    gTrkEnds[7]->SetMarkerColor(kRed);
    gTrkEnds[7]->SetMarkerStyle(20);
    gTrkEnds[7]->SetMarkerSize(0.8);

    TGraph2D* gSpt = new TGraph2D();
    gSpt->SetName("Space Points");

    TGraph2D* gTrkPt = new TGraph2D();
    gTrkPt->SetName("Space Points");
    gTrkPt->SetMarkerStyle(2);
    gTrkPt->SetMarkerColor(kGreen);


    int N_evt=0;
    int N_trk=0;
    int N_trkIn=0;

    int iTrkAll=0;
    int iTrkTru=0;
    int iTrkIn=0;
    int iTrkTruIn=0;
    int iSpt=0;
    int iTrkPt=0;

    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];

        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        yad::Reco R(filename.c_str());
        int n_evt=R.GetEntries();

        N_evt+=n_evt;

        for (int i_evt=0; i_evt < n_evt; i_evt++) {

            R.GetEntry(i_evt);

            N_trk+=n_trk;

            for(int i_pfp=0; i_pfp < R.NPfp ; i_pfp++) {

                int i_trk = R.PfpTrkID->at(i_pfp)
                if (i_trk < 0) {continue;}

                bool is_inside = true;

                for (int i_spt=0; i_spt < R.PfpNSpt->at(i_pfp); i_spt++) {              

                    double X = R.Spt->at(i_pfp)[i_spt].X();     
                    double Y = R.Spt->at(i_pfp)[i_spt].Y();     
                    double Z = R.Spt->at(i_pfp)[i_spt].Z();     

                    gSpt->SetPoint(iSpt++,Y,Z,X);
                
                    bool x_inside = -320 <= X && X <= 350;
                    bool y_inside = -317 <= Y && Y <= 317;
                    bool z_inside = 20 <= Z && Z <= 180;

                    is_inside = is_inside && x_inside && y_inside && z_inside;

                    if (!is_inside) {break;}
                } //end of spt loop 

                double StX =  R.TrkPt->at(i_trk)[0].X();
                double StY =  R.TrkPt->at(i_trk)[0].Y();
                double StZ =  R.TrkPt->at(i_trk)[0].Z();
                double EndX = R.TrkPt->at(i_trk).back().X();
                double EndY = R.TrkPt->at(i_trk).back().Y();
                double EndZ = R.TrkPt->at(i_trk).back().Z();

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

                int iSt=4, iEnd=5;
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

                // for (int i_pt=0; i_pt<TrkNPt->at(i_trk); i_pt++) {
                //     gTrkPt->SetPoint(iTrkPt++,TrkPt[0]->at(i_trk)[i_pt],TrkPt[1]->at(i_trk)[i_pt],TrkPt[2]->at(i_trk)[i_pt]);
                // }

                if (!is_inside) {continue;}

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

        file.Close();
    
    } //end of file loop

    cout << "nEvent=" << N_evt << endl;
    cout << "nTrack=" << N_trk << endl;
    cout << "nInsde=" << N_trkIn << endl;
    cout << "nSpPts=" << iSpt << endl;
    cout << "nTrkPt=" << iTrkPt << endl;


    TCanvas* c1 = new TCanvas("c1","Track Ends");
    c1->Divide(2,2);
    c1->cd(1);
    gTrkEnds[0]->Draw("p");
    gTrkEnds[1]->Draw("samep");
    gSpt->Draw("samep");
    // gTrkPt->Draw("samep");
    c1->cd(2);
    gTrkEnds[2]->Draw("p");
    gTrkEnds[3]->Draw("samep");  
    c1->cd(3);
    gTrkEnds[4]->Draw("p");
    gTrkEnds[5]->Draw("samep");  
    gSpt->Draw("samep");
    c1->cd(4);
    gTrkEnds[6]->Draw("p");
    gTrkEnds[7]->Draw("samep");  
}