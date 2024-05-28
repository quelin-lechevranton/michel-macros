#include "YAD_tools.h"

const size_t n_file=3;
const vector<string> filelist = yad::ReadFileList(n_file,"ijclab.list");

void MichelEnergy() {

    TH1D* hDep = new TH1D("hDep",";Energy (GeV?);count",10,0,0.5);
    TH1D* hCal = new TH1D("hCal",";dEdx (GeV/cm?);count",40,0,40);
    TH1D* hPrtNPt = new TH1D("hPrtNPt",";PrtNPt;count",40,60,100);
    TH1D* hMich0 = new TH1D("hMich",";E (GeV);count",50,0,0.2);
    TH1D* hMich1 = new TH1D("hMich1",";E (GeV);count",50,0,0.2);

    vector<TGraph2D*> gPrt(2);
    vector<size_t> igPrt={0,0};
    gPrt[0] = new TGraph2D();
    gPrt[0]->SetMarkerColor(kAzure+2);
    gPrt[0]->SetMarkerStyle(20);
    gPrt[0]->SetMarkerSize(0.3);

    gPrt[1] = new TGraph2D();
    gPrt[1]->SetMarkerColor(kOrange+7);
    gPrt[1]->SetMarkerStyle(20);
    gPrt[1]->SetMarkerSize(0.3);

    size_t i_file=0;
    size_t N=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {

                if (T.PrtNPt->at(i_prt) < 20) continue;
                // if (T.PrtNPt->at(i_prt) > 500) continue;
                // if (T.PrtPdg->at(i_prt)==11) {hPrtE}

                bool is_inside=true;
                for (size_t i_ppt=0; i_ppt < T.PrtNPt->at(i_prt); i_ppt++) {
                    double X = (*T.PrtX)[i_prt][i_ppt];
                    double Y = (*T.PrtY)[i_prt][i_ppt];
                    double Z = (*T.PrtZ)[i_prt][i_ppt];

                    is_inside = -350 < X && X < 350 &&
                                -350 < Y && Y < 350 &&
                                   0 < Z && Z < 300;
                    
                    if (!is_inside) break;
                }
                if (!is_inside) continue;

                if (T.PrtPdg->at(i_prt)==11) {
                    for (size_t i_ppt=0; i_ppt < T.PrtNPt->at(i_prt); i_ppt++) {
                        double X = (*T.PrtX)[i_prt][i_ppt];
                        double Y = (*T.PrtY)[i_prt][i_ppt];
                        double Z = (*T.PrtZ)[i_prt][i_ppt];

                        gPrt[0]->SetPoint(igPrt[0]++,Y,Z,X);
                        hMich0->Fill((*T.PrtE)[i_prt][i_ppt]);
                    }
                    hMich1->Fill((*T.PrtE)[i_prt][0]);
                }
                if (T.PrtPdg->at(i_prt)==13) {

                    bool is_michel = false;
                    for (size_t i_ppt=0; i_ppt < T.PrtNPt->at(i_prt); i_ppt++) {
                        double X = (*T.PrtX)[i_prt][i_ppt];
                        double Y = (*T.PrtY)[i_prt][i_ppt];
                        double Z = (*T.PrtZ)[i_prt][i_ppt];

                        gPrt[1]->SetPoint(igPrt[1]++,Y,Z,X);
                    }
                }
            }
            
            for (size_t i_dep=0; i_dep < T.NDep; i_dep++) {
                if (T.DepPdg->at(i_dep)!=11) continue;

                hDep->Fill(T.DepE->at(i_dep));
            } //end deposit loop
        } //end event loop
    } //end file loop

    // cout << N << endl;
    TCanvas* c1 = new TCanvas("c1","MichelEnergy");
    c1->Divide(2,1);
    // c1->cd();
    // hMuon->Draw("colZ");
    c1->cd(1);
    gPrt[0]->Draw("p");
    gPrt[1]->Draw("samep");
    // hPrtNPt->Draw("hist");
    // hDep->Draw("hist");
    c1->cd(2);
    hMich0->Draw("hist");
    hMich1->Draw("samehist");
    // hCal->Draw("hist");
}