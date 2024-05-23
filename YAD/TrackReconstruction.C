#include "YADtools.h"

const size_t n_file=1;
const vector<string> filelist = yad::ReadFileList(n_file,"ijclab.list");

double ZenithAngle(double x, double m) {
    return TMath::RadToDeg()*TMath::ACos(x/m);
}

void TrackReconstruction() {

    // Track Length ==================================================
    vector<TH1D*> hTrkLen(4);
    int binTrkLen=50, minTrkLen=0, maxTrkLen=300;

    hTrkLen[0] = new TH1D("hTrackoneLength","",binTrkLen,minTrkLen,maxTrkLen);
    hTrkLen[0]->SetLineColor(kGreen+2);
    hTrkLen[0]->SetLineWidth(2);

    hTrkLen[1] = new TH1D("hMaxLength","",binTrkLen,minTrkLen,maxTrkLen);
    hTrkLen[1]->SetLineColor(kRed+1);
    hTrkLen[1]->SetLineWidth(2);

    hTrkLen[2] = new TH1D("hSumLength","",binTrkLen,minTrkLen,maxTrkLen);
    hTrkLen[2]->SetLineColor(kViolet-5);
    hTrkLen[2]->SetLineWidth(2);

    hTrkLen[3] = new TH1D("hAllLength","",binTrkLen,minTrkLen,maxTrkLen);
    hTrkLen[3]->SetLineColor(kBlack);
    // hTrkLen[3]->SetLineWidth(2);

    THStack* hsTrkLen = new THStack("hsTrkLen","All Tracks - #color[209]{single-track events} - #color[99]{multi-track events (max)} - #color[51]{multi-track events (sum)};Track Length (cm);count");

    for (int i=0; i<3; i++) {hsTrkLen->Add(hTrkLen[i]);}

    // Start Momentum ================================================
    vector<TH1D*> hStartP(4);
    int binStP=40, minStP=0, maxStP=3;

    hStartP[0] = new TH1D("hTrackeverTrueMomentum","",binStP,minStP,maxStP);
    hStartP[0]->SetLineColor(kBlack);
    hStartP[0]->SetLineWidth(1);

    hStartP[1] = new TH1D("hTracklessTrueMomentum","",binStP,minStP,maxStP);
    hStartP[1]->SetLineColor(kBlue-4);
    hStartP[1]->SetLineWidth(2);

    hStartP[2] = new TH1D("hTrackoneTrueMomentum","",binStP,minStP,maxStP);
    hStartP[2]->SetLineColor(kGreen+2);
    hStartP[2]->SetLineWidth(2);

    hStartP[3] = new TH1D("hTrackfulTrueMomentum","",binStP,minStP,maxStP);
    hStartP[3]->SetLineColor(kRed+1);
    hStartP[3]->SetLineWidth(2);

    THStack* hsStartP = new THStack("hsStartP","All Tracks - #color[209]{single-track events} - #color[99]{multi-track events} - #color[57]{trackless events};Generated Momentum (GeV);count");
    for (int i=0; i<hStartP.size(); i++) {hsStartP->Add(hStartP[i]);}
 
    // Zenith Angle ==================================================
    vector<TH1D*> hZen(4);
    int binZen=30, minZen=0, maxZen=90;

    hZen[0] = new TH1D("hTrackeverTrueZen","",binZen,minZen,maxZen);
    hZen[0]->SetLineColor(kBlack);
    hZen[0]->SetLineWidth(1);

    hZen[1] = new TH1D("hTracklessTrueZen","",binZen,minZen,maxZen);
    hZen[1]->SetLineColor(kBlue-4);
    hZen[1]->SetLineWidth(2);

    hZen[2] = new TH1D("hTrackoneTrueZen","",binZen,minZen,maxZen);
    hZen[2]->SetLineColor(kGreen+2);
    hZen[2]->SetLineWidth(2);

    hZen[3] = new TH1D("hTrackfulTrueZen","",binZen,minZen,maxZen);
    hZen[3]->SetLineColor(kRed+1);
    hZen[3]->SetLineWidth(2);

    THStack* hsZen = new THStack("hsZen","All Tracks - #color[209]{single-track events} - #color[99]{multi-track events} - #color[57]{trackless events};Generated Zenith (deg);count");
    for (int i=0; i<hZen.size(); i++) {hsZen->Add(hZen[i]);}

    // TrkLen % StartP ===============================================
    TH2D* gTrkLen_StP = new TH2D("hTrkLen_StP",";Generated Momentum (GeV);Track Length (cm)",binStP,minStP,maxStP,binTrkLen,minTrkLen,maxTrkLen);

    // Zenith Angle % StartP =========================================
    vector<TH2D*> hZen_StP(4);
    hZen_StP[0] = new TH2D("hZen_StP0","All Events;Generated Momentum (GeV);Generated Zenith (deg)",binStP,minStP,maxStP,binZen,minZen,maxZen);
    hZen_StP[1] = new TH2D("hZen_StP1","Trackless Events;Generated Momentum (GeV);Generated Zenith (deg)",binStP,minStP,maxStP,binZen,minZen,maxZen);
    hZen_StP[2] = new TH2D("hZen_StP2","Trackone Events;Generated Momentum (GeV);Generated Zenith (deg)",binStP,minStP,maxStP,binZen,minZen,maxZen);
    hZen_StP[3] = new TH2D("hZen_StP3","Trackful Events;Generated Momentum (GeV);Generated Zenith (deg)",binStP,minStP,maxStP,binZen,minZen,maxZen);

    int nEvent=0;
    int nTrackless=0;
    int nTrackful=0;

    size_t i_file=0;
    for (string filename : filelist) {

        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Reco R(filename.c_str());
        yad::Truth T(filename.c_str());

        int n_evt=R.GetEntries();
        nEvent+=n_evt;


        for (int i_evt=0; i_evt < n_evt; i_evt++) {

            R.GetEntry(i_evt);
            T.GetEntry(i_evt);

            double TruStPx= T.PrtPx->at(0)[0];
            double TruStPy= T.PrtPy->at(0)[0];
            double TruStPz= T.PrtPz->at(0)[0];
            double TruStP = TMath::Sqrt(TruStPx*TruStPx + TruStPy*TruStPy+ TruStPz*TruStPz);
            double TruZen = ZenithAngle(-TruStPx,TruStP);

            hStartP[0]->Fill(TruStP);
            hZen[0]->Fill(TruZen);

            hZen_StP[0]->Fill(TruStP,TruZen);

            if (R.NTrk==0) {
                nTrackless++;

                hStartP[1]->Fill(TruStP);

                hZen[1]->Fill(TruZen);

                // gTrkLen_StP->Fill(TruStP,0.0);
                hZen_StP[1]->Fill(TruStP,TruZen);
            }
            else if (R.NTrk==1) {

                double TrkLen = R.TrkLength->at(0);

                hTrkLen[0]->Fill(TrkLen); 
                hTrkLen[3]->Fill(TrkLen); 

                hStartP[2]->Fill(TruStP);

                hZen[2]->Fill(TruZen);

                gTrkLen_StP->Fill(TruStP,TrkLen);

                hZen_StP[2]->Fill(TruStP,TruZen);
            }
            else {
                nTrackful++;

                double MaxLen=0;
                double MinLen=0;
                double SumLen=0;

                for(int i_trk=0; i_trk < R.NTrk ; i_trk++) {

                    double TrkLen = R.TrkLength->at(i_trk);

                    MaxLen = MaxLen > TrkLen ? MaxLen : TrkLen;
                    MinLen = MinLen > TrkLen ? TrkLen : MinLen;
                    SumLen+= TrkLen;

                    hTrkLen[3]->Fill(TrkLen);

                }

                hTrkLen[1]->Fill(MaxLen);
                hTrkLen[2]->Fill(SumLen);

                hStartP[3]->Fill(TruStP);

                hZen[3]->Fill(TruZen);

                // gTrkLen_StP->Fill(TruStP,SumLen);
                
                hZen_StP[3]->Fill(TruStP,TruZen);
            }
        
        } //end of event loop

    
    } //end of file loop

    cout << "nEvent=" << nEvent << endl;
    cout << "nTrackful=" << nTrackful << " (" << 100.*nTrackful/nEvent << "%)" << endl;
    cout << "nTrackless=" << nTrackless << " (" << 100.*nTrackless/nEvent << "%)" << endl;

    TCanvas* c1 = new TCanvas("c1","Track Length1");
    c1->Divide(2,2);
    c1->cd(1);
    hsTrkLen->Draw("nostack");

    c1->cd(2);
    hsStartP->Draw("nostack");

    c1->cd(3);
    hsZen->Draw("nostack");

    c1->cd(4);
    gTrkLen_StP->Draw("colZ");

    TCanvas* c2 = new TCanvas("c2","Track Length2");
    c2->Divide(2,2);
    c2->cd(1);
    hZen_StP[0]->Draw("colZ");

    c2->cd(2);
    hZen_StP[1]->Draw("colZ");

    c2->cd(3);
    hZen_StP[2]->Draw("colZ");

    c2->cd(4);
    hZen_StP[3]->Draw("colZ");
    // c1->SaveAs("R.TrkLength.pdf");
}