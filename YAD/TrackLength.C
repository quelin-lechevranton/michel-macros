
#include "YADtools.h"

const vector<string> filelist = yad::ReadFileList(1,"jeremy.list");

double ZenithAngle(double x, double m) {
    return TMath::RadToDeg()*TMath::ACos(x/m);
}

void TrackLength() {

    // Track Length ==================================================
    vector<TH1D*> hTrkLen(4);
    int binTrkLen=30, minTrkLen=0, maxTrkLen=700;

    hTrkLen[0] = new TH1D("hTrackoneLength","",binTrkLen,minTrkLen,maxTrkLen);
    hTrkLen[0]->SetLineColor(kGreen+2);
    hTrkLen[0]->SetLineWidth(2);

    hTrkLen[1] = new TH1D("hMaxLength","",binTrkLen,minTrkLen,maxTrkLen);
    hTrkLen[1]->SetLineColor(kRed+1);
    hTrkLen[1]->SetLineWidth(2);

    hTrkLen[2] = new TH1D("hMinLength","",binTrkLen,minTrkLen,maxTrkLen);
    hTrkLen[2]->SetLineColor(kViolet-5);
    hTrkLen[2]->SetLineWidth(2);

    hTrkLen[3] = new TH1D("hSumLength","",binTrkLen,minTrkLen,maxTrkLen);
    hTrkLen[3]->SetLineColor(kViolet-5);
    hTrkLen[3]->SetLineWidth(2);

    THStack* hsTrkLen = new THStack("hsTrkLen",";TrackLength (cm);count");
    for (int i=0; i<hTrkLen.size(); i++) {if (i==2) {continue;} hsTrkLen->Add(hTrkLen[i]);}

    // Start Momentum ================================================
    vector<TH1D*> hStartP(4);
    int binStP=30, minStP=0, maxStP=4;

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

    THStack* hsStartP = new THStack("hsStartP",";Momentum (GeV);count");
    for (int i=0; i<hStartP.size(); i++) {hsStartP->Add(hStartP[i]);}
 
    // Zenith Angle ==================================================
    vector<TH1D*> hZen(5);
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

    hZen[3] = new TH1D("hTrackoneZen","",binZen,minZen,maxZen);
    hZen[3]->SetLineColor(kGreen+2);
    hZen[3]->SetLineWidth(2);
    hZen[3]->SetLineStyle(kDashed);

    hZen[4] = new TH1D("hTrackfulTrueZen","",binZen,minZen,maxZen);
    hZen[4]->SetLineColor(kRed+1);
    hZen[4]->SetLineWidth(2);

    // hZen[4] = new TH1D("hTrackfulMaxZen","",binZen,minZen,maxZen);
    // hZen[4]->SetLineColor(kOrange+7);

    THStack* hsZen = new THStack("hsZen",";Zen (deg);count");
    for (int i=0; i<hZen.size(); i++) {hsZen->Add(hZen[i]);}

    // TrkLen % StartP ===============================================
    TH2D* gTrkLen_StP = new TH2D("hTrkLen_StP",";Momentum (GeV);Track Length (cm)",binStP,minStP,maxStP,binTrkLen,minTrkLen,maxTrkLen);

    // int iTrkLen_StP=0;

    int nEvent=0;
    int nTrackless=0;
    int nTrackful=0;

    for (int i_file=0; i_file < filelist.size(); i_file++) {
        string filename = filelist[i_file];

        cout << "\e[3mOpening file #" << i_file+1 << ": " << filename << "\e[0m" << endl;

        yad::Reco R(filename.c_str());
        yad::Truth T(filename.c_str());

        int n_ev=R.GetEntries();
        nEvent+=n_ev;


        for (int i_ev=0; i_ev < n_ev; i_ev++) {

            R.GetEntry(i_ev);
            T.GetEntry(i_ev);

            double TruStPx= T.PrtStP[0].Px();
            double TruStP = T.PrtStP[0].P();
            double TruZen = ZenithAngle(-TruStPx,TruStP);

            int nTracks=R.NTrk;
            vector<double>* TrackLength=R.TrkLength;

            hStartP[0]->Fill(TruStP);
            hZen[0]->Fill(TruZen);

            if (nTracks==0) {
                nTrackless++;

                hStartP[1]->Fill(TruStP);

                hZen[1]->Fill(TruZen);

                // gTrkLen_StP->Fill(TruStP,0.0);
            }
            else if (nTracks==1) {

                double TrkLen = TrackLength->at(0);
                double TrkStDirX= R.TrkDir[0][0].X();
                double TrkZen = ZenithAngle(-TrkStDirX,1.);

                hTrkLen[0]->Fill(TrkLen); 

                hStartP[2]->Fill(TruStP);

                hZen[2]->Fill(TruZen);
                hZen[3]->Fill(TrkZen);

                gTrkLen_StP->Fill(TruStP,TrkLen);

            }
            else {
                nTrackful++;

                double MaxLen=0;
                double MinLen=0;
                double SumLen=0;

                for(int i_trk=0; i_trk < nTracks ; i_trk++) {

                    double TrkLen = TrackLength->at(i_trk);

                    MaxLen = MaxLen > TrkLen ? MaxLen : TrkLen;
                    MinLen = MinLen > TrkLen ? TrkLen : MinLen;
                    SumLen+= TrkLen;

                }

                hTrkLen[1]->Fill(MaxLen);
                hTrkLen[2]->Fill(MinLen);
                hTrkLen[3]->Fill(SumLen);

                hStartP[3]->Fill(TruStP);

                hZen[4]->Fill(TruZen);

                gTrkLen_StP->Fill(TruStP,SumLen);
            }
        
        } //end of event loop

    
    } //end of file loop

    cout << "nEvent=" << nEvent << endl;
    cout << "nTrackful=" << nTrackful << endl;
    cout << "nTrackless=" << nTrackless << endl;

    TCanvas* c1 = new TCanvas("c1","Track Length");
    c1->Divide(2,2);
    c1->cd(1);
    hsTrkLen->Draw("nostack");

    c1->cd(2);
    hsStartP->Draw("nostack");

    c1->cd(3);
    hsZen->Draw("nostack");

    c1->cd(4);
    gTrkLen_StP->Draw("colZ");

    // c1->SaveAs("TrackLength.pdf");
}