//#include "OmegaDumper_tools.h"
#include "OmegaLight_tools.h"

const size_t n_file=1;
const vector<string> filelist = omega::ReadFileList(n_file,"list/cosmics.list");

const bool v = true;

const omega::Limits det = omega::fiducial;
const omega::Binning bRR = {100,0,300}; //cm
const omega::Binning bdQdx = {50,0,5}; //MeV/cm
const omega::Binning bE = {100,0,100000};
const omega::Binning bBragg = {40,0,20};

const double length_mu_min = 20; //cm
const double n_cal_min = 1; //????

const double dQdx_MIP = 2; //MeV/cm
const double dQdx_min_ratio = 1;
// const double dQdx_min = dQdx_MIP*dQdx_min_ratio;

const double bragg_length = 15; //cm
const double bragg_razor = 6; //MeV/cm

const double length_el_max = 20; //cm
const double end_radius = 30; //cm

void MichelSpectrum(size_t i=0) {
    clock_t start_time=clock();

    TH1D* hE = new TH1D("hE",";SumADC;#",bE.n,bE.min,bE.max);

    TH1D* hBragg = new TH1D("hBragg",";SumADC;#",bBragg.n,bBragg.min,bBragg.max);

    TH2D* hdQdx = new TH2D(
        "hdQdx",
        "after bragg selection;residual range (cm);dE/dx (MeV/cm)",
        bRR.n,bRR.min,bRR.max,
        bdQdx.n,bdQdx.min,bdQdx.max
    );

    TGraph2D* gSpacePoint = new TGraph2D();
    gSpacePoint->SetMarkerStyle(20);
    gSpacePoint->SetMarkerSize(0.2);
    gSpacePoint->SetMarkerColor(kBlack);

    TGraph2D* gTrack = new TGraph2D();
    gTrack->SetMarkerStyle(20);
    gTrack->SetMarkerSize(0.2);
    gTrack->SetMarkerColor(kOrange+9);

    TGraph2D* gEnd = new TGraph2D();
    gEnd->SetMarkerStyle(20);
    gEnd->SetMarkerSize(0.2);
    gEnd->SetMarkerColor(kOrange+7);

    TGraph2D* gEl = new TGraph2D();
    gEl->SetMarkerStyle(20);
    gEl->SetMarkerSize(0.2);
    gEl->SetMarkerColor(kAzure);



    size_t nc=4;
    vector<size_t> c;
    while (nc--) c.push_back(0);

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Reco R(filename.c_str());

        for (size_t i_evt=0; i_evt < R.N; i_evt++) {
        // size_t i_evt=0; {

            if (v) cout << "Event#" << i_evt+1 << "/" << R.N << "\r" << flush;

            R.GetEvtPfp(i_evt);

            if (R.Trk.N < 1) continue;

            struct { double *X,*Y,*Z; } End;
            for (size_t i_pfp=0; i_pfp < R.Pfp.N; i_pfp++) {

                R.GetPfpSpt(i_pfp);
                for (size_t i_spt=0; i_spt < R.Spt.N; i_spt++) {
                    gSpacePoint->AddPoint(
                        *R.Spt.Y[i_spt],
                        *R.Spt.Z[i_spt],
                        *R.Spt.X[i_spt]
                    );
                }

                if (!R.Pfp.isTrk[i_pfp]) continue;
                R.GetPfpTrk(i_pfp);

                for (size_t i_spt=0; i_spt < R.Spt.N; i_spt++) {
                    gTrack->AddPoint(
                        *R.Spt.Y[i_spt],
                        *R.Spt.Z[i_spt],
                        *R.Spt.X[i_spt]
                    );
                }


                // if (R.Trk.Length < length_mu_min) continue;
                if (R.Cal.NPt < n_cal_min) continue;
                
                bool inside = omega::IsInside(
                    R.Trk.X,
                    R.Trk.Y,
                    R.Trk.Z,
                    det
                );

                double avg_dQdx=0;
                for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                    avg_dQdx += *R.Cal.dQdx[i_cal];
                } //end calorimetry loop
                avg_dQdx /= R.Cal.NPt;
                const double dQdx_min = avg_dQdx*dQdx_min_ratio;

                bool uprigth = *R.Trk.X[0] > *R.Trk.X.back();

                double bragg_int=0;
                size_t n_cal_bragg=0;
                if (uprigth) {
                    while (
                        *R.Cal.ResRange[n_cal_bragg++] < bragg_length 
                        && n_cal_bragg < R.Cal.NPt
                    );

                    for (size_t i_cal=0; i_cal < n_cal_bragg; i_cal++) {
                        double dQdx = *R.Cal.dQdx[i_cal]; 
                        if (dQdx < dQdx_min) continue;
                        bragg_int += dQdx;
                    }
                } else {
                    while (
                        *R.Cal.ResRange[R.Cal.NPt-1-n_cal_bragg++] > R.Cal.Range - bragg_length
                        && n_cal_bragg < R.Cal.NPt-1
                    );

                    for (size_t i_cal=R.Cal.NPt-n_cal_bragg; i_cal < R.Cal.NPt; i_cal++) {
                        double dQdx = *R.Cal.dQdx[i_cal]; 
                        if (dQdx < dQdx_min) continue;
                        bragg_int += dQdx;
                    }
                }

                hBragg->Fill(bragg_int/avg_dQdx);

                if (bragg_int/avg_dQdx < bragg_razor) continue;
                c[0]++;

                // for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                //     double dQdx = *R.Cal.dQdx[i_cal]; 
                //     double resrange;
                //     if (!upside_down) resrange = *R.Cal.ResRange[i_cal];
                //     else resrange = R.Cal.Range - *R.Cal.ResRange[i_cal];
                //     else continue;

                //     hdQdx->Fill(resrange,dQdx);
                // } //end calorimetry loop

                if (uprigth) {
                    End.X=R.Trk.X.back();
                    End.Y=R.Trk.Y.back();
                    End.Z=R.Trk.Z.back();
                    c[1]++;
                } else {
                    End.X=R.Trk.X[0];
                    End.Y=R.Trk.Y[0];
                    End.Z=R.Trk.Z[0];
                    c[2]++;
                }
                // gEnd->AddPoint(*End.Y,*End.Z,*End.X);
                for (size_t i_spt=0; i_spt < R.Spt.N; i_spt++) {
                    gEnd->AddPoint(
                        *R.Spt.Y[i_spt],
                        *R.Spt.Z[i_spt],
                        *R.Spt.X[i_spt]
                    );
                }
                
                double E=0;
                for (size_t j_pfp=0; j_pfp < R.Pfp.N; j_pfp++) {

                    if (R.Pfp.isTrk[j_pfp]) {
                        R.GetPfpTrk(j_pfp);
                        if (R.Trk.Length > length_el_max) continue;
                    }
                    R.GetPfpSpt(j_pfp);

                    for (size_t j_spt=0; j_spt < R.Spt.N; j_spt++) {
                        if (omega::Distance(
                            End.X,End.Y,End.Z,
                            R.Spt.X[j_spt],R.Spt.Y[j_spt],R.Spt.Z[j_spt]
                        ) > end_radius) continue;

                        gEl->AddPoint(*R.Spt.Y[j_spt],*R.Spt.Z[j_spt],*R.Spt.X[j_spt]);
                        R.GetSptHit(j_spt);

                        for (size_t j_hit=0; j_hit < R.Hit.N; j_hit++) {
                            if (*R.Hit.Plane[j_hit] != 2) continue;
                            E += *R.Hit.SumADC[j_hit];
                        }
                    }
                } //end michel spectrum pfparticle loop
                if (E>1) hE->Fill(E);
		        else c[3]++;
            } //end muon selection pfparticle loop
        } //end event loop
    } //end file loop
    cout << endl;

    for (size_t a : c) cout << a << endl;

    TCanvas* c0 = new TCanvas("c0","MichelSpectrum");
    c0->cd();
    hE->Draw("hist");

    TLine* l = new TLine(bragg_razor,0,bragg_razor,hBragg->GetMaximum());

    TCanvas* c1 = new TCanvas("c1","MichelSpectrum");
    c1->cd();
    hBragg->Draw("hist");
    l->Draw();
    // hdQdx->SetMinimum(1);
    // gPad->SetLogz();
    // hdQdx->Draw("colz");

    TCanvas* c2 = new TCanvas("c2","MichelSpectrum");
    c2->cd();
    gSpacePoint->Draw("p");

    TCanvas* c3 = new TCanvas("c3","MichelSpectrum");
    c3->cd();
    gTrack->Draw("p");

    TCanvas* c4 = new TCanvas("c4","MichelSpectrum");
    c4->cd();
    gEnd->Draw("p");
    // gEl->Draw("samep");

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}
