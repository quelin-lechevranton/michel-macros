//#include "OmegaDumper_tools.h"
#include "OmegaLight_tools.h"

const size_t n_file=6;
const vector<string> filelist = omega::ReadFileList(n_file,"list/pdvd_cosmics.list");

const bool v = true;

const omega::Limits det = omega::fiducial;
const omega::Binning bRR = {100,0,300}; //cm
const omega::Binning bdEdx = {50,0,5}; //MeV/cm
const omega::Binning bE = {100,0,100000};
const omega::Binning bBragg = {62,-1,30}; //#avg dEdx

const double length_mu_min = 20; //cm
const size_t n_cal_min = 1; 
const size_t n_cal_body_min = 20;

// const double dEdx_MIP = 2; //MeV/cm
const double dEdx_min_ratio = 1.5;
// const double dQdx_min = dQdx_MIP*dQdx_min_ratio;

const double bragg_length = 2; //cm
const double bragg_razor = 10; //MeV/cm

const double length_el_max = 20; //cm
const double end_radius = 30; //cm

void MichelSpectrum_cosmics(size_t i=0) {
    clock_t start_time=clock();

    TH1D* hE = new TH1D("hE",";#AvgdEdx;#",bE.n,bE.min,bE.max);

    TH1D* hBragg = new TH1D("hBragg",";#AvgdEdx;#",bBragg.n,bBragg.min,bBragg.max);

    TH2D* hdEdx = new TH2D(
        "hdQdx",
        "after bragg selection;residual range (cm);dE/dx (MeV/cm)",
        bRR.n,bRR.min,bRR.max,
        bdEdx.n,bdEdx.min,bdEdx.max
    );

    TGraph2D* gEnd = new TGraph2D();
    gEnd->SetMarkerStyle(20);
    gEnd->SetMarkerSize(0.2);
    gEnd->SetMarkerColor(kOrange+7);

    TGraph2D* gEl = new TGraph2D();
    gEl->SetMarkerStyle(20);
    gEl->SetMarkerSize(0.2);
    gEl->SetMarkerColor(kAzure);



    size_t nc=4;
    vector<string> cn = {"tracks","bragg tracks","bragg upright","bragg reverse"};
    vector<size_t> c;
    for (int i=0; i<nc; i++) c.push_back(0);

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Reco R(filename.c_str());

        for (size_t i_evt=0; i_evt < R.N; i_evt++) {
        // size_t i_evt=0; {

            if (v) cout << "Event#" << i_evt+1 << "/" << R.N << "\r" << flush;
            if (v) cout << "evt#" << i_evt+1 << "\r" << flush;

            R.GetEvtPfp(i_evt);

            if (R.Trk.N < 1) continue;
            if (v) cout << "\tn trk: " << R.Trk.N << endl;


            size_t i_trk=1;
            struct { double *X,*Y,*Z; } End;
            for (size_t i_pfp=0; i_pfp < R.Pfp.N; i_pfp++) {
                
                R.GetPfpSpt(i_pfp);

                if (!R.Pfp.isTrk[i_pfp]) continue;
                c[0]++;

                R.GetPfpTrk(i_pfp);
               
                if (v) cout << "\ttrk#" << i_trk++ << "\r" << flush;

                // if (R.Trk.Length < length_mu_min) continue;

                if (R.Cal.NPt < n_cal_min) continue;

                if (v) cout << "\t\tn cal: " << R.Cal.NPt << endl;
                
                // bool inside = omega::IsInside(
                //     R.Trk.X,
                //     R.Trk.Y,
                //     R.Trk.Z,
                //     det
                // );

                bool uprigth = *R.Trk.X[0] > *R.Trk.X.back();
                if (v) {
                    if (uprigth) {
                        cout << "\t\tupright" << endl;
                    } else {
                        cout << "\t\tupside down" << endl;
                    }
                }

                double avg_dEdx=0;
                for (size_t i_cal=0; i_cal < R.Cal.NPt; i_cal++) {
                    avg_dEdx += *R.Cal.dEdx[i_cal];
                } //end calorimetry loop
                avg_dEdx /= R.Cal.NPt;
                const double dEdx_min = avg_dEdx*dEdx_min_ratio;

                if (v) cout << "\t\tavg dEdx: " << avg_dEdx << endl;

                double bragg_int=0;
                size_t n_cal_bragg=0;
                if (uprigth) {
                    while (
                        *R.Cal.ResRange[n_cal_bragg++] < bragg_length 
                        && n_cal_bragg < R.Cal.NPt
                    );

                    for (size_t i_cal=0; i_cal < n_cal_bragg; i_cal++) {
                        double dQdx = *R.Cal.dEdx[i_cal]; 
                        if (dQdx < dEdx_min) continue;
                        bragg_int += dQdx;
                    }
                } else {
                    while (
                        *R.Cal.ResRange[R.Cal.NPt-1-n_cal_bragg++] > R.Cal.Range - bragg_length
                        && n_cal_bragg < R.Cal.NPt-1
                    );

                    for (size_t i_cal=R.Cal.NPt-n_cal_bragg; i_cal < R.Cal.NPt; i_cal++) {
                        double dQdx = *R.Cal.dEdx[i_cal]; 
                        if (dQdx < dEdx_min) continue;
                        bragg_int += dQdx;
                    }
                }

                hBragg->Fill(bragg_int/avg_dEdx);

                if (v) cout << "\t\tbragg int: " << bragg_int << endl;

                if (bragg_int/avg_dEdx < bragg_razor) continue;
                c[1]++;

                if (v) cout << "\t\t\e[92mdead muon\e[0m" << endl;

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
                    c[2]++;
                } else {
                    End.X=R.Trk.X[0];
                    End.Y=R.Trk.Y[0];
                    End.Z=R.Trk.Z[0];
                    c[3]++;
                }


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
                hE->Fill(E);
            } //end muon selection pfparticle loop
        } //end event loop
    } //end file loop
    cout << endl;

    cout << "results: " << endl;
    for (size_t i=0; i<nc; i++) cout << "\t\t" << cn[i] << ": " << c[i] << endl;

    TCanvas* c0 = new TCanvas("c0","MichelSpectrum");
    c0->cd();
    hE->Draw("hist");

    TLine* l = new TLine(bragg_razor,0,bragg_razor,hBragg->GetMaximum());
    l->SetLineColor(kViolet);
    l->SetLineWidth(2);

    TCanvas* c1 = new TCanvas("c1","MichelSpectrum");
    c1->cd();
    gPad->SetLogy();
    hBragg->Draw("hist");
    l->Draw();
    // hdQdx->SetMinimum(1);
    // gPad->SetLogz();
    // hdQdx->Draw("colz");

    TCanvas* c2 = new TCanvas("c4","MichelSpectrum");
    c2->cd();
    gEnd->Draw("p");
    gEl->Draw("samep");
    
    c0->SaveAs("out/MichelSpectrum_cosmics_el.pdf");
    c1->SaveAs("out/MichelSpectrum_cosmics_bragg.pdf");
    c2->SaveAs("out/MichelSpectrum_cosmics_3d.pdf");

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}
