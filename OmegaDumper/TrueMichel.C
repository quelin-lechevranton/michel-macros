#include "OmegaLight_tools.h"

const size_t n_file=37;
const vector<string> filelist = omega::ReadFileList(n_file,"list/muplus.list");

const bool v = false;

const double dep_step = 0.03;

const double r_coincidence_rough =1.; //cm
const size_t step_rough = 10;
const size_t window_rough = 4;
const double r_coincidence =1.; //cm

const size_t n_least_deposits=20;
const size_t n_dep_bragg=20;
const double dEdx_min_ratio = 1;
const double bragg_razor=12;

const omega::Limits det = omega::fiducial;

const omega::Binning bRR = {100,0,200};
const omega::Binning bdEdx = {50,0,5};
const omega::Binning bBragg = {100,0,20}; 
const omega::Binning bMichel = {35,0,70};

vector<size_t> n_c = {0,0,0,0,0,0,0,0,0};
#define COND(i) {n_c[i]++; if (v) cout << "\t\t\e[91mc" << i << "\e[0m" << endl; continue;} i_c++;

#define FILLHBRAGG2D(i) for (size_t i_dep=0; i_dep <= i_max ; i_dep++) {hBragg2D[i]->Fill((i_max-i_dep)*dep_step, *T.Dep.E[i_dep]/dep_step);}

void TrueMichel() {
    clock_t start_time=clock();

    TColor color;

    vector<TH2D*> hBragg2D(2);
    hBragg2D[0] = new TH2D(
        "hBragg2D0",
        "Muon Energy Loss (before selection);Residual Range (cm);dEdx (MeV/cm)",
        bRR.n, bRR.min, bRR.max,
        bdEdx.n, bdEdx.min, bdEdx.max
    );
    hBragg2D[1] = new TH2D(
        "hBragg2D1",
        "Muon Energy Loss (after selection);Residual Range (cm);dEdx (MeV/cm)",
        bRR.n, bRR.min, bRR.max,
        bdEdx.n, bdEdx.min, bdEdx.max
    );
    vector<TH1D*> hBragg(2);
    hBragg[0] = new TH1D(
        "hBragg0",
        "before selection;avg dEdx per deposit (MeV/cm);#",
        bBragg.n, bBragg.min, bBragg.max
    );
    hBragg[1] = new TH1D(
        "hBragg1",
        "after selection;avg dEdx per deposit (MeV/cm);#",
        bBragg.n, bBragg.min, bBragg.max
    );
    TH1D* hMichel = new TH1D(
        "hMichel",
        "Michel Spectrum;Total Deposited Energy (MeV);#",
        bMichel.n, bMichel.min, bMichel.max
    );

    size_t  N_evt=0,
            N_prt=0;
    size_t  N_mu_out=0;
    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Truth T(filename.c_str());

        N_evt+=T.N;
        for (size_t i_evt=0; i_evt < T.N; i_evt++) {

            if (!v) cout << "Event#" << i_evt+1 << "/" << T.N << "\r" << flush;
            if (v) cout << "evt#" << i_evt+1 << "\r";

            T.GetEvt(i_evt);

            N_prt+=T.Prt.N;

            for (size_t i_prt=0; i_prt < T.Prt.N; i_prt++) {
                T.GetPrt(i_prt);
                if (T.Prt.Pdg!=13 && T.Prt.Pdg!=-13) continue;
                if (!omega::IsInside(T.Prt.X,T.Prt.Y,T.Prt.Z,det)) N_mu_out++;
            }

            for (size_t i_prt=0; i_prt < T.Prt.N; i_prt++) {

                if (v) cout << "\tprt# " << i_prt << "\r" << flush;
                size_t i_c=0;

                T.GetPrt(i_prt);
                T.GetPrtDep(i_prt);

                if (T.Prt.Pdg != 11 && T.Prt.Pdg != -11) COND(i_c)
                if (T.Dep.N < n_least_deposits) COND(i_c)
                // bool inside = omega::IsInside(T.Dep.X, T.Dep.Y, T.Dep.Z, det);
                // if (!inside) COND(i_c)

                if (T.Prt.isOrphelin) COND(i_c) 
                T.GetPrtMom(i_prt);
                T.GetMomDep(i_prt);
                bool inside = omega::IsInside(T.Dep.X, T.Dep.Y, T.Dep.Z, det);
                if (!inside) COND(i_c)
                if (T.Mom.Pdg != 13 && T.Mom.Pdg != -13) COND(i_c)

                size_t n_max_rough=0, i_max_rough=0;
                for (size_t i_dep=0; i_dep<T.Dep.N; i_dep+=step_rough) {
                    size_t n_coincidence=0;
                    for (size_t i_ppt=0; i_ppt<T.Prt.NPt; i_ppt+=step_rough) {
                        double dist = omega::Distance(
                            T.Prt.X[i_ppt],T.Prt.Y[i_ppt],T.Prt.Z[i_ppt],
                            T.Dep.X[i_dep],T.Dep.Y[i_dep],T.Dep.Z[i_dep]
                        );
                        if (dist < r_coincidence_rough) n_coincidence++;
                    }
                    if (n_coincidence >= n_max_rough) {
                        n_max_rough = n_coincidence;
                        i_max_rough = i_dep;
                    }
                }

                int i_dep_ini = i_max_rough-window_rough*step_rough;
                i_dep_ini = i_dep_ini > 0 ? i_dep_ini : 0;
                int i_dep_end = i_max_rough+window_rough*step_rough;
                i_dep_end = i_dep_end > T.Dep.N ? T.Dep.N : i_dep_end;
                size_t n_max=0, i_max=0;
                for (size_t i_dep=i_dep_ini; i_dep<i_dep_end; i_dep++) {
                    size_t n_coincidence=0;
                    for (size_t i_ppt=0; i_ppt<T.Prt.NPt; i_ppt++) {
                        double dist = omega::Distance(
                            T.Prt.X[i_ppt],T.Prt.Y[i_ppt],T.Prt.Z[i_ppt],
                            T.Dep.X[i_dep],T.Dep.Y[i_dep],T.Dep.Z[i_dep]
                        );
                        if (dist < r_coincidence) n_coincidence++;
                    }
                    if (n_coincidence >= n_max) {
                        n_max = n_coincidence;
                        i_max = i_dep;
                    }
                }
                
                int n_dep_body = i_max - n_dep_bragg;
                if (n_dep_body <= 0) COND(i_c)
                double avg_body_dEdx=0;
                for (size_t i_dep=0; i_dep<n_dep_body; i_dep++) {
                    avg_body_dEdx += *T.Dep.E[i_dep];
                }
                avg_body_dEdx /= dep_step * n_dep_body;
                double bragg_int=0;
                for (size_t i_dep=n_dep_body; i_dep < i_max; i_dep++) {
                    double dEdx = *T.Dep.E[i_dep] / dep_step;
                    bragg_int += dEdx;
                }

                hBragg[0]->Fill(bragg_int/n_dep_bragg);
                FILLHBRAGG2D(0)

                if (bragg_int/n_dep_bragg < bragg_razor) COND(i_c)

                hBragg[1]->Fill(bragg_int/n_dep_bragg);
                FILLHBRAGG2D(1)

                T.GetPrtDep(i_prt);
                double E=0;
                for (size_t i_dep=0; i_dep < T.Dep.N ; i_dep++) {
                    E+=*T.Dep.E[i_dep];
                }
                hMichel->Fill(E);

            } //end particle loop
        } //end event loop
    } //end file loop
    cout << endl;

    TLine* l = new TLine(bragg_razor,0,bragg_razor,hBragg[0]->GetMaximum());
    l->SetLineColor(kViolet); 
    l->SetLineWidth(2);

    TCanvas* c1 = new TCanvas("c1","TrueMichel - Bragg");
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetLogz();
    hBragg2D[0]->Draw("colz");
    c1->cd(2);
    gPad->SetLogz();
    hBragg2D[1]->Draw("colz");
    c1->cd(3);
    gPad->SetLogy();
    hBragg[0]->Draw("hist");
    l->Draw();
    c1->cd(4);
    hBragg[1]->Draw("hist");

    TCanvas* c2 = new TCanvas("c2","TrueMichel - Michel Spectrum");
    c2->cd();
    hMichel->Draw("hist");


    cout << N_evt << " events treated in " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;

    size_t N_mu_in = N_evt-N_mu_out;

    cout << "particle total: " << N_prt << endl; 
    cout << "nbr of muon inside: " << N_mu_in << " (" << 100.*(N_mu_in)/N_evt << "%)" << endl;
    cout << "number of fails per condition:" << endl;
    size_t rem_prt=N_prt;
    for (int i=0; i<n_c.size(); i++) {
        cout << "c" << i << ": " << n_c[i] << "\r";
        if (i) cout << "\t\t(" << 100.*n_c[i]/rem_prt << "%)" << "\r";
        rem_prt -= n_c[i];
        cout << "\t\t\t\t rem: " << rem_prt << "\r";
        cout << "\t\t\t\t\t\t(" << 100.*rem_prt/N_prt << "%)" << "\r";
        cout << "\t\t\t\t\t\t\t(" << 100.*rem_prt/N_evt << "% of evt) (" << 100.*rem_prt/N_mu_in << "% of mu in)";
        cout << endl;
    }
}
