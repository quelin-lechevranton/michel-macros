#include "YAD_tools.h"


#define FILLSPECTRUM(i) totalE=0;for (size_t i_dep=0; i_dep < n_el_dep ; i_dep+=1) {totalE+=(*T.DepE)[i_prt][i_dep];}hElE[i]->Fill(totalE);

// analysis parameter
vector<vector<double>> detlim={{-350, 350}, {-350,350}, {0,200}}; //cm
double coincidence_radius_gros=1.; //cm
size_t n_step_gros=10; //distance per step : 0.03 cm
double coincidence_radius_fin=1.; //cm
size_t n_scan_fin=4*n_step_gros;
size_t n_least_deposits=20;
size_t n_last_deposits=20;
double threshold_dEdx=3.5;


const size_t n_file=30;
const vector<string> filelist = yad::readFileList(n_file,"list/jeremy.list");


void fillElectronSpectrum(size_t,size_t);

void TrueMichel_v2() {

    clock_t start_time=clock();

    TColor color;

    // TGraph2D *gMu = new TGraph2D();
    // size_t igMu = 0;
    // // gMu->SetMarkerColor(kOrange+7);
    // gMu->SetMarkerColor(color.GetColor("#f35c22"));
    // gMu->SetMarkerStyle(20);
    // gMu->SetMarkerSize(0.3);

    // TGraph2D *gMux = new TGraph2D();
    // size_t igMux = 0;
    // gMux->SetMarkerColor(kRed);
    // gMux->SetMarkerStyle(20);
    // gMux->SetMarkerSize(0.3);
    
    // TGraph2D *gPrt = new TGraph2D();
    // size_t igPrt = 0;
    // gPrt->SetMarkerColor(kPink+7);
    // gPrt->SetMarkerStyle(20);
    // gPrt->SetMarkerSize(0.3);

    // TGraph2D *gEl = new TGraph2D();
    // size_t igEl = 0;
    // // gEl->SetMarkerColor(kAzure+2);
    // gEl->SetMarkerColor(color.GetColor("#436188"));
    // gEl->SetMarkerStyle(20);
    // gEl->SetMarkerSize(0.3);

    TH1D* hNDep = new TH1D("hNDep","Electron Deposit Number;nDep;#",100,0,1000);

    vector<TH2D*> hdEdx(2);
    hdEdx[0] = new TH2D("hdEdx0","Muon Energy Loss (before selection);Residual Range (cm);dEdx (MeV/cm)",100,0,200,50,0,5);
    hdEdx[1] = new TH2D("hdEdx1","Muon Energy Loss;Residual Range (cm);dEdx (MeV/cm)",100,0,200,50,0,5);

    vector<TH1D*> hAvgdEdx(2);
    hAvgdEdx[0] = new TH1D("hAvgdEdx0","Muon Energy Loss (before selection);dEdx (MeV/cm);#",30,0,30);
    hAvgdEdx[1] = new TH1D("hAvgdEdx1","Muon Energy Loss;dEdx (MeV/cm);#",30,0,30);

    size_t nElE=6;
    vector<TH1D*> hElE(nElE);
    hElE[0] = new TH1D("hElE0","Michel Spectrum;Total Deposited Energy (MeV);#",35,0,70);
    hElE[1] = new TH1D("hElE1","Electron Spectrum;Total Deposited Energy (MeV);#",35,0,70);
    hElE[2] = new TH1D("hElE2","NoFewDep Electron Spectrum;Total Deposited Energy (MeV);#",35,0,70);
    hElE[3] = new TH1D("hElE3","FromMu Electron Spectrum;Total Deposited Energy (MeV);#",35,0,70);
    hElE[4] = new TH1D("hElE4","NoFewDepMu Electron Spectrum;Total Deposited Energy (MeV);#",35,0,70);
    hElE[5] = new TH1D("hElE5","Inside Electron Spectrum;Total Deposited Energy (MeV);#",35,0,70);
    // THStack* hsElE= new THStack("hsElE",";Total Deposited Energy (MeV);#");
    for (size_t iElE=0; iElE < nElE; iElE++) {
        hElE[iElE]->SetLineWidth(2);
        hElE[iElE]->SetLineColor(color.GetColor("#436188"));
        // hElE[iElE]->SetLineColor(color.GetColor(iElE+2));
        // hsElE->Add(hElE[iElE]);
    }


    

    size_t N_evt=0;
    size_t N_mu_inside=0;
    size_t  N_prt=0,
            N_not_el=0,
            N_few_el_NDep=0,
            // N_orph=0,
            N_not_from_mu=0,
            N_few_mu_NDep=0,
            N_outside=0,
            N_no_bragg=0,
            N_mich=0;
    

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        // yad::Reco R(filename.c_str());


        size_t n_evt = T.GetEntries();
        N_evt+=n_evt;
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            T.GetEntry(i_evt);
            // R.GetEntry(i_evt);

            N_prt+=T.NPrt;

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {
                if (T.PrtMomID->at(i_prt)!=-1) continue;
                if (T.PrtPdg->at(i_prt)!=13 && T.PrtPdg->at(i_prt)!=-13) continue;
                if (yad::isInside(
                    T.DepX->at(i_prt),
                    T.DepY->at(i_prt),
                    T.DepZ->at(i_prt),
                    detlim[0][0], detlim[0][1],
                    detlim[1][0], detlim[1][1],
                    detlim[2][0], detlim[2][1]
                )) {
                    N_mu_inside++;
                }
            }

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {


                size_t n_el_dep = T.PrtNDep->at(i_prt);
                double totalE;
                if (T.PrtPdg->at(i_prt)!=11 && T.PrtPdg->at(i_prt)!=-11) {N_not_el++; continue;} //electrons only
                FILLSPECTRUM(1)

                if (n_el_dep < n_least_deposits) {N_few_el_NDep++; continue;} //enough electron deposits
                FILLSPECTRUM(2)

                int i_mom = T.PrtMomID->at(i_prt);
                // if (i_mom==-1) {N_orph++; continue;} //no orphelin electrons
                if (T.PrtPdg->at(i_mom)!=13) {N_not_from_mu++; continue;} //electrons coming from muons only
                FILLSPECTRUM(3)

                size_t n_mu_dep = T.PrtNDep->at(i_mom);
                if (n_mu_dep < n_last_deposits) {N_few_mu_NDep++; continue;} //enough deposits to check Bragg peak
                FILLSPECTRUM(4)


                if (!yad::isInside(
                    T.DepX->at(i_prt),
                    T.DepY->at(i_prt),
                    T.DepZ->at(i_prt),
                    detlim[0][0], detlim[0][1],
                    detlim[1][0], detlim[1][1],
                    detlim[2][0], detlim[2][1]
                )) {N_outside++; continue;}
                FILLSPECTRUM(5)


                //now we want to find the closest muon deposit to the electron track
                //first step is rough-and-ready second is more precise our the result of first step
                
                size_t n_max_gros=0, i_max_gros=0;

                //computing the number of electron deposit near every deposit of the muon mother
                for (size_t i_mu_dep=0; i_mu_dep < n_mu_dep; i_mu_dep+=n_step_gros) {

                    double x_mu = (*T.DepX)[i_mom][i_mu_dep];
                    double y_mu = (*T.DepY)[i_mom][i_mu_dep];
                    double z_mu = (*T.DepZ)[i_mom][i_mu_dep];

                    size_t n_coincidence=0;
                    for (size_t i_el_dep=0; i_el_dep < n_el_dep; i_el_dep+=n_step_gros) {
                        
                        double x_el = (*T.DepX)[i_prt][i_el_dep];
                        double y_el = (*T.DepY)[i_prt][i_el_dep];
                        double z_el = (*T.DepZ)[i_prt][i_el_dep];

                        double dist= yad::distance(x_el,y_el,z_el,x_mu,y_mu,z_mu);
                        if (dist < coincidence_radius_gros) {
                            n_coincidence++;
                        }
                    }

                    if ( n_coincidence > n_max_gros ) {
                        n_max_gros = n_coincidence;
                        i_max_gros = i_mu_dep;
                    }
                } //end deposit loop

                // more precise count
                size_t i_mu_dep_ini = i_max_gros-n_scan_fin > 0 ? i_max_gros-n_scan_fin : 0;
                size_t i_mu_dep_fin = i_max_gros+n_scan_fin > n_mu_dep ? n_mu_dep : i_max_gros+n_scan_fin;
                size_t n_max=0, i_max=0;
                for (size_t i_mu_dep = i_mu_dep_ini; i_mu_dep < i_mu_dep_fin; i_mu_dep++) {
                    
                    double x_mu = (*T.DepX)[i_mom][i_mu_dep];
                    double y_mu = (*T.DepY)[i_mom][i_mu_dep];
                    double z_mu = (*T.DepZ)[i_mom][i_mu_dep];

                    size_t n_coincidence=0;
                    for (size_t i_el_dep=0; i_el_dep < n_el_dep; i_el_dep++) {
                        
                        double x_el = (*T.DepX)[i_prt][i_el_dep];
                        double y_el = (*T.DepY)[i_prt][i_el_dep];
                        double z_el = (*T.DepZ)[i_prt][i_el_dep];

                        double dist= yad::distance(x_el,y_el,z_el,x_mu,y_mu,z_mu);
                        if (dist < coincidence_radius_fin) {
                            n_coincidence++;
                        }
                    }
                    if ( n_coincidence > n_max ) {
                        n_max = n_coincidence;
                        i_max = i_mu_dep;
                    }
                } //end deposit loop
                // cout << "evt#" << i_evt << ":prt#" << i_prt << ":i_max=" << i_max << endl;

                //compute total range and distances between points of muon track until electron emission
                double mu_range=0;
                vector<double> distances;
                double X0 = (*T.DepX)[i_mom][0];
                double Y0 = (*T.DepY)[i_mom][0];
                double Z0 = (*T.DepZ)[i_mom][0];
                for (size_t i_dep=0; i_dep <= i_max ; i_dep++) {

                    double X = (*T.DepX)[i_mom][i_dep];
                    double Y = (*T.DepY)[i_mom][i_dep];
                    double Z = (*T.DepZ)[i_mom][i_dep];

                    // gMu->SetPoint(igMu++,Y,Z,X);

                    double dist = yad::distance(X,Y,Z,X0,Y0,Z0);
                    distances.push_back(dist);
                    mu_range+=dist;
                    X0=X,Y0=Y,Z0=Z;
                }

                //check if muon is dying before emission
                double avg_last_dEdx=0;
                for (size_t i_dep=i_max-n_last_deposits; i_dep <= i_max ; i_dep++) {
                    double E = (*T.DepE)[i_mom][i_dep];
                    double dist=distances[i_dep];

                    avg_last_dEdx += E/dist/n_last_deposits;
                }
                // if (n_el_dep < 300) {
                //     cout << avg_last_dEdx << endl;
                // } else {
                //     cout << "\t" << avg_last_dEdx << endl;
                // }

                double mu_res_range = mu_range;
                for (size_t i_dep=1; i_dep <= i_max ; i_dep++) {

                    double E = (*T.DepE)[i_mom][i_dep];

                    double dist=distances[i_dep];
                    mu_res_range-=dist;
                    hdEdx[0]->Fill(mu_res_range, E/dist);
                }

                hAvgdEdx[0]->Fill(avg_last_dEdx);
                if (avg_last_dEdx<threshold_dEdx) {N_no_bragg++; continue;}
                hAvgdEdx[1]->Fill(avg_last_dEdx);

                N_mich++;

                //plot dEdx vs. ResidualRange until electron emission
                mu_res_range = mu_range;
                for (size_t i_dep=1; i_dep <= i_max ; i_dep++) {

                    double E = (*T.DepE)[i_mom][i_dep];

                    double dist=distances[i_dep];
                    mu_res_range-=dist;
                    hdEdx[1]->Fill(mu_res_range, E/dist);
                }

                // //plot all muon deposits before electron emission
                // for (size_t i_dep=0; i_dep <= i_max ; i_dep+=1) {

                //     double X = (*T.DepX)[i_mom][i_dep];
                //     double Y = (*T.DepY)[i_mom][i_dep];
                //     double Z = (*T.DepZ)[i_mom][i_dep];

                //     gMu->SetPoint(igMu++,Y,Z,X);
                //     // if (i_dep <= i_max-n_last_deposits) {
                //     //     gMu->SetPoint(igMu++,Y,Z,X);
                //     // } else {
                //     //     gMux->SetPoint(igMux++,Y,Z,X);
                //     // }
                // }

                //plot all prt points
                // for (size_t i_ppt=0; i_ppt < T.PrtNPt->at(i_mom); i_ppt++) {
                //     double X = (*T.PrtX)[i_mom][i_ppt];
                //     double Y = (*T.PrtY)[i_mom][i_ppt];
                //     double Z = (*T.PrtZ)[i_mom][i_ppt];

                //     gPrt->SetPoint(igPrt++,Y,Z,X);
                // } //end particlepoint loop

                //plot all electron deposits + total deposited energy
                // for (size_t i_dep=0; i_dep < n_el_dep ; i_dep+=1) {

                //     double X = (*T.DepX)[i_prt][i_dep];
                //     double Y = (*T.DepY)[i_prt][i_dep];
                //     double Z = (*T.DepZ)[i_prt][i_dep];

                //     // gEl->SetPoint(igEl++,Y,Z,X);
                // }

                FILLSPECTRUM(0);
                hNDep->Fill(n_el_dep);
            } //end particle loop
        } //end event loop
    } //end file loop

   

    TCanvas* c1 = new TCanvas("c1","TrueMichel_v2");
    c1->cd();
    // gMu->Draw("p");
    // gEl->Draw("samep");
    // // gPrt->Draw("samep");
    // // gMux->Draw("samep");

    TCanvas* c2 = new TCanvas("c2","TrueMichel_v2");
    // c2->Divide(2,2);
    // c2->cd(1);
    // hNDep->Draw("hist");
    // c2->cd(2);
    // hdEdx->Draw("colZ");
    // c2->cd(3);
    c2->Divide(3,2);
    for (size_t i=0; i<6; i++) {
        c2->cd(i+1);
        hElE[i]->Draw("hist");
    }
    // c2->cd();
    // hsElE->Draw();

    TCanvas* c3 = new TCanvas("c3","TrueMichel_v2");
    c3->Divide(2,1);
    c3->cd(1);
    hdEdx[0]->Draw("colZ");
    c3->cd(2);
    hdEdx[1]->Draw("colZ");

    TCanvas* c4 = new TCanvas("c4","TrueMichel_v2");
    c4->Divide(2,1);
    c4->cd(1);
    hAvgdEdx[0]->Draw("hist");
    c4->cd(2);
    hAvgdEdx[1]->Draw("hist");

    
    
    size_t rem_prt=N_prt;
    cout << "N_prt: " << N_prt << endl;

    cout << "\tN_not_el: " << N_not_el << " - " << 100.*N_not_el/rem_prt << "%" << endl;
    rem_prt-=N_not_el;
    cout << "\t\tremaining particles: " << rem_prt << endl;

    // cout << "\tN_orph: " << N_orph << " - " << 100.*N_orph/rem_prt << "%" << endl;
    // rem_prt-=N_orph;
    // cout << "\t\tremaining particles: " << rem_prt << endl;

    cout << "\tN_not_from_mu: " << N_not_from_mu << " - " << 100.*N_not_from_mu/rem_prt << "%" << endl;
    rem_prt-=N_not_from_mu;
    cout << "\t\tremaining particles: " << rem_prt << endl;

    cout << "\tN_few_el_NDep: " << N_few_el_NDep << " - " << 100.*N_few_el_NDep/rem_prt << "%" << endl;
    rem_prt-=N_few_el_NDep;
    cout << "\t\tremaining particles: " << rem_prt << endl;

    cout << "\tN_few_mu_NDep: " << N_few_mu_NDep << " - " << 100.*N_few_mu_NDep/rem_prt << "%" << endl;
    rem_prt-=N_few_mu_NDep;
    cout << "\t\tremaining particles: " << rem_prt << endl;

    cout << "\tN_outside: " << N_outside << " - " << 100.*N_outside/rem_prt << "%" << endl;
    rem_prt-=N_outside;
    cout << "\t\tremaining particles: " << rem_prt << endl;

    cout << "\tN_no_bragg: " << N_no_bragg << " - " << 100.*N_no_bragg/rem_prt << "%" << endl;
    rem_prt-=N_no_bragg;
    cout << "\t\tremaining particles: " << rem_prt << endl;

    cout << "N_mu_inside: " << N_mu_inside << " - " << 100.*N_mu_inside/N_evt << "% of events" << endl;
    cout << "N_mich: " << N_mich << " - " << 100.*N_mich/N_prt << "% of MCParticles or " << 100.*N_mich/N_evt << "% of events or " << 100.*N_mich/N_mu_inside << "% of inside muons" << endl;

    // c1->SaveAs("out/TrueMichel_v2Graph2D.root");
    // c2->SaveAs("out/TrueMichel_v2Hist.root");

    cout << N_evt << " events treated in " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}