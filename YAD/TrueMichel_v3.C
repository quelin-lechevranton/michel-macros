#include "YAD_tools.h"

vector<vector<double>> detlim={{-320, 350}, {-317,317}, {20,280}}; //cm

const size_t n_file=30;
const vector<string> filelist = yad::readFileList(n_file,"list/jeremy.list");


void TrueMichel_v3() {
    clock_t start_time=clock();

    TColor color;

    size_t  N_evt=0,
            N_prt=0;
    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());

        size_t n_evt = T.GetEntries();
        N_evt+=n_evt;
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            T.GetEntry(i_evt);

            // N_prt+=T.NPrt;

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


                if (!yad::isInside(T.DepX->at(i_prt),T.DepY->at(i_prt),T.DepZ->at(i_prt))) {N_outside++; continue;}
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

                if (avg_last_dEdx<threshold_dEdx) {N_no_bragg++; continue;}

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