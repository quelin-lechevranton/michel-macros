#include "OmegaLight_tools.h"

const size_t n_file=3;
const vector<string> filelist = omega::ReadFileList(n_file,"list/jeremy.list");

void Template() {
    clock_t start_time=clock();

    /* graph declarations */

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        omega::Truth T(filename.c_str());
        omega::Reco R(filename.c_str());

        for (size_t i_evt=0; i_evt < R.N; i_evt++) {

            cout << "Event#" << i_evt+1 << "/" << R.N << "\r" << flush;

            T.GetEvt(i_evt);
            
            for (size_t i_prt=0; i_prt < T.Prt.N; i_prt++) {
                
                T.GetPrt(i_prt);

                T.GetPrtDep(i_prt);

                int p = T.Prt.Pdg;


                for (size_t i_dep=0; i_dep < T.Dep.N; i_dep++) {
                    double e = *T.Dep.E[i_dep];
                } //end deposit loop
            } //end mcparticle loop

            // R.GetEvtPfp(i_evt);

            // for (size_t i_pfp=0; i_pfp < R.Pfp.N; i_pfp++) {

            //     cout << "\t" << i_pfp+1 << "/" << R.Pfp.isTrk.size() << endl;

            //     if (!R.Pfp.isTrk[i_pfp]) continue;
            //     R.GetPfpTrk(i_pfp);

            //     R.GetPfpSpt(i_pfp);
                
            //     for (size_t i_spt=0; i_spt < R.Spt.N; i_spt++) {

            //     } //end spacepoint loop
            // } //end pfparticle loop
        } //end event loop
    } //end file loop
    cout << endl;

    // TCanvas* c1 = new TCanvas("c1","Template");
    // c1->cd();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}