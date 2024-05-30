#include "YAD_tools.h"

const size_t n_file=1;
const vector<string> filelist = yad::readFileList(n_file,"list/ijclab.list");

void Template() {
    clock_t start_time=clock();

    /* graph declarations */

    size_t i_file=0;
    for (string filename : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntries();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            cout << "Event#" << i_evt+1 << "/" << n_evt << "\r" << flush;

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {

                for (size_t i_ppt=0; i_ppt < T.PrtNPt->at(i_prt); i_ppt++) {

                } //end particlepoint loop

                for (size_t i_dep=0; i_dep < T.PrtNDep->at(i_prt); i_dep++) {

                } //end deposit loop
            } //end particle loop

            for (size_t i_pfp=0; i_pfp < R.NPfp; i_pfp++) {

                for (size_t i_clu=0; i_clu < R.PfpNClu->at(i_pfp); i_clu++) {

                } //end cluster loop

                for (size_t i_spt=0; i_spt < R.PfpNSpt->at(i_pfp); i_spt++) {

                } //end spacepoint loop
            } //end pfp loop

            for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {

                for (size_t i_tpt=0; i_tpt < R.TrkNPt->at(i_trk); i_tpt++) {

                } //end trackpoint loop

                for (size_t i_cal=0; i_cal < R.TrkCalNPt->at(i_trk); i_cal++) {

                } //end calorimetry loop
            } //end track loop
 
            for (size_t i_shw=0; i_shw < R.NShw; i_shw++) {

            } //end shower loop
        } //end event loop
    } //end file loop
    cout << endl;

    TCanvas* c1 = new TCanvas("c1","Template");
    c1->cd();

    cout << "total time of execution: " << static_cast<double>(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;
}