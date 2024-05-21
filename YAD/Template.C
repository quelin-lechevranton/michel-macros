#include "YADtools.h"

const size_t n_file=1;
const vector<string> filelist = yad::ReadFileList(n_file,"ijclab.list");

void Template() {

    size_t i_file=0;
    for (string file : filelist) {
        
        cout << "\e[3mOpening file #" << ++i_file << "/" << n_file << ": " << filename << "\e[0m" << endl;

        yad::Truth T(filename.c_str());
        yad::Reco R(filename.c_str());

        size_t n_evt = R.GetEntris();
        for (size_t i_evt=0; i_evt < n_evt; i_evt++) {

            T.GetEntry(i_evt);
            R.GetEntry(i_evt);

            for (size_t i_prt=0; i_prt < T.NPrt; i_prt++) {

            } //end mcparticle loop

            for (size_t i_dep=0; i_dep < T.NDep; i_dep++) {

            } //end deposit loop

            for (size_t i_pfp=0; i_pfp < R.NPfp; i_pfp++) {

                for (size_t i_clu=0; i_clu < R.PfpNClu->at(i_pfp); i_clu++) {

                } //end cluster loop

                for (size_t i_spt=0; i_spt < R.PfpNSpt->at(i_pfp); i_spt++) {

                } //end spacepoint loop
            } //end pfp loop

            for (size_t i_trk=0; i_trk < R.NTrk; i_trk++) {

                for (size_t i_tpt=0; i_tpt < R.TrkNPt->at(i_trk); i_tpt++) {

                } //end trackpoint loop
            } //end track loop
 
            for (size_t i_shw=0; i_shw < R.NSHw; i_shw++) {

            } //end shower loop
        } //end event loop
    } //end file loop
}