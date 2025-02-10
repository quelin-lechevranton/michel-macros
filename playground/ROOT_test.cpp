#include <TFile.h>
#include <TTree.h>
#include <ROOT/RVec.hxx>
#include <Math/GenVector/PositionVector3D.h>
#include <iostream>
#include <vector>


struct S {
    std::vector<Int_t> a;
    std::vector<Float_t> b;
};

namespace geo {
    using Point_t = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::GlobalCoordinateSystemTag>;
}

void ROOT_test() {
    TFile* f = new TFile("ROOT_test.root", "RECREATE");
    TTree* t= new TTree("t", "");

    // ROOT::RVec<S> v;
    S v;
    Int_t n;

    v.a.resize(20);
    v.b.resize(20);

    t->Branch("n", &n, "n/I");
    t->Branch("v", &v, "a[n]/I:b[n]/F");

    for (int i=0; i<10; i++) {
        n = 2*i;

        for (int j=0; j<n; j++) {
            v.a.push_back(j);
            v.b.push_back(j*1.1);
        }

        std::cout << i << " n:" << n << std::endl;

        t->Fill();
    }

    // std::vector<geo::Point_t> v;

    // t->Branch("v", &v);

    // for (int i=0; i<10; i++) {
    //     v.push_back(geo::Point_t(i, i, i));
    //     t->Fill();
    // }



    t->Write();
    t->Print();
    f->Close();
}