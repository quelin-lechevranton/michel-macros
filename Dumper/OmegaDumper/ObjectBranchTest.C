typedef struct {
	double X;
	double Y;
	double Z;
} Point;

void ObjectBranchTest() {
	TFile *f = new TFile("ObjectBranchTest.root","recreate");
	TTree *t = new TTree("Tree","a tree");

	Point pt={1,2,3};

	t->Branch("Branch",&pt,"X/D:Y/D:Z/D");
	t->Fill();
	t->Write();
}
