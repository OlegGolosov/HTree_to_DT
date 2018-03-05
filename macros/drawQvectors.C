//@brief Scripts for drawing trig variables
//@author T.Mamontova
//@date 20.11.17
//@param tree without .root
//@return pics with histos

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>

void drawQnvectors(char* input) {

TCanvas* canv;
TLegend* leg;
TChain *ch = new TChain;
char picName[200];
char trees[200];
sprintf(trees, "%s.root/tree", input);
ch->Add(trees);

const Short_t maxNTracks = 200;
Float_t phi[maxNTracks];
Short_t nTracks;
Float_t Q1_x;
Float_t Q1_y;

ch->SetBranchAddress("phi", phi);
ch->SetBranchAddress("nTracks", &nTracks);

TH2F* hQ1 = new TH2F("hQ1", "Q1_x vs Q1_y ; Q1_x; Q1_y; nEvents", 100, -10, 10, 100, -10, 10);



int entry=0;
while(entry < ch->GetEntries()) {
	Q1_x=0;
	Q1_y=0;
    ch->GetEntry(entry);
	for (Int_t i = 0; i< nTracks; i++) {
		Q1_x+=TMath::Cos(phi[i]);
		Q1_y+=TMath::Sin(phi[i]);
	}
	hQ1->Fill(Q1_x, Q1_y);
    entry++;
}


canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hQ1->Draw();
sprintf(picName, "../results/Q1.png");
canv->SaveAs(picName);
sprintf(picName, "../results/Q1.C");
canv->SaveAs(picName);
delete canv;
delete leg;



}


