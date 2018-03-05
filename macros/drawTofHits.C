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

void drawTofHits(char* input) {

TCanvas* canv;
TLegend* leg;
TChain *ch = new TChain;
char picName[200];
char trees[200];
sprintf(trees, "%s.root/tree", input);
ch->Add(trees);

const Short_t nTriggers = 4;
Bool_t trigInd[nTriggers];
Int_t nTofHitsCut;

ch->SetBranchAddress("trigInd", trigInd);
ch->SetBranchAddress("nTofHitsCut", &nTofHitsCut);

TH1F* hTOFPT2 = new TH1F("hTof", "Tof Hits with PT; nEvents; nHits", 120, 0, 120);
TH1F* hTOFPT3 = new TH1F("hTof", "Tof Hits with PT; nEvents; nHits", 120, 0, 120);



int entry=0;
int count;
while(entry < ch->GetEntries()) {
    ch->GetEntry(entry);
    count = 0;
    if (trigInd[1]) hTOFPT2->Fill(nTofHitsCut);
	if (trigInd[2]) hTOFPT3->Fill(nTofHitsCut);	 
    entry++;
}


canv = new TCanvas();
hTOFPT2->SetLineColor(2);
hTOFPT3->SetLineColor(3);
hTOFPT3->SetStats(0);
hTOFPT2->SetStats(0);
hTOFPT2->Draw();
hTOFPT3->Draw("same");



leg = new TLegend(0.7,0.7,0.9,0.9);
leg->AddEntry(hTOFPT2, "PT2", "l");
leg->AddEntry(hTOFPT3, "PT3", "l");
leg->Draw("same");

sprintf(picName, "../results/TOF.png");
canv->SaveAs(picName);
sprintf(picName, "../results/TOF.C");
canv->SaveAs(picName);
delete canv;
delete leg;



}


