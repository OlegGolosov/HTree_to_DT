//@brief Scripts for drawing trig variables
//@author T.Mamontova
//@date 20.11.17
//@param tree without .root
//@return pics with histos

#include "TTree.h"
#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include <TNtuple.h>
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"

void drawHistoForTrigger(char* input) {

TChain *ch = new TChain;
char picName[200];
char trees[200];
sprintf(trees, "%s.root/tree", input);
ch->Add(trees);

const Short_t nTriggers = 4;
Bool_t trigInd[nTriggers];

ch->SetBranchAddress("trigInd", trigInd);

TH1F* hTrig = new TH1F("hCuts", "Event cuts; ; Events", nTriggers+1, 0, nTriggers+1);


int entry=0;
int count;
while(entry < ch->GetEntries()) {
    ch->GetEntry(entry);
    count = 0;
    for (Int_t i = 0; i < nTriggers; i++){
        if(trigInd[i]){
	 hTrig->Fill(i+0.5);
	 count++;}
		}
 
       if(count == 0) hTrig->Fill(nTriggers+0.5);
	
    
    entry++;
}

string cutName[nTriggers+1] = {"PT1","PT2","PT3","PT4", "No one"};

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
for (Int_t i=0; i < nTriggers+1; i++) { 
    hTrig->GetXaxis()->SetBinLabel(i+1,cutName[i].c_str());
}
hTrig->SetFillColor(38);
hTrig->SetBarWidth(0.9);
hTrig->SetStats(0);
hTrig->Draw("b");
sprintf(picName, "../results/Triggers.png");
canv->SaveAs(picName);
sprintf(picName, "../results/Triggers.C");
canv->SaveAs(picName);
delete canv;
delete leg;



}


