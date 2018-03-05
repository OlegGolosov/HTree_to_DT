#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>

void drawMass(char* input){

    TCanvas* canv;
    TLegend* leg;	
	TChain *ch = new TChain;
    char picName[200];
    char trees[200];
    sprintf(trees, "%s.root/tree", input);
    ch->Add(trees);

    TH2F *hMTheta = new TH2F("hMTheta", "Mass vs Theta; m [Mev/c^2]; #theta[rad]; nTracks", 1000, 0, 4000, 100, 0, 1.6);
    TH2F *hPhiTheta = new TH2F ("hPhiTheta", "Phi vs Theta; #phi [rad]; #theta [rad]; nTracks", 100, 0, 6.2, 100, 0, 1.6);

    Short_t nTracks;
    const Short_t maxNTracks = 200;
    Float_t pt[maxNTracks];
    Float_t theta[maxNTracks];
    Float_t mass[maxNTracks];
    Float_t phi[maxNTracks];

    ch->SetBranchAddress("pt", pt);
    ch->SetBranchAddress("mass", mass);
    ch->SetBranchAddress("theta", theta);
    ch->SetBranchAddress("nTracks", &nTracks);
    ch->SetBranchAddress("phi", phi);


    int entry=0;
    while(entry < ch->GetEntries()) {
        ch->GetEntry(entry);
        for (Int_t i = 0; i< nTracks; i++) {
            hMTheta->Fill(mass[i], theta[i]);
            hPhiTheta->Fill(phi[i], theta[i]);
        }
        entry++;
}


    canv = new TCanvas();
    leg = new TLegend(0.7,0.7,0.9,0.9);
    hMTheta->Draw("colz");
    hMTheta->SetStats(false);
    sprintf(picName, "./pics/MassVsTheta.png");
    canv->SaveAs(picName);
    sprintf(picName, "./pics/MassVsTheta.C");
    canv->SaveAs(picName);
    delete canv;
    delete leg;

    canv = new TCanvas();
    leg = new TLegend(0.7,0.7,0.9,0.9);
    hPhiTheta->Draw("colz");
    hPhiTheta->SetStats(false);
    sprintf(picName, "./pics/PhiVsTheta.png");
    canv->SaveAs(picName);
    sprintf(picName, "./pics/PhiVsTheta.C");
    canv->SaveAs(picName);
    delete canv;
    delete leg;
}
