#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>

void drawVertex(char* input){

    TCanvas* canv;
    TLegend* leg;	
	TChain *ch = new TChain;
    char picName[200];
    char trees[200];
    sprintf(trees, "%s.root/tree", input);
    ch->Add(trees);

    TH1F *hvZ = new TH1F("hvZ", "Z coordinate for vertex; v_Z[mm]; nEvents", 1000, -100, 50);
    TH1F *hvZCut = new TH1F("hvZCut", "Z coordinate for vertex with cuts; v_Z[mm]; nEvents", 1000, -100, 50);
    TH2F *hvXvY = new TH2F ("hvXvY", "Vertex distribution in xy plane; v_X[mm]; v_Y[mm]; nEvents", 100, -10, 10, 100, -10, 10);
    TH2F *hvXvYCut = new TH2F ("hvXvY", "Vertex distribution in xy plane with cuts; v_X[mm]; v_Y[mm]; nEvents", 100, -10, 10, 100, -10, 10);
    TH1F *hvChi2 = new TH1F("hvChi2", "#chi^2 distribution; #chi^2 ; nEvents", 100, 0, 100);
    TH1F *hvChi2Cut = new TH1F("hvChi2Cut", "#chi^2 distribution with; #chi^2 ; nEvents", 100, 0, 100);

    const Short_t nCuts = 8;
    Bool_t cuts[nCuts];
    Float_t vX;      
    Float_t vY;      
    Float_t vZ;      
    Float_t vChi2;

    ch->SetBranchAddress("vX", &vX);
    ch->SetBranchAddress("vY", &vY);
    ch->SetBranchAddress("vZ", &vZ);
    ch->SetBranchAddress("vChi2", &vChi2);
    ch->SetBranchAddress("cuts", cuts);

    Int_t entry =0;

    while (entry < ch->GetEntries()){
        ch->GetEntry(entry);
        hvZ->Fill(vZ);
        hvXvY->Fill(vX, vY);
        hvChi2->Fill(vChi2);

        if (cuts[0] && cuts[1] && cuts[2] && cuts[3] && cuts[4] && cuts[5] && cuts[6] && cuts[7] &&
            vZ<0 && vZ > -60 && 0.5<vChi2 && vChi2<40 && TMath::Sqrt(vX*vX+vY*vY)<3)
        {
            hvZCut->Fill(vZ);
            hvXvYCut->Fill(vX, vY);
            hvChi2Cut->Fill(vChi2);
        }
        entry++;
    }

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hvZ->SetLineColor(2);
hvZCut->SetLineColor(3);
hvZ->SetStats(0);
hvZCut->SetStats(0);
hvZ->Draw();
hvZCut->Draw("same");
leg->AddEntry(hvZ, "without cuts", "l");
leg->AddEntry(hvZCut, "with cuts", "l");
leg->Draw("same");

sprintf(picName, "../results/vZ.png");
canv->SaveAs(picName);
sprintf(picName, "../results/vZ.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hvXvY->SetStats(kFALSE);
hvXvY->Draw("colz");
sprintf(picName, "../results/vXvY.png");
canv->SaveAs(picName);
sprintf(picName, "../results/vXvY.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hvXvYCut->SetStats(kFALSE);
hvXvYCut->Draw("colz");
sprintf(picName, "../results/vXvYCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/vXvYCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hvChi2->SetLineColor(2);
hvChi2Cut->SetLineColor(3);
hvChi2->SetStats(0);
hvChi2Cut->SetStats(0);
hvChi2->Draw();
hvChi2Cut->Draw("same");
leg->AddEntry(hvChi2, "without cuts", "l");
leg->AddEntry(hvChi2Cut, "with cuts", "l");
leg->Draw("same");
sprintf(picName, "../results/Chi2.png");
canv->SaveAs(picName);
sprintf(picName, "../results/Chi2.C");
canv->SaveAs(picName);
delete canv;
delete leg;

}
