//@brief Scripts for drawing track variables
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
#include <TMath.h>
#include <iostream>

void drawTrackvariables(char* input) {

TCanvas* canv;
TLegend* leg;
TChain *ch = new TChain;
char picName[200];
char trees[200];
sprintf(trees, "%s.root/tree", input);
ch->Add(trees);

const Short_t maxNTracks = 200;
Short_t nTracks;
Float_t p[maxNTracks];
Float_t DCAxy[maxNTracks];
Float_t DCAz[maxNTracks];
Float_t chi2all[maxNTracks];
Float_t chi2inner[maxNTracks];
Float_t chi2outer[maxNTracks];
Float_t metaMatchRadius[maxNTracks];
Float_t pCorr[maxNTracks];
Float_t metaDx[maxNTracks];
Float_t metaDy[maxNTracks];
Short_t mdcSecId[maxNTracks];
Float_t vZ;
Float_t P;
Float_t Delta;
Float_t metaMass[maxNTracks];
Float_t theta[maxNTracks];
Float_t pt[maxNTracks];
Float_t phi[maxNTracks];
Float_t px_sum;
Float_t py_sum;


ch->SetBranchAddress("nTracks", &nTracks);
ch->SetBranchAddress("metaMatchRadius", metaMatchRadius);
ch->SetBranchAddress("metaDx", metaDx);
ch->SetBranchAddress("metaDy", metaDy);
ch->SetBranchAddress("chi2all", chi2all);
ch->SetBranchAddress("chi2inner", chi2inner);
ch->SetBranchAddress("DCAxy", DCAxy);
ch->SetBranchAddress("DCAz", DCAz);
ch->SetBranchAddress("mdcSecId", mdcSecId);
ch->SetBranchAddress("pCorr", pCorr);
ch->SetBranchAddress("p", p);
ch->SetBranchAddress("vZ", &vZ);
ch->SetBranchAddress("metaMass", metaMass);
ch->SetBranchAddress("theta", theta);
ch->SetBranchAddress("pt", pt);
ch->SetBranchAddress("phi", phi);


TH1F* hMetaR = new TH1F("hMetaR", "Distance of the outer segment to the meta hit ; r[mm];nTracks", 1000, 0, 100);
TH1F* hMetaDx = new TH1F("hMetaDx", "Distance in X of the outer segment to the meta hit; x[mm]; nTracks", 1000, 0, 100);
TH1F* hMetaDy = new TH1F("hMetaDy", "Distance in Y of the outer segment to the meta hit; y[mm]; nTracks", 1000, 0, 100);
TH1F* hChi2all = new TH1F("hChi2all", "Tracking #chi^2; #chi^2; nTracks", 1000, 0, 100);
TH1F* hChi2in = new TH1F("hChi2in", "MDC inner segment fitter #chi^2; #chi^2; nTracks", 1000, 0, 100);
TH1F* hChi2out = new TH1F("hChi2out", "MDC outer segment fitter #chi^2; #chi^2; nTracks", 1000, 0, 100);
TH1F* hDCAxy = new TH1F("hDCAxy", "Distance of closest point to beamline; r[mm]; nTracks", 100, -50, 50);
TH1F* hDCAz = new TH1F("hDCAz", "Z coordinate of closest point to beamline - Zvertex[mm]; z-Vz[mm]; nTracks", 100, -50, 50);
TH1F* hSecId = new TH1F("hSecId", "Sector information from MDC; nSector; nTracks", 6, -0.5, 5.5);
TH2F* hPcorr = new TH2F("hPcorr", "Abs(p-p_{corr}) vs p; p [Mev/c]; #Delta p [MeV/c]; nTracks", 100, 0, 2000, 1000, 0, 100);
TH1F* hMassTOF = new TH1F("hMassTOF", "Mass for TOF; mass [MeV/c^2]; nTracks", 4000, 0, 4000);
TH1F* hMassRPC = new TH1F("hMassRPC", "Mass for RPC; mass [MeV/c^2]; nTracks", 4000, 0, 4000);
TH1F* hMass = new TH1F("hMass", "Mass histo; mass[MeV/c^2]; nTracks", 4000, 0, 4000);
TH2F* hPtPolar = new TH2F("hPtPolar", "Polar histo for Pt; p_x [Mev/c]; p_y[Mev/c]", 1000, -2000, 2000, 1000, -2000, 2000);
TH2F* hPtSum = new TH2F("hPtSum", "Polar histo for p_T sum; p_x[MeV/c]; p_y[MeV/c]", 1000, -1000, 1000, 1000, -1000, 1000);




int entry=0;
std::cout << "nEntries = "<< ch->GetEntries() << std::endl;
while(entry < ch->GetEntries()) {

    ch->GetEntry(entry);
    //std::cout << "nTracks = "<<nTracks << std::endl;
    px_sum = 0;
    py_sum = 0;
    for (Short_t i = 0; i < nTracks; i++){
        hMetaR->Fill(metaMatchRadius[i]);
        hMetaDx->Fill(metaDx[i]);
        hMetaDy->Fill(metaDy[i]);
        hChi2all->Fill(chi2all[i]);
        hChi2in->Fill(chi2inner[i]);
        hChi2out->Fill(chi2outer[i]);
        hDCAxy->Fill(DCAxy[i]);
        hDCAz->Fill(DCAz[i]-vZ);
        hSecId->Fill(mdcSecId[i]);
        P = p[i];
        Delta = TMath::Abs(p[i]-pCorr[i]);
        hPcorr->Fill(P, Delta);
        hMass->Fill(metaMass[i]);
        //std::cout << theta[i] <<std::endl;
        if (1.2 < theta[i]< 1.3) 
            {
                //std::cout << theta[i] << "" << metaMass[i] <<std::endl;
                hMassTOF->Fill(metaMass[i]);
            }

        if (0.4 < theta[i] < 0.5) hMassRPC->Fill(metaMass[i]);
        hPtPolar->Fill(pt[i]*TMath::Cos(phi[i]), pt[i]*TMath::Sin(phi[i]));
        px_sum =px_sum+pt[i]*TMath::Cos(phi[i]);
        py_sum =py_sum+pt[i]*TMath::Sin(phi[i]);

    }
    hPtSum->Fill(px_sum, py_sum);
    entry++;
}


canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hMetaR->Draw();
sprintf(picName, "../results/MetaR.png");
canv->SaveAs(picName);
sprintf(picName, "../results/MetaR.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hMetaDx->SetLineColor(2);
hMetaDy->SetLineColor(3);
hMetaDx->SetStats(0);
hMetaDy->SetStats(0);
hMetaDx->Draw();
hMetaDy->Draw("same");
leg->AddEntry(hMetaDx, "MetaDx", "lp");
leg->AddEntry(hMetaDy, "MetaDy", "lp");
leg->Draw("same");
sprintf(picName, "../results/MetaDxy.png");
canv->SaveAs(picName);
sprintf(picName, "../results/MetaDxy.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hChi2all->SetLineColor(2);
hChi2in->SetLineColor(3);
hChi2out->SetLineColor(4);
hChi2all->SetStats(0);
hChi2in->SetStats(0);
hChi2out->SetStats(0);
hChi2all->Draw();
hChi2in->Draw("same");
hChi2out->Draw("same");
leg->AddEntry(hChi2all, "Chi2all", "l");
leg->AddEntry(hChi2in, "Chi2inner", "l");
leg->AddEntry(hChi2out, "Chi2out", "l");
leg->Draw("same");

sprintf(picName, "../results/Chi2.png");
canv->SaveAs(picName);
sprintf(picName, "../results/Chi2.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hDCAxy->Draw();
sprintf(picName, "../results/DCAxy.png");
canv->SaveAs(picName);
sprintf(picName, "../results/DCAxy.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hDCAz->Draw();
sprintf(picName, "../results/DCAz.png");
canv->SaveAs(picName);
sprintf(picName, "../results/DCAz.C");
canv->SaveAs(picName);
delete canv;
delete leg;


canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hSecId->Draw();
sprintf(picName, "../results/SecId.png");
canv->SaveAs(picName);
sprintf(picName, "../results/SecID.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hPcorr->Draw("colz");
hPcorr->SetStats(false);
sprintf(picName, "../results/PCorr.png");
canv->SaveAs(picName);
sprintf(picName, "../results/PCorr.C");
canv->SaveAs(picName);
delete canv;
delete leg;


canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hMassRPC->SetLineColor(2);
hMassTOF->SetLineColor(3);
hMassRPC->SetStats(0);
hMassTOF->SetStats(0);
hMassRPC->Draw();
hMassTOF->Draw("same");
leg->AddEntry(hMassRPC, "Meta mass for RPC", "lp");
leg->AddEntry(hMassTOF, "Meta mass for TOF", "lp");
leg->Draw("same");
sprintf(picName, "../results/MassMeta.png");
canv->SaveAs(picName);
sprintf(picName, "../results/MassMeta.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hPtPolar->Draw("colz");
hPtPolar->SetStats(false);
sprintf(picName, "../results/Pt_polar.png");
canv->SaveAs(picName);
sprintf(picName, "../results/Pt_polar.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hPtSum->Draw("colz");
hPtSum->SetStats(false);
sprintf(picName, "../results/Pt_sum.png");
canv->SaveAs(picName);
sprintf(picName, "../results/Pt_sum.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hMetaR->Draw();
sprintf(picName, "../results/Mass.png");
canv->SaveAs(picName);
sprintf(picName, "../results/Mass.C");
canv->SaveAs(picName);
delete canv;
delete leg;


delete hMetaR;
delete hMetaDx;
delete hMetaDy;
delete hChi2all;
delete hChi2in;
delete hChi2out;
delete hDCAxy;
delete hDCAz;
delete hSecId;
delete hPcorr; 
delete hMassTOF;
delete hMassRPC;
delete hPtPolar;
delete hPtSum;
delete hMass;

}

