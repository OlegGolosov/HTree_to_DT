//@brief Scripts for drawing correlation variables
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

void DrawCorrelationHits(char* input) {

TCanvas* canv;
TLegend* leg;
TChain *ch = new TChain;
char picName[200];
char trees[200];
sprintf(trees, "%s/tree", input);
ch->Add(trees);

const Short_t nCuts = 8;
Bool_t cuts[nCuts];
Int_t nRpcClustCut;                      
Int_t nTofHitsCut;
Short_t nWallHitsTot;
Float_t wallChargeTot;
Short_t nTracks;
Float_t vX;      
Float_t vY;      
Float_t vZ;      
Float_t vChi2;


ch->SetBranchAddress("cuts", cuts);
ch->SetBranchAddress("nRpcClustCut", &nRpcClustCut);
ch->SetBranchAddress("nTofHitsCut", &nTofHitsCut);
ch->SetBranchAddress("nWallHitsTot", &nWallHitsTot);
ch->SetBranchAddress("wallChargeTot", &wallChargeTot);
ch->SetBranchAddress("nTracks", &nTracks);
ch->SetBranchAddress("vX", &vX);
ch->SetBranchAddress("vY", &vY);
ch->SetBranchAddress("vZ", &vZ);
ch->SetBranchAddress("vChi2", &vChi2);

TH2F* hTOFFWhits = new TH2F("hTOFFWhits", "TOF hits vs FW Hits; Tof Hits; FW hits; nEntries", 100, 0, 100, 120, 0, 120);
TH2F* hTOFFWhitsCut = new TH2F("hTOFFWhitsCut", "TOF hits vs FW Hits with cuts; TOF Hits; FW hits; nEntries", 100, 0, 100, 120, 0, 120);
TH2F* hTOFFWcharge = new TH2F("hTOFFWcharge", "TOF hits vs FW charge; TOF hits; FW charge; nEntries", 100, 0, 100, 500, 0, 500);
TH2F* hTOFFWchargeCut = new TH2F("hTOFFWchargeCut", "TOF hits vs FW charge with cuts; TOF hits; FW charge; nEntries", 100, 0, 100, 500, 0, 500);
TH2F* hTOFMDC = new TH2F("hTOFMDC", "TOF hits vs MDC tracks; TOF hits; MDC tracks; nEntries", 100, 0, 100, 120, 0, 120);
TH2F* hTOFMDCCut = new TH2F("hTOFMDCCut", "TOF hits vs MDC tracks with cuts; TOF hits; MDC tracks; nEntries", 100, 0, 100, 120, 0, 120);
TH2F* hTOFRPC = new TH2F("hTOFRPC", "TOF hits vs RPC hits; TOF hits; RPC hits; nEntries", 100, 0, 100, 150, 0, 150);
TH2F* hTOFRPCCut = new TH2F("hTOFRPCCut", "TOF hits vs RPC hits with cuts; TOF hits; RPC hits; nEntries", 100, 0, 100, 150, 0, 150);

TH2F* hFWhitscharge = new TH2F("hFWhitscharge", "FW hits vs FW charge; FW hits; FW charge; nEntries ", 120, 0, 120, 500, 0, 500);
TH2F* hFWhitschargeCut = new TH2F("hFWhitschargeCut", "FW hits vs FW charge cuts; FW hits; FW charge; nEntries", 120, 0, 120, 500, 0, 500);
TH2F* hFWhitsMDC = new TH2F("hFWhitsMDC", "FW hits vs MDC tracks; FW hits; MDC tracks; nEntries", 120, 0, 120, 120, 0, 120);
TH2F* hFWhitsMDCCut = new TH2F("hFWhitsMDCCut", "FW hits vs MDC tracks with cuts; FW hits; MDC tracks; nEntries", 120, 0, 120, 120, 0, 120);
TH2F* hFWhitsRPC = new TH2F("hFWhitsRPC", "FW hits vs RPC hits; FW hits; RPC hits; nEntries", 120, 0, 120, 150, 0, 150);
TH2F* hFWhitsRPCCut = new TH2F("hFWhitsRPCCut", "FW hits vs RPC hits with cuts; FW hits; RPC hits; nEntries", 120, 0, 120, 150, 0, 150);

TH2F* hFWchargeMDC = new TH2F("hFWchargeMDC", "FW charge vs MDC tracks; FW charge; MDC tracks; nEntries", 500, 0, 500, 120, 0, 120);
TH2F* hFWchargeMDCCut = new TH2F("hFWchargeMDCCut", "FW charge vs MDC tracks with cuts; FW charge; MDC tracks; nEntries", 500, 0, 500, 120, 0, 120);
TH2F* hFWchargeRPC = new TH2F("hFWchargeRPC", "FW charge vs RPC hits; FW charge; RPC hits; nEntries", 500, 0, 500, 150, 0, 150);
TH2F* hFWchargeRPCCut = new TH2F("hFWchargeRPCCut", "FW charge vs RPC hits with cuts; FW charge; RPC hits; nEntries", 500, 0, 500, 150, 0, 150);

TH2F* hMDCRPC = new TH2F("hMDCRPC", "MDC hits vs RPC hits; MDC hits; RPC hits; nEntries", 120, 0, 120, 150, 0, 150);
TH2F* hMDCRPCCut = new TH2F("hMDCRPCCut", "MDC hits vs RPC hits with cuts; MDC hits; RPC hits; nEntries", 120, 0, 120, 150, 0, 150);

TH1F* hTR = new TH1F("hTR", "TOF hits + RPC hits with cuts", 250, 0, 250);

Int_t entry =0;
// Bool_t check;
while (entry < ch->GetEntries()){
    ch->GetEntry(entry);
    hTOFFWhits->Fill(nTofHitsCut, nWallHitsTot);
    hTOFFWcharge->Fill(nTofHitsCut, wallChargeTot);
    hTOFMDC->Fill(nTofHitsCut, nTracks);
    hTOFRPC->Fill(nTofHitsCut, nRpcClustCut);
    hFWhitscharge->Fill(nWallHitsTot, wallChargeTot);
    hFWhitsMDC->Fill(nWallHitsTot, nTracks);
    hFWhitsRPC->Fill(nWallHitsTot, nRpcClustCut);
    hFWchargeMDC->Fill(wallChargeTot, nTracks);
    hFWchargeRPC->Fill(wallChargeTot, nRpcClustCut);
    hMDCRPC->Fill(nTracks, nRpcClustCut);

    if (cuts[0] && cuts[1] && cuts[2] && cuts[3] && cuts[4] && cuts[5] && cuts[6] && cuts[7] &&
        vZ<0 && vZ > -60 && 0.5<vChi2 && vChi2<40 && TMath::Sqrt(vX*vX+vY*vY)<3)
    {
        hTOFFWhitsCut->Fill(nTofHitsCut, nWallHitsTot);
        hTOFFWchargeCut->Fill(nTofHitsCut, wallChargeTot);
        hTOFMDCCut->Fill(nTofHitsCut, nTracks);
        hTOFRPCCut->Fill(nTofHitsCut, nRpcClustCut);
        hFWhitschargeCut->Fill(nWallHitsTot, wallChargeTot);
        hFWhitsMDCCut->Fill(nWallHitsTot, nTracks);
        hFWhitsRPCCut->Fill(nWallHitsTot, nRpcClustCut);
        hFWchargeMDCCut->Fill(wallChargeTot, nTracks);
        hFWchargeRPCCut->Fill(wallChargeTot, nRpcClustCut);
        hMDCRPCCut->Fill(nTracks, nRpcClustCut);
        hTR->Fill(nTofHitsCut+nRpcClustCut);
    }
    entry++;
    
}


canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hTOFFWhits->SetStats(kFALSE);
hTOFFWhits->Draw("colz");
sprintf(picName, "../results/TOFFWhits.png");
canv->SaveAs(picName);
sprintf(picName, "../results/TOFFWhits.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hTOFFWhitsCut->SetStats(kFALSE);
hTOFFWhitsCut->Draw("colz");
sprintf(picName, "../results/TOFFWhitsCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/TOFFWhitsCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hTOFFWcharge->SetStats(kFALSE);
hTOFFWcharge->Draw("colz");
sprintf(picName, "../results/TOFFWcharge.png");
canv->SaveAs(picName);
sprintf(picName, "../results/TOFFWcharge.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hTOFFWchargeCut->SetStats(kFALSE);
hTOFFWchargeCut->Draw("colz");
sprintf(picName, "../results/TOFFWchargeCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/TOFFWchargeCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hTOFMDC->SetStats(kFALSE);
hTOFMDC->Draw("colz");
sprintf(picName, "../results/TOFMDC.png");
canv->SaveAs(picName);
sprintf(picName, "../results/TOFMDC.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hTOFMDCCut->SetStats(kFALSE);
hTOFMDCCut->Draw("colz");
sprintf(picName, "../results/TOFMDCCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/TOFMDCCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hTOFRPC->SetStats(kFALSE);
hTOFRPC->Draw("colz");
sprintf(picName, "../results/TOFRPC.png");
canv->SaveAs(picName);
sprintf(picName, "../results/TOFRPC.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hTOFRPCCut->SetStats(kFALSE);
hTOFRPCCut->Draw("colz");
sprintf(picName, "../results/TOFRPCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/TOFRPCCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWhitscharge->SetStats(kFALSE);
hFWhitscharge->Draw("colz");
sprintf(picName, "../results/hFWhitscharge.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWhitscharge.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWhitschargeCut->SetStats(kFALSE);
hFWhitschargeCut->Draw("colz");
sprintf(picName, "../results/hFWhitschargeCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWhitschargeCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWhitsMDC->SetStats(kFALSE);
hFWhitsMDC->Draw("colz");
sprintf(picName, "../results/hFWhitsMDC.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWhitsMDC.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWhitsMDCCut->SetStats(kFALSE);
hFWhitsMDCCut->Draw("colz");
sprintf(picName, "../results/hFWhitsMDCCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWhitsMDCCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWhitsRPC->SetStats(kFALSE);
hFWhitsRPC->Draw("colz");
sprintf(picName, "../results/hFWhitsRPC.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWhitsRPC.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWhitsRPCCut->SetStats(kFALSE);
hFWhitsRPCCut->Draw("colz");
sprintf(picName, "../results/hFWhitsRPCCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWhitsRPCCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWchargeMDC->SetStats(kFALSE);
hFWchargeMDC->Draw("colz");
sprintf(picName, "../results/hFWchargeMDC.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWchargeMDC.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWchargeMDCCut->SetStats(kFALSE);
hFWchargeMDCCut->Draw("colz");
sprintf(picName, "../results/hFWchargeMDCCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWchargeMDCCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWchargeRPC->SetStats(kFALSE);
hFWchargeRPC->Draw("colz");
sprintf(picName, "../results/hFWchargeRPC.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWchargeRPC.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hFWchargeRPCCut->SetStats(kFALSE);
hFWchargeRPCCut->Draw("colz");
sprintf(picName, "../results/hFWchargeRPCCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hFWchargeRPCCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hMDCRPC->SetStats(kFALSE);
hMDCRPC->Draw("colz");
sprintf(picName, "../results/hMDCRPC.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hMDCRPC.C");
canv->SaveAs(picName);
delete canv;
delete leg;

canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hMDCRPCCut->SetStats(kFALSE);
hMDCRPCCut->Draw("colz");
sprintf(picName, "../results/hMDCRPCCut.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hMDCRPCCut.C");
canv->SaveAs(picName);
delete canv;
delete leg;


canv = new TCanvas();
leg = new TLegend(0.7,0.7,0.9,0.9);
hTR->SetStats(kFALSE);
hTR->Draw();
sprintf(picName, "../results/hTR.png");
canv->SaveAs(picName);
sprintf(picName, "../results/hTRC.C");
canv->SaveAs(picName);
delete canv;
delete leg;
}


