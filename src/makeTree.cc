//Alexander Zaytsev, Elizaveta Zherebtsova, Tatiana Mamontova
//a.zaytsev@gsi.de
//October 2016, November 2017

#include "hades.h"
#include "hloop.h"
#include "htime.h"
#include "hcategory.h"
#include "hcategorymanager.h"
#include "hparticlecand.h"
#include "hparticletracksorter.h"
#include "hparticlebooker.h"
#include "hparticletool.h"
#include "hparticledef.h"
#include "hparticleevtinfo.h"
#include "henergylosscorrpar.h"
#include "hphysicsconstants.h"
#include "hwallhit.h"
#include "walldef.h"
#include "hruntimedb.h"
#include "hrun.h"
#include "heventheader.h"

#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TRandom2.h"

#include "mhwalldivider.h"

#include <iostream>
#include <vector>

using namespace std;


const Float_t D2R = TMath::DegToRad(); //deg to rad
const Float_t Y_BEAM = 0.740151;  //beam rapidity for 1.23GeV/c nucleon

//  infileList : comma seprated file list "file1.root,file2.root" or "something*.root"
//  outfile    : optional (not used here) , used to store hists in root file
//  nEvents    : number of events to processed. if  nEvents < entries or < 0 the chain will be processed

Int_t makeTree(TString infileList, TString outfile, Int_t nEvents=-1)
{
    Bool_t isSimulation = kFALSE;

    // create loop object and hades
    HLoop loop(kTRUE);

    // list of all files with working sectors
    if(!isSimulation) loop.readSectorFileList("/lustre/nyx/hades/dst/apr12/gen8/sector_selection/FileListHadron.list",kFALSE,kFALSE);

    // reading input files and declaring containers
    Bool_t ret = kFALSE;
    if(infileList.Contains(",")){
        ret = loop.addMultFiles(infileList); // file1,file2,file3
    } else{
        ret = loop.addFiles(infileList); // myroot*.root
    }

    if(ret == 0) {
        cerr<<"READBACK: ERROR : cannot find inputfiles : "<<infileList.Data()<<endl;
        return 1;
    }

    //if(!loop.setInput("")) {   // all input categories
    if(!loop.setInput("-*,+HParticleCand,+HParticleEvtInfo,+HWallHit")) {
        cerr<<"READBACK: ERROR : cannot read input !"<<endl;
        exit(1);
    }
    // read all categories
    loop.printCategories();
    loop.printChain();

    //parameters
    HEnergyLossCorrPar dEdxCorr;
    dEdxCorr.setDefaultPar("apr12");

    // input data
    HCategory* candCat    = (HCategory*)HCategoryManager::getCategory(catParticleCand);
    HCategory* evtInfoCat = (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo);
    HCategory* wallCat    = (HCategory*)HCategoryManager::getCategory(catWallHit);


    //#######################################################################
    //#######################################################################

    // ###### ######  #####  #####
    //   ##   ##  ##  ##     ##
    //   ##   #####   #####  #####
    //   ##   ## ##   ##     ##
    //   ##   ##  ##  #####  #####
    //-----------------------------------------------------------------------



    const Int_t nRings = 9;                  //number of rings in FW
    Short_t nWallHitsTot;                    //total number of FW hits
    const Short_t nCuts = 8;                 //number of cuts for event
    const Int_t nModules = 304;              //number of modules in FW
    const Int_t nTriggers = 4;               //triggers PT1-4


    //Hits
    Int_t nRpcClust;                         //number of RPC clusters
    Int_t nRpcClustCut;                      //number of RPC clusters with time cut
    Int_t nRpcHits;                          //number of RPC hits with time cut
    Int_t nRpcHitsCut;                       //number of RPC hits with time cut
    Int_t nTofHits;                          //number of TOF hits
    Int_t nTofHitsCut;                       //number of TOF hits with time cut
    Int_t primaryTracks;                     //number of primary tracks
    Int_t selectedTracks;                    //number of selected tracks
    Float_t psiEPa;                          //EPa angle from FW
    Float_t psiEPb;                          //EPb angle from FW
    Int_t nA;                                //n hits in sample A
    Int_t nB;                                //n hits in sample A
    Bool_t trigInd[nTriggers];               //type of trigger
    Short_t runId;                           //run number
    Float_t runTime;                         //run time

    //FW hits
    const Int_t maxNWallHits = 200;          //maximal number of FW hits
    Short_t wallModuleIndex[maxNWallHits];   //index of hit module
    Float_t wallHitTime[maxNWallHits];       //time of FW hit
    Float_t wallHitCharge[maxNWallHits];     //charge of FW hit
    Float_t wallHitDistance[maxNWallHits];   //distance to FW module (may be used together with hit time to reject secondary particles)
    Short_t wallHitRing[maxNWallHits];       //wall ring of a hit
    Float_t wallHitPhi[maxNWallHits];        //phi angle of hit module in rad
    Float_t wallHitX[maxNWallHits];          //hit module X
    Float_t wallHitY[maxNWallHits];          //hit module Y
    Float_t wallHitZ[maxNWallHits];          //hit module Z
    Float_t wallHitEta[maxNWallHits];        //eta of FW hit
    Bool_t  isWallHitOk[maxNWallHits];       //output of hit cuts
    Float_t wallChargeTot;                   //total charge FW
    Float_t wallChargeTot_ring[nRings];      //FW rings charges (1-9)
    Float_t wallChargeTot_mod[nModules];     //FW modules charges



    //MDC
    Short_t nTracks = 0;    //number of all tracks in MDC
    Short_t nProtons = 0;   //number of protons
    Short_t nPiP = 0;       //number of pi+
    Short_t nPiM = 0;       //number of pi-

    //Vertex
    Float_t vX;      //x coordinate of primary vertex
    Float_t vY;      //y coordinate of primary vertex
    Float_t vZ;      //z coordinate of primary vertex
    Float_t vChi2;   //vertex reconstruction chi^2

    //Time
    Int_t time;

    //Cuts
    Bool_t cuts[nCuts]; //trigger; vertexClust; vertexCand; goodSTART;
                        //noPileUpSTART; noVETO; goodSTARTVETO; goodSTARTMETA

    //tracks
    const Short_t maxNTracks = 200;
    Short_t pid[maxNTracks];
    Float_t p[maxNTracks];
    Float_t phi[maxNTracks];
    Float_t theta[maxNTracks];
    Float_t pt[maxNTracks];
    Float_t rapidity[maxNTracks];
    Float_t eta[maxNTracks];
    Float_t metaBeta[maxNTracks];
    Float_t metaMass[maxNTracks];
    Short_t charge[maxNTracks];
    Float_t mdcdEdx[maxNTracks];
    Float_t tofdEdx[maxNTracks];
    Float_t DCAxy[maxNTracks];
    Float_t DCAz[maxNTracks];
    Short_t mdcNhits[maxNTracks];
    Short_t mdcNhitsInner[maxNTracks];
    Short_t mdcNhitsOuter[maxNTracks];
    Float_t chi2all[maxNTracks];
    Float_t chi2inner[maxNTracks];
    Float_t chi2outer[maxNTracks];
    Float_t metaQ[maxNTracks];
    Float_t metaMatchRadius[maxNTracks];
    Float_t pCorr[maxNTracks];
    Float_t pt_corr[maxNTracks];
    Float_t rapidity_corr[maxNTracks];
    Float_t metaDx[maxNTracks];
    Float_t metaDy[maxNTracks];
    Short_t mdcSecId[maxNTracks];

    TFile* out = new TFile(outfile.Data(),"RECREATE");
    out->cd();

    TTree* tree = new TTree("tree", "HADES Au+Au 1.23 GeV 8gen tree for flow analysis");


    tree->Branch("nRpcClust",          &nRpcClust,         "nRpcClust/I");
    tree->Branch("nRpcClustCut",       &nRpcClustCut,      "nRpcClustCut/I");
    tree->Branch("nRpcHits",           &nRpcHits,          "nRpcHits/I");
    tree->Branch("nRpcHitsCut",        &nRpcHitsCut,       "nRpcHitsCut/I");
    tree->Branch("nTofHits",           &nTofHits,          "nTofHits/I");
    tree->Branch("nTofHitsCut",        &nTofHitsCut,       "nTofHitsCut/I");
    tree->Branch("primaryTracks",      &primaryTracks,     "primaryTracks/I");
    tree->Branch("selectedTracks",     &selectedTracks,    "selectedTracks/I");
    tree->Branch("psiEPa",             &psiEPa,            "psiEPa/F");
    tree->Branch("psiEPb",             &psiEPb,            "psiEPb/F");
    tree->Branch("nA",                 &nA,                "nA/I");
    tree->Branch("nB",                 &nB,                "nB/I");
    tree->Branch("trigInd",            &trigInd,           TString::Format("trigInd[%i]/O", nTriggers));
    tree->Branch("runId",              &runId,             "runId/S");
    tree->Branch("nWallHitsTot",       &nWallHitsTot,      "nWallHitsTot/S");
    tree->Branch("runTime",            &runTime,           "runTime/F");

    tree->Branch("cuts",             cuts,            TString::Format("cuts[%i]/O", nCuts));

    //FW Hits
    tree->Branch("wallModuleIndex", wallModuleIndex, "wallModuleIndex[nWallHitsTot]/S");
    tree->Branch("wallHitTime",     wallHitTime,     "wallHitTime[nWallHitsTot]/F");
    tree->Branch("wallHitCharge",   wallHitCharge,   "wallHitCharge[nWallHitsTot]/F");
    tree->Branch("wallHitDistance", wallHitDistance, "wallHitDistance[nWallHitsTot]/F");
    tree->Branch("wallHitRing",     wallHitRing,     "wallHitRing[nWallHitsTot]/S");
    tree->Branch("wallHitPhi",      wallHitPhi,      "wallHitPhi[nWallHitsTot]/F");
    tree->Branch("wallHitX",        wallHitX,        "wallHitX[nWallHitsTot]/F");
    tree->Branch("wallHitY",        wallHitY,        "wallHitY[nWallHitsTot]/F");
    tree->Branch("wallHitZ",        wallHitZ,        "wallHitZ[nWallHitsTot]/F");
    tree->Branch("wallHitEta",      wallHitEta,      "wallHitEta[nWallHitsTot]/F");
    tree->Branch("isWallHitOk",     isWallHitOk,     "isWallHitOk[nWallHitsTot]/O");
    tree->Branch("wallChargeTot_mod",   wallChargeTot_mod,      TString::Format("wallChargeTot_mod[%i]/F", nModules));
    tree->Branch("wallChargeTot",       &wallChargeTot,         "wallChargeTot/F");
    tree->Branch("wallChargeTot_ring",  wallChargeTot_ring,     TString::Format("wallChargeTot_ring[%i]/F", nRings));

    //MDC
    tree->Branch("nTracks",  &nTracks,  "nTracks/S");
    tree->Branch("nProtons", &nProtons, "nProtons/S");
    tree->Branch("nPiP",     &nPiP,     "nPiP/S");
    tree->Branch("nPiM",     &nPiM,     "nPiM/S");

    //Vertex
    tree->Branch("vX",    &vX,    "vX/F");
    tree->Branch("vY",    &vY,    "vY/F");
    tree->Branch("vZ",    &vZ,    "vZ/F");
    tree->Branch("vChi2", &vChi2, "vChi2/F");

    //Time
    tree->Branch("time", &time, "time/I");

    //Tracks
    tree->Branch("pid",      pid,      "pid[nTracks]/S");
    tree->Branch("p",        p,        "p[nTracks]/F");
    tree->Branch("phi",      phi,      "phi[nTracks]/F");
    tree->Branch("theta",    theta,    "theta[nTracks]/F");
    tree->Branch("pt",       pt,       "pt[nTracks]/F");
    tree->Branch("rapidity", rapidity, "rapidity[nTracks]/F");
    tree->Branch("eta",      eta,      "eta[nTracks]/F");
    tree->Branch("metaBeta", metaBeta, "metaBeta[nTracks]/F");
    tree->Branch("metaMass", metaMass, "metaMass[nTracks]/F");
    tree->Branch("charge",   charge,   "charge[nTracks]/S");
    tree->Branch("mdcdEdx",  mdcdEdx,  "mdcdEdx[nTracks]/F");
    tree->Branch("tofdEdx",  tofdEdx,  "tofdEdx[nTracks]/F");
    tree->Branch("DCAxy",    DCAxy,    "DCAxy[nTracks]/F");
    tree->Branch("DCAz",     DCAz,     "DCAz[nTracks]/F");
    tree->Branch("mdcNhits", mdcNhits, "mdcNhits[nTracks]/S");
    tree->Branch("mdcNhitsInner", mdcNhitsInner, "mdcNhitsInner[nTracks]/S");
    tree->Branch("mdcNhitsOuter", mdcNhitsOuter, "mdcNhitsOuter[nTracks]/S");
    tree->Branch("chi2all",  chi2all,  "chi2all[nTracks]/F");
    tree->Branch("chi2inner",chi2inner,"chi2inner[nTracks]/F");
    tree->Branch("chi2outer",chi2outer,"chi2outer[nTracks]/F");
    tree->Branch("metaQ",    metaQ,    "metaQ[nTracks]/F");
    tree->Branch("metaMatchRadius",  metaMatchRadius,  "metaMatchRadius[nTracks]/F");
    tree->Branch("pCorr",    pCorr,   "pCorr[nTracks]/F");
    tree->Branch("pt_corr",  pt_corr,  "pt_corr[nTracks]/F");
    tree->Branch("rapidity_corr", rapidity_corr,"rapidity_corr[nTracks]/F");
    tree->Branch("metaDx",   metaDx,    "metaDx[nTracks]/F");
    tree->Branch("metaDy",   metaDy,    "metaDy[nTracks]/F");
    tree->Branch("mdcSecId", mdcSecId,  "mdcSecId[nTracks]/F");


    //#######################################################################
    //#######################################################################

    //##      #######  #######  #####
    //##      ##   ##  ##   ##  ##  ##
    //##      ##   ##  ##   ##  #####
    //##      ##   ##  ##   ##  ##
    //######  #######  #######  ##
    //-----------------------------------------------------------------------

    Int_t entries = loop.getEntries();
    if(nEvents < entries && nEvents >= 0 ) entries = nEvents;
    TString filename;
    Int_t sectors [6];
    MHWallDivider* divider = new MHWallDivider();
    TRandom2 *randomGenerator = new TRandom2();
    Float_t rndm;

    for (Int_t i = 0; i < entries; i++) {
        Int_t nbytes =  loop.nextEvent(i);             // get next event. categories will be cleared before
        if(nbytes <= 0) { cout<<nbytes<<endl; break; } // last event reached
        if(i%5000 == 0) cout<<"event "<<i<<endl;

        loop.getSectors(sectors); // fill sector array

        Int_t g, day, hour, minute;
        TString* be = new TString("be");

        if(loop.isNewFile(filename)){
            if(!isSimulation) filename = HTime::stripFileName(filename,kTRUE,kFALSE);
            HTime::splitFileName(filename,*be,g,day,hour,minute,g,g,kFALSE);
            time = day*24*60 + hour*24 + minute;
        }

        //-------------------------------------------------
        // summary event info object
        HParticleEvtInfo* evtInfo=0;
        evtInfo = HCategoryManager::getObject(evtInfo,evtInfoCat,0 );

        //cut info
        cuts[0] = evtInfo->isGoodEvent(Particle::kGoodTRIGGER);
        cuts[1] = evtInfo->isGoodEvent(Particle::kGoodVertexClust);
        cuts[2] = evtInfo->isGoodEvent(Particle::kGoodVertexCand);
        cuts[3] = evtInfo->isGoodEvent(Particle::kGoodSTART);
        cuts[4] = evtInfo->isGoodEvent(Particle::kNoPileUpSTART);
        cuts[5] = evtInfo->isGoodEvent(Particle::kNoVETO);
        cuts[6] = evtInfo->isGoodEvent(Particle::kGoodSTARTVETO);
        cuts[7] = evtInfo->isGoodEvent(Particle::kGoodSTARTMETA);

        //get Run number
        runId = gHades->getCurrentEvent()->getHeader()->getEventRunNumber();
        runTime = gHades->getCurrentEvent()->getHeader()->getTime();

        //get type of trigger

        for (Int_t k = 0; k < nTriggers; k++){
            trigInd[k] = kFALSE;
            if (gHades->getCurrentEvent()->getHeader()->isTBit(k+11)) trigInd[k] = kTRUE;
        }

        //get primary vertex
        HVertex vertexReco = gHades->getCurrentEvent()->getHeader()->getVertexReco();

        vX = vertexReco.getX();
        vY = vertexReco.getY();
        vZ = vertexReco.getZ();
        vChi2 = vertexReco.getChi2();

       //Track Primary and Hits

        psiEPa = evtInfo -> getPhiA () * D2R;
        psiEPb = evtInfo -> getPhiB () * D2R;
        nA = evtInfo -> getNA ();
        nB = evtInfo -> getNB ();
        nRpcClust = evtInfo->getSumRpcMult();
        nRpcClustCut = evtInfo->getSumRpcMultCut();
        nRpcHits = evtInfo->getSumRpcMultHit();
        nRpcHitsCut = evtInfo->getSumRpcMultHitCut();
        nTofHits = evtInfo->getSumTofMult();
        nTofHitsCut = evtInfo->getSumTofMultCut();
        primaryTracks=evtInfo->getSumPrimaryParticleCandMult();
        selectedTracks=evtInfo->getSumSelectedParticleCandMult();

        // loop over FW hits
        HWallHit* wallHit = 0;

        nWallHitsTot = wallCat->getEntries();
        if (nWallHitsTot>maxNWallHits) {
            cerr << "ERROR: nWallHitsTot="<<nWallHitsTot<<" > maxNWallHits="<<maxNWallHits<<endl;
            return 1;
        }

        Short_t subevent;
        Short_t ring;
        Float_t psi;
        Short_t randomSubEvent;
        Float_t hit_beta;
        Float_t thetaFW;

        wallChargeTot=0;
        for (Short_t j=0; j<nRings; j++){
            wallChargeTot_ring[j] = 0;
        }
        for (Short_t j=0; j<nModules; j++){
            wallChargeTot_mod[j] = 0;
        }


        for(Short_t j=0; j<nWallHitsTot; j++) {
            wallHit = HCategoryManager::getObject(wallHit,wallCat,j);

            wallModuleIndex[j] = wallHit->getCell();
            wallHitTime[j] = wallHit->getTime();
            wallHitCharge[j] = wallHit->getCharge();
            wallHitDistance[j] = wallHit->getDistance();
            isWallHitOk[j] = kFALSE;
            wallHitRing[j] = divider->GetRing(wallModuleIndex[j]);
            ring = wallHitRing[j];
            if (ring==-1) {
                cerr << "Error in short MHWallDivider::GetRing(short i="<<wallModuleIndex[j]<<"): it's returned -1" << endl;
                return 2;
            }
            subevent = divider->GetSubEvent(wallModuleIndex[j]);

            psi = wallHit->getPhi() * D2R;
            wallHitPhi[j] = psi;
            wallHit -> getXYZLab(wallHitX [j], wallHitY [j], wallHitZ [j]);
            thetaFW = wallHit->getTheta();
            wallHitEta[j] = -TMath::Log(TMath::Tan(thetaFW/2.));


            // std::cout << "Ind" << wallModuleIndex[j] << std::endl;
            // std::cout << "Time" << wallHitTime[j] << std::endl;
            // std::cout << "Charge" << wallHitCharge[j] << std::endl;
            // std::cout << "ring" << ring << std::endl;
            // std::cout << "phi" << wallHitPhi[j] << std::endl;
            //cuts by B.Kardan
            hit_beta = wallHitDistance[j]/wallHitTime[j]/299.792458;
            if ( (ring<=4            && wallHitCharge[j]>80 && hit_beta>0.84 && hit_beta<1.) ||
                 ((ring==5||ring==6) && wallHitCharge[j]>85 && hit_beta>0.85 && hit_beta<1.) ||
                 (ring>6             && wallHitCharge[j]>86 && hit_beta>0.80 && hit_beta<1.) ) {

                isWallHitOk[j] = kTRUE;

                wallChargeTot+=wallHitCharge[j];

//                //Charge for rings
//                for (Short_t k = 0; k < nRings; k++){
//                    if (ring == k) wallChargeTot_ring[k]+=wallHitCharge[j];
//                }
//
//                //Charge for modules
//                for (Int_t k = 1; k <= nModules; k++){
//                    if(wallModuleIndex[j] == k) wallChargeTot_mod[k-1]+=wallHitCharge[j];
//                }
            }//B.Kardan cuts for wall hits


            //Charge for rings
            for (Short_t k = 0; k < nRings; k++){
                if (ring == k) wallChargeTot_ring[k]+=wallHitCharge[j];
            }

            //Charge for modules
            for (Int_t k = 1; k <= nModules; k++){
                if(wallModuleIndex[j] == k) wallChargeTot_mod[k-1]+=wallHitCharge[j];
            }

        }//loop over wall hits


        // loop over particle candidates in event
        nTracks = 0;
        nProtons = 0;
        nPiM = 0;
        nPiP = 0;
        if(!candCat) continue;
        Int_t size = candCat->getEntries();
        HParticleCand* cand=0;
        Int_t itr;
        for(Int_t j = 0; j < size; j++) {
            cand = HCategoryManager::getObject(cand,candCat,j);
            if(!cand) continue;
            if(!loop.goodSector(cand->getSector())) { continue;}  // skip inactive sectors
            if(!cand->isFlagBit(kIsUsed)) continue;
            if( cand->getMomentum() == cand->getMomentumOrg() ) continue; //skip tracks with too high pt
            itr = nTracks;
            pid[itr] = cand->getPID();
            p[itr] = cand->getMomentum();
            theta[itr] = cand->getTheta() * D2R;
            pt[itr] = p[itr] * TMath::Sin( theta[itr] );
            metaMass[itr] = cand->getMass();
            Float_t pz = p[itr] * TMath::Cos( theta[itr] );
            Float_t e_lab = TMath::Sqrt( metaMass[itr]*metaMass[itr] + p[itr]*p[itr] );
            rapidity[itr] = 0.5 * TMath::Log( (e_lab + pz)/(e_lab - pz) ) - Y_BEAM;
            eta[itr] = -TMath::Log(TMath::Tan(theta[itr]/2));

            if (pid[itr]>=0) {
                pCorr[itr] = cand->getCorrectedMomentumPID(pid[itr]);        // retrieve corrected mom
                cand->setMomentum(pCorr[itr]);                                   // write it back
                cand->calc4vectorProperties(HPhysicsConstants::mass(pid[itr]));   // sync with lorentz vector
            }

            pCorr[itr] = cand->getMomentum();
            pt_corr[itr] = pCorr[itr] * TMath::Sin( theta[itr] );
            phi[itr] = cand->getPhi()*D2R;
            metaBeta[itr] = cand->getBeta();
            charge[itr] = cand->getCharge();
            mdcNhitsInner[itr] = cand->getNLayer(0);
            mdcNhitsOuter[itr] = cand->getNLayer(1);
            mdcNhits[itr] = cand->getNLayer(2);//inner+outer
            mdcdEdx[itr] = cand->getMdcdEdx();
            tofdEdx[itr] = cand->getTofdEdx();
            DCAxy[itr] = cand->getR();
            DCAz[itr] = cand->getZ();
            chi2all[itr] = cand->getChi2();
            chi2inner[itr] = cand->getInnerSegmentChi2();
            chi2outer[itr] = cand->getOuterSegmentChi2();
            metaQ[itr] = cand->getMetaMatchQuality();
            metaDx[itr] = cand->getRkMetaDx();
            metaDy[itr] = cand->getRkMetaDy();
            metaMatchRadius[itr] = cand->getMetaMatchRadius();
            mdcSecId[itr] = cand->getSector();
            pz = pCorr[itr] * TMath::Cos( theta[itr] );
            Float_t pidMass = (pid[itr]>=0) ? HPhysicsConstants::mass(pid[itr]) : metaMass[itr];
            e_lab = TMath::Sqrt( pidMass*pidMass + pCorr[itr]*pCorr[itr] );
            rapidity_corr[itr] = 0.5 * TMath::Log( (e_lab + pz)/(e_lab - pz) ) - Y_BEAM;

            if( pid[itr]==14 ) nProtons++;
            if( pid[itr]==8 )  nPiP++;
            if( pid[itr]==9 )  nPiM++;

            nTracks++;
            if (nTracks >= maxNTracks) {
                cerr << "ERROR: nTracks > maxNTracks=="<<maxNTracks<<endl;
                return 3;
            }

            //track cuts for MDC Q-vectors
            /*if ( pid[itr]==14 && pt_corr[itr]>250 && pt_corr[itr]<1700 &&
                 TMath::Abs(DCAz[itr]-vZ)<15. && TMath::Abs(DCAxy[itr])<15. &&
                 beta[itr]<1 && metaMass[itr]>600 && metaMass[itr]<1200) */

        } // end cand loop

        tree->Fill();

    } // end eventloop

    tree->Write();
    //out->Save();
    out->Close();

    delete gHades;
    return 0;
}

int main(int argc, char **argv)
{
    TROOT Analysis("Analysis","compiled analysis macro");

    // argc is the number of arguments in char* array argv
    // CAUTION: argv[0] contains the progname
    // argc has to be nargs+1

    TString nevts ;
    TString filenumber;
    switch (argc) {
    case 4:       // just inputfile name + nArgs
        nevts  = argv[3];
        return makeTree(TString(argv[1]),TString(argv[2]), nevts.Atoi());
        break;
    default:
        cerr<<"ERROR: loopDST() : WRONG NUMBER OF ARGUMENTS! TString infile="",TString outfile="", nevents=-1"<<endl;
        return 1; // fail
    }
}
