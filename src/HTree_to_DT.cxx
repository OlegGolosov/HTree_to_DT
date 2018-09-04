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
#include "TLorentzVector.h"

#include <iostream>

#include "DataTreeEvent.h"
//#include "DataTreeTrack.h"
//#include "DataTreeTOFHit.h"
#include "mhwalldivider.h"


using namespace std;

static const int nSimpleTriggers = 15;
const int [nSimpleTriggers] triggerMap =
{
	Particle::kGoodVertexClust, 
	Particle::kGoodVertexCand, 
	Particle::kGoodSTART,
	Particle::kNoPileUpSTART,
	Particle::kNoPileUpMETA,
	Particle::kNoPileUpMDC,
	Particle::kNoFlashMDC,
	Particle::kGoodMDCMult,
	Particle::kGoodMDCMIPSMult,
	Particle::kGoodLepMult,
	Particle::kGoodTRIGGER,
	Particle::kGoodSTART2,
	Particle::kNoVETO,
	Particle::kGoodSTARTVETO,
	Particle::kGoodSTARTMETA
};

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

    //Time
    Int_t time;

		DataTreeEvent DTEvent = new DataTreeEvent ();
    TFile* out = new TFile(outfile.Data(),"RECREATE");
    out->cd();

    TTree* tree = new TTree("DataTree", "HADES Au+Au 1.23 GeV 9gen tree for flow analysis");
    tree->Branch("DTEvent", &DTEvent);

    Int_t entries = loop.getEntries();
    if(nEvents < entries && nEvents >= 0 ) entries = nEvents;
    TString filename;
    Int_t sectors [6];
    MHWallDivider* divider = new MHWallDivider();
		
    for (Int_t i = 0; i < entries; i++) {
				DTEvent -> ClearEvent();
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
        evtInfo = HCategoryManager::getObject (evtInfo, evtInfoCat, 0);

        //get type of trigger
				HEventHeader *header = gHades->getCurrentEvent()->getHeader();
        for (Int_t k = 0; k < nSimpleTriggers; k++) {
						DTEvent -> AddTrigger ();
						DTEvent -> GetTrigger (k) -> SetIsFired (evtInfo -> isGoodEvent (triggerMap [k]));
				}
        for (Int_t k = nSimpleTriggers; k < kNtriggers; k++){
						DTEvent -> AddTrigger ();
						if (header->isTBit(k+11)) DTEvent -> GetTrigger (k) -> SetIsFired(kTRUE);
        }
				
        //get Run number
				DTEvent -> SetRunId ( header->getEventRunNumber() );
				DTEvent -> SetEventId ( header->getEventSeqNumber() );
				DTEvent -> SetEventTimestamp ( header->getTime() );

        //get primary vertex
        HVertex vertexReco = header->getVertexReco();
				DTEvent -> SetVertexPosition(vertexReco.getX(), vertexReco.getY(), vertexReco.getZ(), EnumVertexType::kReconstructedVertex);
				DTEvent -> SetVertexQuality(vertexReco.getChi2(), EnumVertexType::kReconstructedVertex);

       //Track Primary and Hits
//        nRpcClust = evtInfo->getSumRpcMult();
//        nRpcClustCut = evtInfo->getSumRpcMultCut();
//        nRpcHits = evtInfo->getSumRpcMultHit();
//        nTofHits = evtInfo->getSumTofMult();
//        nRpcHitsCut = evtInfo->getSumRpcMultHitCut();
//        nTofHitsCut = evtInfo->getSumTofMultCut();
//        primaryTracks=evtInfo->getSumPrimaryParticleCandMult();
//        selectedTracks=evtInfo->getSumSelectedParticleCandMult();

        // loop over FW hits
        Short_t ring;
        Float_t psi;
        Float_t wallHitBeta, wallHitX, wallHitY, wallHitZ, wallHitPhi, wallHitTheta;
				ushort wallModuleIndex, ring, nWallHitsTot = wallCat->getEntries();
				float wallHitTime, wallHitDistance, wallChargeTot = 0.;
				short wallHitCharge, isWallHitOk;
				HWallHit* wallHit = 0;
				
				wallHit = HCategoryManager::getObject(wallHit,wallCat,0);
				wallHit -> getXYZLab(&wallHitX, &wallHitY, &wallHitZ);
				DTEvent -> SetPsdPosition (0., 0., wallHitZ);
        for(Short_t j=0; j<nWallHitsTot; j++) {
						isWallHitOk = -1;
            wallHit = HCategoryManager::getObject(wallHit,wallCat,j);
            wallModuleIndex = wallHit->getCell();
            wallHit -> getXYZLab(&wallHitX, &wallHitY, &wallHitZ);
            wallHitTime = wallHit->getTime();
            wallHitCharge = wallHit->getCharge();
            wallHitDistance = wallHit->getDistance();
						wallHitPhi = wallHit->getPhi() * D2R;
						wallHitTheta = wallHit->getTheta() * D2R;
            ring = divider -> GetRing(wallModuleIndex);
						wallHitX = wallHitDistance * cos (wallHitTheta) * cos (wallHitPhi);
						wallHitY = wallHitDistance * cos (wallHitTheta) * sin (wallHitPhi);
						wallHitZ = wallHitDistance * cos (wallHitTheta);
            if (ring==-1) {
                cerr << "Error in short MHWallDivider::GetRing(short i="<<wallModuleIndex[j]<<"): it returned -1" << endl;
                return 2;
            }
						//cuts by B.Kardan
            wallHitBeta = wallHitDistance / wallHitTime / 299.792458;
            if ( (ring <= 4            		 && wallHitCharge > 80 && hit_beta > 0.84 && hit_beta < 1.) ||
                 ((ring == 5 || ring == 6) && wallHitCharge > 85 && hit_beta > 0.85 && hit_beta < 1.) ||
                 (ring > 6                 && wallHitCharge > 86 && hit_beta > 0.80 && hit_beta < 1.) ) {

                isWallHitOk = 1;
                wallChargeTot += wallHitCharge;
            }//B.Kardan cuts for wall hits
						
						DTEvent -> AddPSDModule(1);
						DTEvent -> GetPSDModule(j) -> SetId (ring * isWallHitOk);
						DTEvent -> GetPSDModule(j) -> SetPosition(wallHitX, wallHitY, wallHitBeta);
						DTEvent -> GetPSDModule(j) -> SetEnergy(wallHitCharge);
        }//loop over wall hits
				DTEvent -> SetPSDEnergy(wallChargeTot);

        // loop over particle candidates in event
        if(!candCat) continue;
        Int_t size = candCat->getEntries();
        HParticleCand* cand=0;
        Int_t itr, charge;
				DTTrack *track;
				DataTreeTOFHit *hit;
				TLorentzVector trackPar;
				float p, theta, pt, eta, phi, mass;
				
        for(Int_t j = 0; j < size; j++) {
            cand = HCategoryManager::getObject(cand,candCat,j);
            if(!cand) continue;
            if(!loop.goodSector(cand->getSector())) { continue;}  // skip inactive sectors
            if(!cand->isFlagBit(kIsUsed)) continue;
            if( cand->getMomentum() == cand->getMomentumOrg() ) continue; //skip tracks with too high pt ???
            pid = cand -> getPID();
//            chi2inner[itr] = cand->getInnerSegmentChi2();
//            chi2outer[itr] = cand->getOuterSegmentChi2();
//            mdcSecId[itr] = cand->getSector();
						
            if (pid >= 0) {
								mass = HPhysicsConstants::mass(pid);
                p = cand->getCorrectedMomentumPID(pid);        // retrieve corrected mom
                cand->setMomentum(p);                                   // write it back
                cand->calc4vectorProperties(mass);   // sync with lorentz vector
            }
						else {
							mass = cand -> getMass(); // META mass
							p = cand -> getMomentum();
						}
						
            theta = cand -> getTheta() * D2R;
            phi = cand->getPhi()*D2R;
            pt = p * TMath::Sin( theta );
            eta = -TMath::Log(TMath::Tan(theta/2.));
						
						trackPar.SetPtEtaPhiM (pt, eta, phi, mass);
						DTEvent -> AddVertexTrack();
						track = DTEvent -> GetLastVertexTrack ();
						track -> SetMomentum (trackPar);
            track -> SetCharge (cand -> getCharge ());
            track -> SetChi2 (cand -> getChi2 ());
            track -> SetNDF (1.);
            track -> SetNumberOfHits(cand -> getNLayer (0), HADES_constants::kMDC_in);
            track -> SetNumberOfHits(cand -> getNLayer (1), HADES_constants::kMDC_out);
            track -> SetNumberOfHits(cand -> getNLayer (2), HADES_constants::kMDC_all);
            track -> SetdEdx(cand -> getMdcdEdx (), HADES_constants::kMDC_all);
            track -> SetdEdx(cand -> getTofdEdx (), HADES_constants::kMETA);
            track -> SetDCA (cand -> getR (), cand -> getR (), cand -> getZ () - vertexReco.getZ());
						
						DTEvent -> AddTOFHit();
						hit = DTEvent -> GetLastTOFHit();
						hit -> AddRecoTrackId (itr);
						if (cand -> getSystem () == 0) hit -> SetStatus (kRPC);
						else hit -> SetStatus (kTOF);
						hit -> SetTime (cand -> getDistanceToMetaHit () / cand -> getBeta () / 299.792458);
						hit -> SetPathLength(cand -> getDistanceToMetaHit ());
						hit -> SetPosition (cand -> getRkMetaDx (), cand -> getRkMetaDy (), cand -> getMetaMatchRadius ()); // META match qa - NOT POSITION!!!
						hit -> SetCharge(cand -> getCharge ());
						hit -> SetSquaredMass(cand -> getMass2 ());
						hit -> SetSquaredMassError(cand -> getMetaMatchQuality ());
						
            itr++;
        } // end cand loop
				
				cout << "sumTOFRPC = " << evtInfo -> getSumRpcMult() << "/tsumTOFRPCCut = " << evtInfo -> getSumRpcMultCut() << "/tnTOFHits = " << DTEvent -> GetNTOFHits () << endl; 
				cout << "nTracks = " << evtInfo -> getSumPrimaryParticleCandMult() << "/tnTracksSel = " << evtInfo -> getSumSelectedParticleCandMult () << "/tnTOFHits = " << DTEvent -> GetNVertexTracks () << endl;
				
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
