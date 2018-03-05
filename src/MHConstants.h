//Alexander Zaytsev
//a.zaytsev@gsi.de
//October 2016

#ifndef MHConstants_h
#define MHConstants_h

#include <iostream>

using namespace std;

//from treeMaker
const Short_t nWallSubEvents = 6;        //number of FW subevents
const Short_t nWallRandomSubEvents = 10; //number of random FW subevents
const Short_t nMdcSubEvents = 3;         //number of MDC subevents 
const Short_t nHarmonics = 1;            //number of harmonics for q-vectors

const Short_t maxNWallHits = 200;
const Short_t maxNTracks = 200;

/*const Int_t N_MULT_TYPES = 3; //nTracks;   nProtons+nPiP+nPiM;   wallQw[0]
const Double_t multBinning[N_MULT_TYPES][3] = {40, 0, 120, //nTracks
                                               45, 0, 90,  //nProtons+nPiP+nPiM
                                               45, 0, 90}; //wallQw[0]
const int      nBins[N_MULT_TYPES] = {40, 45, 45}; //CAUTION: it must be equal to multBinning[i][0]
const int      nBinsMaximal = 45; //CAUTION: must be equal to maximal number in nBins[N_MULT_TYPES]
const int      multTypeForFlow = 0;*/

//temporary solution
const Int_t    N_MULT_TYPES = 1; //nTracks;
const Double_t multBinning[N_MULT_TYPES][3] = {40, 0, 120};
const int      nBins[N_MULT_TYPES] = {40}; //CAUTION: it must be equal to multBinning[i][0]
const int      nBinsMaximal = 40; //CAUTION: must be equal to maximal number in nBins[N_MULT_TYPES]
const int      multTypeForFlow = 0;


const Double_t ptBinning[3] = {20, 200, 1700}; //bin width: 75 MeV/c
const Double_t yBinning[3] = {20, -0.7, 0.8};


enum QSubEvent {
    eWall,
    eWallRing7_9,
    eWallRing1_6,
    eWallRing1_5,
    eWallRing1_4,
    eWallRing1_3,
    eWallRandomA,
    eWallRandomB,
    eWallSadovskyA,
    eWallSadovskyB,
    eWallChequerwiseA,
    eWallChequerwiseB,
    eWallRing7_9RandomA,
    eWallRing7_9RandomB,
    eWallRing1_6RandomA,
    eWallRing1_6RandomB,
    eMdcNegRapidity,
    eMdcHighNegRapidity,
    eMdcPosRapidity,
    N_Q_SUBEVENTS
};
//ATTENTION: change it if you're going to change enum QSubEvent
const QSubEvent firstMdcSub = eMdcNegRapidity;
const QSubEvent lastMdcSub = eMdcPosRapidity;

enum UType {
    eProtons,
    ePiP,
    ePiM,
    N_U_TYPES
};


void GetSubEventName(int isub, char* subName) {    
    switch (isub) {
    case eWall:
        sprintf(subName, "Wall");
        return;
    case eWallRing7_9:
        sprintf(subName, "WallRing7_9");
        return;
    case eWallRing1_6:
        sprintf(subName, "WallRing1_6");
        return;    
    case eWallRing1_5:
        sprintf(subName, "WallRing1_5");
        return;
    case eWallRing1_4:
        sprintf(subName, "WallRing1_4");
        return;
    case eWallRing1_3:
        sprintf(subName, "WallRing1_3");
        return;
    case eWallRandomA:
        sprintf(subName, "WallRandomA");
        return;
    case eWallRandomB:
        sprintf(subName, "WallRandomB");
        return;
    case eWallSadovskyA:
        sprintf(subName, "WallSadovskyA");
        return;
    case eWallSadovskyB:
        sprintf(subName, "WallSadovskyB");
        return;   
    case eWallChequerwiseA:
        sprintf(subName, "WallChequerwiseA");
        return;
    case eWallChequerwiseB:
        sprintf(subName, "WallChequerwiseB");
        return;
    case eWallRing7_9RandomA:
        sprintf(subName, "WallRing7_9RandomA");
        return;
    case eWallRing7_9RandomB:
        sprintf(subName, "WallRing7_9RandomB");
        return;
    case eWallRing1_6RandomA:
        sprintf(subName, "WallRing1_6RandomA");
        return;
    case eWallRing1_6RandomB:
        sprintf(subName, "WallRing1_6RandomB");
        return;
    case eMdcNegRapidity:
        sprintf(subName, "MdcNegRapidity");
        return;
    case eMdcHighNegRapidity:
        sprintf(subName, "MdcHighNegRapidity");
        return;
    case eMdcPosRapidity:
        sprintf(subName, "MdcPosRapidity");
        return;
    default:
        cerr << "ERROR: QSubEvent=="<<isub<<" >= N_Q_SUBEVENTS in GetSubEventName()!" << endl;
    }
    sprintf(subName, "ERROR");
    return;
}

void GetSubEventTag(int isub, char* subName) {    
    switch (isub) {
    case eWall:
        sprintf(subName, "FW");
        return;
    case eWallRing7_9:
        sprintf(subName, "Ring7-9");
        return;
    case eWallRing1_6:
        sprintf(subName, "Ring1-6");
        return;    
    case eWallRing1_5:
        sprintf(subName, "Ring1-5");
        return;
    case eWallRing1_4:
        sprintf(subName, "Ring1-4");
        return;
    case eWallRing1_3:
        sprintf(subName, "Ring1-3");
        return;
    case eWallRandomA:
        sprintf(subName, "Random");
        return;
    case eWallRandomB:
        sprintf(subName, "Random");
        return;
    case eWallSadovskyA:
        sprintf(subName, "Sadovsky");
        return;
    case eWallSadovskyB:
        sprintf(subName, "Sadovsky");
        return;   
    case eWallChequerwiseA:
        sprintf(subName, "Chequerwise");
        return;
    case eWallChequerwiseB:
        sprintf(subName, "Chequerwise");
        return;
    case eWallRing7_9RandomA:
        sprintf(subName, "Ring7_9Random");
        return;
    case eWallRing7_9RandomB:
        sprintf(subName, "Ring7_9Random");
        return;
    case eWallRing1_6RandomA:
        sprintf(subName, "Ring1_6Random");
        return;
    case eWallRing1_6RandomB:
        sprintf(subName, "Ring1_6Random");
        return;
    case eMdcNegRapidity:
        sprintf(subName, "Mdc-");
        return;
    case eMdcHighNegRapidity:
        sprintf(subName, "Mdc--");
        return;
    case eMdcPosRapidity:
        sprintf(subName, "Mdc+");
        return;
    default:
        cerr << "ERROR: QSubEvent=="<<isub<<" >= N_Q_SUBEVENTS in GetSubEventName()!" << endl;
    }
    sprintf(subName, "ERROR");
    return;
}


void GetUTypeName(int utype, char* typeName) {
    switch (utype) {
    case eProtons:
        sprintf(typeName, "uProt");
        return;
    case ePiP:
        sprintf(typeName, "uPiP");
        return;
    case ePiM:
        sprintf(typeName, "uPiM");
        return;
    default:
        cerr << "ERROR:  UType >= N_U_TYPES in GetUTypeName()!" << endl;
    }
    sprintf(typeName, "ERROR");
    return;
}


void GetMultTypeName(int imt, char* typeName) {
    switch (imt) {
    case 0:
        sprintf(typeName, "NTracks");
        return;
    case 1:
        sprintf(typeName, "NPiProtons");
        return;
    case 2:
        sprintf(typeName, "NWallHits");
        return;
    default:
        cerr << "ERROR:  imt is too big in GetMultTypeName()!" << endl;
    }
    sprintf(typeName, "ERROR");
    return;
}


bool IsUQCorrelated(int utype, int isub) {
    if (isub >= firstMdcSub) return true; //if you decided to change it, then you need to change "firstMdcSub-behavior" in for-loops (e.g. in WriteHistos)
    return false;
}

void Swap(int* a, int* b) {
    *a = *a + *b;
    *b = *a - *b;
    *a = *a - *b;
}

bool IsQQCorrelated(int isub, int isub2) {
    if (isub==isub2) return true;
    
    if (isub>isub2) Swap(&isub,&isub2);
    
    //isub2>isub now
    switch (isub2) {
    case eWall:
        return true;
    case eWallRing7_9:
        return true;
    case eWallRing1_6:
    case eWallRing1_5:
    case eWallRing1_4:
    case eWallRing1_3:
        if (isub==eWallRing7_9) return false;
        else return true;
    case eWallRandomA:
        return true;
    case eWallRandomB:
        if (isub==eWallRandomA) return false;
        else return true;
    case eWallSadovskyA:
        return true;
    case eWallSadovskyB:
        if (isub==eWallSadovskyA) return false;
        else return true;
    case eWallChequerwiseA:
        return true;
    case eWallChequerwiseB:
        if (isub==eWallChequerwiseA) return false;
        else return true;
    case eWallRing7_9RandomA:
        return true;
    case eWallRing7_9RandomB:
        if (isub==eWallRing7_9RandomA) return false;
        else return true;
    case eWallRing1_6RandomA:
        return true;
    case eWallRing1_6RandomB:
        if (isub==eWallRing1_6RandomA) return false;
        else return true;
    case eMdcNegRapidity:
    case eMdcHighNegRapidity:
        if (isub<eWallRandomA) return false;
        else return true;   
    case eMdcPosRapidity:
        if (isub<eWallRandomA || isub>=firstMdcSub) return false;
        else return true;
    default:
        cerr << "ERROR: isub2=="<<isub2<<" >= N_Q_SUBEVENTS in IsQQCorrelated()!" << endl;
    }    
    return true;
}

void ArrangeAscending(int* a, int* b, int* c) {
    if (*a > *b) Swap(a,b);
    if (*b > *c) Swap(b,c);
    if (*a > *b) Swap(a,b); 
}

#endif