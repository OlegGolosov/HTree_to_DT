//Alexander Zaytsev
//a.zaytsev@gsi.de
//June 2016

#ifndef MHQVector_h
#define MHQVector_h

#include "TROOT.h"

class MHQVector : public TObject {  
  public:
    Float_t x;
    Float_t y;
    Float_t weight;


    //__________________________________________________________________
    //methods

    MHQVector() {
        Init();
    }

    
    ~MHQVector() {
    }
    
    
    MHQVector(Float_t _x, Float_t _y, Float_t _w) {
        x = _x;
        y = _y;
        weight = _w;
    }
    

    void Init() {
        x = 0;
        y = 0;
        weight = 0;
    }
    
    
    ClassDef(MHQVector, 1);
    
};
#endif
