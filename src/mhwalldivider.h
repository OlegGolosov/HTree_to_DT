//Alexander Zaytsev
//a.zaytsev@gsi.de
//October 2016

#ifndef MHWallDivider_h
#define MHWallDivider_h

class MHWallDivider {   
  public:
    MHWallDivider() {
    
    }

    short GetRing(short i) {
        if ( i < 0 || i > 301) return -1;
        //beam hole
        if (i==65 || i==66 || i==77 || i==78) return -1;

        //small modules
        if ( i<144 ) {
            if ((i>=52 && i<=55) ||
                (i>=88 && i<=91) ||
                i==64 || i==67 ||
                i==76 || i==79) return 0;
            if ((i>=39 && i<=44) ||
                (i>=99 && i<=104) ||
                (i>=51 && i<=87 && (i-51)%12==0) ||
                (i>=56 && i<=92 && (i-56)%12==0)) return 1;
            if ((i>=26 && i<=33) ||
                (i>=110 && i<=117) ||
                (i>=38 && i<=98 && (i-38)%12==0) ||
                (i>=45 && i<=105 && (i-45)%12==0)) return 2;
            if ((i>=13 && i<=22) ||
                (i>=121 && i<=130) ||
                (i>=25 && i<=109 && (i-25)%12==0) ||
                (i>=34 && i<=118 && (i-34)%12==0)) return 3;
            else return 4;
        }
        //middle modules
        else if (i<208) {
            if ((i<=162 && i>=155) ||
                (i<=196 && i>=189) ||
                (i<=186 && i>=165 && (i-165)%4<2) ) return 5;
            else return 6;
        }
        //large modules
        else if (i<302) {
            if ((i<=238 && i>=232) ||
                (i<=279 && i>=273) ||
                (i<=268 && i>=243 && (i-243)%6<2) ) return 7;
            if ((i<=228 && i>=220) ||
                (i<=291 && i>=283) ||
                (i<=272 && i>=239 && (i-242)%3==0) ||
                 i==231 || i==280 ) return 8;
            else return 9;            
        }
        
        return -1;
    }

    short GetSubEvent(short i) {
        short ring = GetRing(i);
        if (ring<=0) return -1;

        switch (ring) {
            case 6: return 2; //WallRing1_6
            case 5: return 3; //WallRing1_5
            case 4: return 4; //WallRing1_4
            case 3:
            case 2:
            case 1: return 5; //WallRing1_3
            default: return 1; //WallRing7_9
        }

        return -1;
    }

};
#endif
