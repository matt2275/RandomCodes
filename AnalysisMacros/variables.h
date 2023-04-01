#ifndef VARIABLES_H
#define VARIABLES_H

#include <vector>
 using namespace std;
// define your own namespace to hold constants
namespace constants
{
    // constants have internal linkage by default
    double pi { 3.14159 };
    constexpr double avogadro { 6.0221413e23 };
    constexpr double myGravity { 9.2 }; // m/s^2 -- gravity is light on this planet
    vector<int> testvec{1,2,3};
    int testarray[] = {1,2,3};
    map<string,vector<double>> mapvectors={
       {"TEST",{1,2,3}},
       {"TEST2",{1.1,2.2,3.3}},
       {"TestBin",{0,1,2,3,7,10}},
       {"gentrk1_pt",{0,1,2,3,7,10}},
       {"gentrk1_eta",{0,1,2,3,7,10}},
       {"gentrk1_phi",{0,1,2,3,7,10}},
       {"gentrk2_pt",{0,1,2,3,7,10}},
       {"gentrk2_eta",{0,1,2,3,7,10}},
       {"gentrk2_phi",{0,1,2,3,7,10}},
       {"eff_eta",{-2.4,-2.2,-2.0,-1.8,-1.6,-1.2,-.6,0,.6,1.2,1.6,1.8,2.0,2.2,2.4}},
       {"abs_eff_eta",{0,.6,1.2,1.6,1.8,2.0,2.2,2.4}},
       {"eff_pt",{0,.5,1,1.5,2,2.5,3.0,3.5,4,5.0,6,8,12,16,20,25}},
       // {"eff_pt",{0,.2,.4,.6,.8,1.0,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.4,5.0,6,8,12,16,20,25}},
       {"eff_phi",{-3.2,-2.8,-2.4,-2.0,-1.6,-1.2,-.8,-.4,0,.4,.8,1.2,1.6,2.0,2.4,2.8,3.2}},
       
    };
    // mapvectors["test"] = {1.1,4.2,3.2};
    // ... other related constants
}
#endif