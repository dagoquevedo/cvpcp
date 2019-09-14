/*
 File       : main.cpp
 Project    : Heuristic for CVPCP
 Developer  : Dagoberto R. Quevedo-Orozco
 Date       : Nov 13 2013
 Institution: PISIS - FIME - UANL
 */


#include "global.h"
#include "readInstance.h"
#include "IGVND.h"

using namespace std;

int main (int argc, const char * argv[])
{
    if(argc < 5) {
        cout<<"Missing parameters.\n";
        return -1;
    }
    
    if(!read(argv[1])) {
        cout<<"File not found.\n";
        return -1;
    }
    
    r_max = atoi(argv[2]);
    alpha = atof(argv[3]);
    e     = atoi(argv[4]);
    
    if (r_max == 0 || alpha == 0 || e == 0) {
        cout<<"Invalid parameters.\n"<<endl;
        return -1;
    }
    
    time_start = getCPUTime();
    
    seed = get_seed();
    srand(seed);
    
    result = IGVND(r_max, alpha, e, 1, 1, 1);
    
    time_stop = getCPUTime();
    
    print(result);
    
    if(argc > 5)
        save(argv[5], result);
    
    finalize();
    
    return 0;
}
