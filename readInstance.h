/*
 File       : readInstance.h
 Project    : Heuristic for CVPCP
 Developer  : Dagoberto Quevedo
 Date       : Mar 26 2013
 Institution: PISIS - FIME - UANL
 */


#ifndef CVPCP_readInstance_h
#define CVPCP_readInstance_h

#include "global.h"

void    size(bool);
bool    read(const char * path);

//read data instance
bool read(const char * path)
{
    ifstream	data;
    int         v1,i,j,type_read;
    float       v2;
    
    D_MAX = DBL_MIN;
    type_read = 0;
    
    data.open(path);
    
    if (!data.is_open())
        return false;
    
    data >> type;
    data >> id;
    data >> n;
    data >> p;
    data >> bestObj;
    
    if (type == 1)
        type_read = 1;
    else if(type == 2)
        type_read = 2;
    else if(type == 3 or type == 8)
        type_read = 3;
    else if (type >= 4 and type <= 7)
        type_read = 4;
    
    switch (type_read) {
        //Beasley
        case 1:
            data >> v1;
            size(true);
            
            for (i = 0; i < n; ++i)
                s[i] = v1;
            
            for (i = 0; i < n; ++i) {
                data >> v1;
                data >> cxy[i][0];
                data >> cxy[i][1];
                data >> w[i];
            }
            
            for (i = 0; i < n; ++i)
                for (j = i; j < n; ++j) {
                    d[i][j] = i != j ? (int)euclidean(cxy[i][0], cxy[j][0], cxy[i][1], cxy[j][1]) : 0;
                    d[j][i] = d[i][j];
                    if (d[i][j] > D_MAX) D_MAX = d[i][j];
                }
            break;
            
        //GalvaoReVelle
        case 2:
            size(false);
            
            for (i = 0; i < n; ++i){
                data >> v2;
                s[i] = (int)v2;
            }
            
            for (i = 0; i < n; ++i){
                data >> v2;
                w[i]=(int)v2;
            }
            
            for (i = 0; i < n; ++i)
                for (j = 0; j < n; ++j){
                    data >> v2;
                    d[i][j] = (int)v2;
                    if (d[i][j] > D_MAX) D_MAX = d[i][j];
                }
            break;
        
        //Lorena and Reinelt
        case 3:
            size(true);
            
            for (i = 0; i < n; ++i) {
                data >> cxy[i][0];
                data >> cxy[i][1];
                data >> s[i];
                data >> w[i];
            }
            
            for (i = 0; i < n; ++i)
                for (j=i; j < n; ++j) {
                    d[i][j] = i!=j ? (int)euclidean(cxy[i][0], cxy[j][0], cxy[i][1], cxy[j][1]) : 0;
                    d[j][i] = d[i][j];
                    if (d[i][j] > D_MAX) D_MAX = d[i][j];
                }
            
            break;   
        
        //OR-Library A-D
        case 4:
            size(false);
            
            for (i = 0; i < n; ++i)
                for (j = 0; j < n; ++j){
                    data >> d[i][j];
                    if (d[i][j] > D_MAX) D_MAX = d[i][j];
                }
            
            for (i = 0; i < n; ++i)
                data >> w[i];
            
            for (i = 0; i < n; ++i)
                data >> s[i];
    }
    data.close();
    return true;
}

void size(bool coor)
{    
    int     i;
    
    d = new int*[n];
    s = new int [n];
    w = new int [n];
    V = vector<int>(n);
    
    for(i = 0; i < n; ++i){
        V[i] = i;
        d[i] = new int[n];
    }
    
    if(coor){
        cxy = new int*[n];
        for(i = 0; i < n; ++i) 
            cxy[i] = new int[2];
    }
}

#endif
