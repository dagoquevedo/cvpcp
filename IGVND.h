/*
 File       : IGVND.h
 Project    : Heuristic for CVPCP
 Developer  : Dagoberto Quevedo
 Date       : Apr 13 2013
 Institution: PISIS - FIME - UANL
 */

#ifndef CVPCP_IGVND_h
#define CVPCP_IGVND_h

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

#include "global.h"

//Functions
double   gamma(int, int);
bool     isBetter(solution, solution);
double   phi(int, int, solution);
double   r(int, solution);

double          a(vector<int>);
int             c(solution, int);
double          f(solution, int);
bool            feasible(solution);
int             l(int, solution);
int             omega(short, int, solution, int);
vector<int>     B(solution);
vector<int>     J(solution);
bool            duplicate(solution S);

//Methods
solution        Construction(int, vector<int>, bool);
solution        IGLS(float, solution);
solution        VND(solution S);
solution        N_1(solution S);
solution        N_2(solution S);
solution        Shake(solution, int);
void            updateYield(solution &);
void            centerAllSubset(solution &);
solution        IGVND(int, float, int, bool, bool, bool);

#pragma mark Funtions
    // gamma(i,j): Greedy function
    double gamma(int i, int j)
    {
        return s[j] * d[i][j];
    }

    // r(k) residual capacity of k partition
    double r(int k, solution S)
    {
        return s[S.P[k]] - S.A[k];
    }

    double phi(int j, int k, solution S)
    {
        return max(d[j][S.P[k]]/float(D_MAX + 1),-1 * r(k, S) + w[j]);
    }

    /* Omega function: k more near to j
        1: \arg\min_{k\in K} \phi(j,k)
        2: \arg\min_{k\in K} \d_{jc(k))
    */
    int omega(short type, int j, solution S, int index = 1)
    {
        map<int,double> Q;
        int k, n, q;
        double D, v_min = DBL_MAX;
        
        for (k = 0; k < S.P.size(); ++k)
        {
            switch (type) {
                case 1:
                    // phi(j,k)
                    D = max(d[j][S.P[k]]/float(D_MAX + 1),-1 * (s[S.P[k]] - S.A[k]) + w[j]);
                    break;
                case 2:
                    D = d[j][S.P[k]];
                    break;
                default:
                    break;
            }
            
            if (index == 1) {
                if (D < v_min) {
                    v_min = D;
                    q = k;
                }
            }
            else
                Q[k] = D;
        }
        
        if (index == 1)
            return q;
        else
        {
            vector<pair<int, double> > B = sortMapDouble(Q);
        
            for (n = 0; n < (int)B.size(); ++n)
                if (n == index - 1)
                    return B[n].first;
        
            return B[n].first;
        }
    }

    // c(k): get center node in k-partition
    int c(solution S, int k)
    {
        int     d_min = INT_MAX, d_max = 0, v = -1;
        double  A = S.A[k];
        
        for (vector<int>::iterator j = S.X[k].begin(); j != S.X[k].end(); ++j)
        {
            d_max = 0;
            if (s[*j] - A >= 0)
                for (vector<int>::iterator i = S.X[k].begin(); i != S.X[k].end(); ++i)
                    if (d[*j][*i] > d_max)
                        d_max = d[*j][*i];
            
            if (d_max != 0)
                if (d_max < d_min) {
                    d_min = d_max;
                    v = *j;
                }
        }
        
        if (v != -1)
            return v;
        else
            return S.P[k];
    }

    // f(S) function: evaluation solution S
    double f(solution S, int k = -1)
    {
        double f_max = 0;
        double D;
        vector<int>::iterator j;
        if (k == -1)
            for (k = 0; k < p; ++k)
            {
                for (j = S.X[k].begin(); j != S.X[k].end(); ++j)
                {
                    D = d[*j][S.P[k]];
                    if(D > f_max)
                        f_max = D;
                }
            }
        else
            for (j = S.X[k].begin(); j != S.X[k].end(); ++j)
            {
                D = d[*j][S.P[k]];
                if(D > f_max)
                    f_max = D;
            }
        
        return f_max;
    }

    // B(S): bottleneck subsets in S,
    vector<int> B(solution S)
    {
        vector<int> K;
        for (int k = 0; k < p; ++k)
            if (f(S, k) == S.f)
                K.push_back(k);
        
        return K;
    }

    // J(S): the demand nodes with maximum distance from the active center node in S
    vector<int> J(solution S)
    {
        vector<int> U;
        for (int k = 0; k < p; ++k)
            for (vector<int>::iterator j = S.X[k].begin(); j != S.X[k].end(); ++j)
                if (d[*j][S.P[k]] == S.f)
                    U.push_back(*j);

        return U;        
    }

    // l(j): k-partition of j node
    int l(int j, solution S)
    {
        for (int k = 0; k < p; ++k)
            if(indexOf(S.X[k], j) != -1)
                return k;
        
        return -1;
    }

    // feasible(S): check the feasibility solution
    bool feasible(solution S)
    {
        for (int k = 0; k < p; ++k)
            if (r(k, S) < 0)
                return false;
    
        return true;
    }

    // Update status of a partition S
    void updateYield(solution &S)
    {
        S.f = f(S);
        S.F = feasible(S);
        S.B = B(S);
        S.J = J(S);
    }

#pragma mark Methods

    // Construction stages
    solution Construction(int p, vector<int> V, bool partial=false)
    {
        int         i, k;
        float       sum_prob = 0.0;
        float       sum_mont = 0.0, pos;
        vector<int> P,U,W;
        solution    S;
        
        map<int,double> G;
        vector<int>::iterator j;
        
        S.A = vector<int>(p);
        
        for (k = 0; k < p; ++k)
            S.A[k] = 0;
        
        
        //Stage A
        copy(V.begin(), V.end(), back_inserter(U));
        i=random_choice(U);
        P.push_back(i);
        U.erase(find(U.begin(), U.end(), i));
        
        for (j = U.begin(); j != U.end(); ++j) {
            G[*j] = gamma(i,*j);
            sum_prob += G[*j];
        }

        while (P.size() < p) {
            sum_mont = 0.0;
            pos = random_float(0.0, 1.0);
            
            for (j=U.begin(); j != U.end(); ++j){
                sum_mont+= G[*j]/sum_prob;
                if(sum_mont >= pos || (pos-sum_mont)/pos < exp(-6)){
                    i = *j;
                    break;
                }
            }

            P.push_back(i);
            U.erase(find(U.begin(), U.end(), i));
            G.erase(i);
            
            sum_prob = 0.0;
            for (j=U.begin(); j != U.end(); ++j) {
                G[*j] = min(G[*j], gamma(i, *j));
                sum_prob += G[*j];
            }
        }
        
        U.clear();
        G.clear();
        
        
        //Stage B
        S.X=vector<vector<int> >(p);
        S.P = P;

        for (k = 0; k < p; ++k)
        {
            S.X[k].push_back(S.P[k]);
            S.A[k] += w[S.P[k]];
        }
        
        W = setDifference(V,P);

        for (j = W.begin(); j != W.end(); ++j) {
            k = omega(1, *j, S);
            S.X[k].push_back(*j);
            S.P[k] = c(S,k);
            S.A[k] += w[*j];
        }

        
        W.clear();
        
        if (!partial)
            updateYield(S);
    
        return S;
    }

    // IGLS
    solution IGLS(float alpha, solution S)
    {
        int                     i, k ,m , N;
        vector<int>             W;
        map<int,double>         U;
        map<int,double>         T;
        vector<int>::iterator   j;
        float                   sum_t,sum_mont,pos;
        
        for (k=0; k<p; ++k) {
            copy(S.X[k].begin(), S.X[k].end(), back_inserter(W));
            
            W.erase(find(W.begin(), W.end(), S.P[k]));
            
            m = 0;
            N = (int)W.size();
            
            if (N > 0) {
                sum_t = 0.0;
                for (j = W.begin(); j != W.end(); ++j) {
                    T[*j]  = d[*j][S.P[k]] > 0 ? d[*j][S.P[k]] : 1;
                    sum_t += T[*j];
                }
                
                while (m/float(N) < alpha && W.size() > 0) {
                    i = -1;
                    pos = random_float(0.0, 1.0);
                    sum_mont = 0.0;
                    
                    for (j = W.begin(); j != W.end(); ++j) {
                        sum_mont += T[*j]/sum_t;
                        if(sum_mont >= pos || (pos - sum_mont)/pos < exp(-6)){
                            i=*j;
                            break;
                        }
                    }
                    sum_t -= T[i];
                    
                    U[i] = d[i][S.P[k]];
                    S.X[k].erase(find(S.X[k].begin(), S.X[k].end(), i));
                    W.erase(find(W.begin(), W.end(), i));
                    S.A[k] -= w[i];
                    
                    m++;
                };
                T.clear();
            }
            W.clear();
        }
    
        vector<pair<int, double> > Q=sortMapDouble(U);
        
        for (m = (int)Q.size() - 1; m >= 0 ; --m) {
            i = Q[m].first;
            k = omega(1, i, S);
            S.X[k].push_back(i);
            S.A[k] += w[i];
        }
        
        U.clear();
        Q.clear();
        
        centerAllSubset(S);
        updateYield(S);
        return S;
    }

    // VND
    solution VND(solution S)
    {
        int k = 1;
        solution So;
        
        while (k <= 2) {
            switch (k) {
                case 1:
                    So = N_1(S);
                    break;
                case 2:
                    So = N_2(S);
                    break;
            }
            if (isBetter(So, S)) {
                S = So;
                k = 1;
            }
            else
                k ++;
        }
        return S;
    }

    // N_1: Reinsertion
    solution N_1(solution S)
    {
        vector<int> Y;
        double      psi, psi_max;
        int         q_,k,q;
        
        Y = S.J;
        
        //Neighborhood by bottleneck
        for (vector<int>::iterator i = Y.begin(); i != Y.end(); ++i) {
            psi_max = DBL_MIN;
            
            q_ = -1;
            k  = l(*i, S);
            
            for (q = 0; q < p; ++q) {
                if (k != q && d[*i][S.P[q]] < S.f)
                {
                    if (r(q,S) - w[*i] >= 0) {
                        psi = d[*i][S.P[k]] - d[*i][S.P[q]];
                        if(psi > psi_max) {
                            psi_max = psi;
                            q_ = q;
                        }
                    }
                }
            }
            
            //Move
            if (q_!=-1) {
                
                S.X[k].erase(find(S.X[k].begin(), S.X[k].end(), *i));
                S.X[q_].push_back(*i);
                S.A[k]  -= w[*i];
                S.A[q_] += w[*i];
            }
        }
        
        Y.clear();
        
        centerAllSubset(S);
        updateYield(S);
        return S;
    }

    // Exchange
    solution N_2(solution S)
    {
        vector<int> W,Y;
        double      psi, psi_max;
        int         j_,q_,k,q;

        Y = S.J;
        
        //Neighborhood by bottleneck
        for (vector<int>::iterator i = Y.begin(); i != Y.end(); ++i) {
            psi_max = DBL_MIN;
            j_ = -1;
            q_ = -1;
            
            k=l(*i, S);
            for (q = 0; q < p; ++q) {
                if(k != q && d[*i][S.P[q]] < S.f) {
                    //W=setDifference(S.X[q], Y);
                    copy(S.X[q].begin(), S.X[q].end(), back_inserter(W));
                    W.erase(find(W.begin(), W.end(), S.P[q]));
                    
                    for (vector<int>::iterator j = W.begin(); j!= W.end(); ++j) {
                        if (d[*j][S.P[k]] < S.f && (r(k,S) - w[*i] + w[*j]) >= 0
                                                && (r(q,S) - w[*j] + w[*i]) >= 0) {
                            
                            psi = (d[*i][S.P[k]] - d[*i][S.P[q]])+
                                  (d[*j][S.P[q]] - d[*j][S.P[k]]);
                            
                            if(psi > psi_max) {
                                psi_max = psi;
                                q_ = q;
                                j_ = *j;
                            }
                        }
                    }
                    W.clear();
                }
            }
            
            //Move
            if (q_ != -1) {
                S.X[k].erase(find(S.X[k].begin(),S.X[k].end(),*i));
                S.X[q_].erase(find(S.X[q_].begin(),S.X[q_].end(),j_));
                S.X[k].push_back(j_);
                S.X[q_].push_back(*i);
                
                S.A[k]  -= w[*i];
                S.A[k]  += w[j_];
                S.A[q_] -= w[j_];
                S.A[q_] += w[*i];
            }
        }
        
        Y.clear();
        centerAllSubset(S);
        updateYield(S);
        return S;
    }

    // Shake
    solution Shake(solution S, int e)
    {
        solution    Sp;
        set<int>    K;
        vector<int> W;
        vector<int> Y = S.J;
        int         c, q;
        
        vector<int>::iterator       i;
        set<int>::reverse_iterator  k;
        
        for (i = Y.begin(); i != Y.end(); ++i) {
            c = omega(2, *i, S, 1);
            if (c == l(*i, S)) {
                K.insert(c);
                for (q = 2; q <= e; ++q) {
                    c = omega(2, *i, S, q);
                    K.insert(c);
                }

                for (k = K.rbegin(); k != K.rend(); ++k) {
                    W=setUnion(W, S.X[*k]);
                    S.X.erase(S.X.begin() + *k);
                    S.P.erase(S.P.begin() + *k);
                    S.A.erase(S.A.begin() + *k);
                }
                    
                //Solver a e-center problem
                Sp=Construction((int)K.size(), W, true);
                for (q = 0; q < Sp.X.size(); ++q) {
                    S.X.push_back(Sp.X[q]);
                    S.P.push_back(Sp.P[q]);
                    S.A.push_back(Sp.A[q]);
                }
                    
                W.clear();
                K.clear();
            }
        }
        
        Y.clear();
        updateYield(S);
        return S;
    }

    void centerAllSubset(solution &S)
    {
        for (int k = 0; k < p; ++k)
            S.P[k] = c(S, k);
    }

#pragma mark Improvement criteria

    // Evaluation function
    bool isBetter(solution So, solution S)
    {
        return (So.f < S.f || (So.f == S.f  && So.B.size() <= S.B.size()
                                            && So.J.size() <  S.J.size()))
                                            && So.F;
    }

    bool duplicate(solution S)
    {
        for (int k = 0; k < p; ++k) {
            for (int l = 0; l < p; ++l) {
                if (k != l) {
                    if (setIntersection(S.X[k], S.X[l]).size()>0) {
                        return true;
                    }
                }
            }            
        }
        return false;
    }

#pragma mark Main IGVND

    // Main IGVND function
    solution IGVND(int r_max, float alpha, int e, bool a_IGLS = true, bool a_VND = true, bool a_Shake = true)
    {
        solution S, S_best;
        
        S = Construction(p, V);

        if (a_VND)
            S = VND(S);
        
        int  i = 0;
        i_best = i;
        S_best = S;
        
        if(__DEBUG__)
            print(S);
        
        while (i++ < r_max) {
            if (a_IGLS)
                S = IGLS(alpha, S);
            if (a_VND)
                S = VND(S);
            
            if (isBetter(S, S_best)) {
                S_best = S;
                i_best = i;
                if(__DEBUG__)
                    print(S);
            }
            else {
                if (a_Shake)
                    S = Shake(S, e);
            }
        }
        return S_best;
    }

#endif
