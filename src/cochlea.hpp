#pragma once

#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

typedef double flotante;
typedef vector<flotante> flotante_vector;
typedef flotante* flotante_puntero;

class Cochlea_t {

    public:
        Cochlea_t(){}
        
        void gaussian_elimination_init();
        void set_arrays();
        
        void setup();
        
        void operator()( flotante_vector &X , flotante_vector &dXdt , flotante t );
        
        flotante_vector g;
        flotante_vector a;
        flotante_vector b;
        flotante_vector c;
        flotante_vector p;
        flotante_vector k;

        flotante_vector X;
        flotante_puntero X_t;
        
        flotante_puntero stimulus, weight;
        flotante_puntero ddivm, sdivm; 
        
        flotante dt;
        
        int N, n_t, n_osc;
        
        flotante alpha, beta, asb, gamma, delta, eps,zeta, N_delta, eta, fluid, nu, m0;
    
        flotante F0;
        flotante dv;
            
};

struct Observer_t {


    int c,k;
    flotante_puntero &m_X_t;
    flotante_puntero &m_times;
    int m_n_eq;
    int m_dec;

    Observer_t( flotante_puntero &X_t , flotante_puntero &times, int n_eq, int decimate=1 )
              : m_X_t( X_t ) , m_times( times ) , m_dec(decimate), m_n_eq(n_eq) {c = 0;k=0; }

    void operator()( const flotante_vector &X , flotante t )
    {
        if (k%m_dec==0)
        {
            // cout << k << " " << c << endl;
            m_times[c] = t;
            for(int i = 0; i<m_n_eq; i++)
            {
                m_X_t[c*m_n_eq+i] = X[i];
            }
            c++;
        }
        k++;
        
//            cout << c << endl;
    }

};

struct Observer_variable_t {

    vector<flotante_vector> &m_X_t;
    vector<flotante> &m_times;
    
    Observer_variable_t( vector<flotante_vector> &X_t , vector<flotante> &times ): m_X_t( X_t ) , m_times( times ){ }

    void operator()( const flotante_vector &X , flotante t )
    {
        m_times.push_back(t);
        m_X_t.push_back(X);
    }

};

inline flotante dampingNL(flotante x,flotante gam)
{

    flotante aux = sqrt(0.5/gam);
    flotante aux2 = fabs(x);

    if (aux2>2*aux)
        return 1;
    if (aux2<aux)
        return gam*x*x;
    else
        return 1-gam*(aux2-2*aux)*(aux2-2*aux);

}

inline flotante dampingNL2(flotante x,flotante gam, flotante delt)
{

    flotante xp = sqrtf( 1.0/ ( -(gam+delt) + (gam+delt)*(gam+delt)/delt ) );
    flotante x0 = (gam/delt+1)*xp;

    flotante aux2 = fabs(x);

    if (aux2>xp)
        return 1;
    if (aux2<xp)
        return gam*x*x;
    else
        return 1-delt*(aux2-x0)*(aux2-x0);

}

inline flotante clip(flotante n, flotante lower, flotante upper) {
  return max(lower, min(n, upper));
}

inline flotante lin_interp(flotante_puntero S, flotante k)
{
    size_t k_1  = floor(k);
    size_t k_2  = k_1 + 1;
    flotante ee = k - k_1;
    return S[k_1]*(1-ee) + S[k_2]*ee;
}


    
//template <typename T>
//T clip(const T& n, const T& lower, const T& upper) {
//  return std::max(lower, std::min(n, upper));
//}
    
