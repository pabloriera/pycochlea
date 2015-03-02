
#include "cochlea.hpp"

void Cochlea_t::setup()
{
    N = n_osc;
    set_arrays();
    gaussian_elimination_init();
    N_delta = N*delta;

}

void Cochlea_t::gaussian_elimination_init()
{
    b[0] = -alpha-1;
    c[0] = alpha;


    //alpha = mass/2/dx/density
    //beta = dx/ height

    //asb = mass height
    //      -----------
    //      dx^2 demsity

    asb =  alpha/beta;

    for (int i = 1; i < N; i++)
    {
        a[i] = asb;
        b[i] = -2*asb-1;
        c[i] = asb;
    }
    
    b[N-1] = 1.0;
    a[N-1] = 0.0;


    c[0] = c[0] / b[0];
 
    /* loop from 1 to N - 1 inclusive */
    for (int i = 1; i < N; i++){
        k[i] = 1.0f / (b[i] - a[i] * c[i - 1]);
        c[i] = c[i] * k[i];
    }
 
}

void Cochlea_t::operator() (  flotante_vector &X , flotante_vector &dXdt , flotante t )
{
    /*
     * X  = vector with current state for in out
     * dX/dT = derivative vector
     * t    = time
     
     * Systems equations
     * X[even] = Y; X[odd] = V
     * dY/dt = V;
     * dV/dt = q - g;
     */

    flotante _X,X_;


    //linear interpolation on FORCE or STIMULUS
    F0 = lin_interp( stimulus, t / dt );
    
    //COMPUTE g = d*V + s*Y    

    for(int i=0;i< N;i++)
    {
//        g[i] = F0*weight[i] + (dv + dvy2*dampingNL2( X[i], gamma, delta ) ) * X[i+N] * ddivm[i] +  X[i] * sdivm[i] ;
        
//        _X =  -lin_interp( X.data(), clip(i - N*delta,0,N-1) )*zet;
//        X_ =  lin_interp( X.data(), clip(i + N*delta,0,N-1) )*eps;

        if(i - N_delta >=0 )
            _X =  lin_interp( X.data(), i - N_delta )*zeta;
        else
            _X = 0;
        
        if( i + N_delta < N-1 )
            X_ = lin_interp( X.data(), i + N_delta )*eps;
        else
        {
            X_ = 0;
            // cout << i << " " << N_delta <<" "<< endl;
        }
            
            // X_ = X[i]*eps*(N-i)/(float)N ;
            // X_ = lin_interp( X.data(),  N-1 )*eps*(N-i)/(float)N ;
        

        if (nu!=0)
            g[i] = X[i+N] * ddivm[i] +  X[i] * sdivm[i] + nu*( tanh(_X/nu) + tanh( X_/nu ) )*sdivm[i] ;
        else if (eta==0)
            g[i] = (  1 + gamma*X[i]*X[i] ) * X[i+N]* ddivm[i]  +  ( X[i] + _X + X_ ) * sdivm[i];
        else
            g[i] = (  1+ gamma*dampingNL( X[i], eta )  ) * X[i+N]* ddivm[i]  +  ( X[i] + _X + X_ ) * sdivm[i];
        

        // g[i] = (  1+ gamma*dampingNL( X[i], eta )  ) * X[i+N]* ddivm[i]  +  ( X[i] + _X + X_ ) * sdivm[i];
        // g[i] = (  1+ gamma*dampingNL( X[i], eta )  ) * X[i+N]* ddivm[i] +  X[i] * sdivm[i] + fluid*( tanh(_X/fluid) + tanh( X_/fluid ) )*sdivm[i] ;

        // g[i] = X[i+N]* ddivm[i]  +  ( X[i] + (_X + X_)*( 1 - gamma*dampingNL( X[i], eta ) ) ) * sdivm[i];

        // g[i] = X[i+N] * ddivm[i] +  X[i] * sdivm[i] + gamma*( tanh(_X/gamma) + tanh( X_/gamma ) )*sdivm[i] ;




        // g[i] = F0*weight[i] + X[i+N] * ddivm[i] +  X[i] * sdivm[i];
        
        // g[i] = X[i+N] * ddivm[i] +  X[i] * sdivm[i] + gamma*( tanh(_X*sdivm[i]/gamma) + tanh( X_*sdivm[i]/gamma ) ) ;

        

        // g[i] = F0*weight[i] + (1 + dv*dampingNL( X[i], gamma ) - dv*dampingNL( _X+X_, gamma )  ) * X[i+N] * ddivm[i] + ( X[i] ) * sdivm[i];
        // g[i] = F0*weight[i] + ( dv + dvy2*X[i]*X[i] ) * X[i+N] * ddivm[i] +  X[i] * sdivm[i] ;
    }
    
    if(fluid)
    {
        //COPMUTE p solving Ap = (q-g);
        p[0] = (F0*weight[0]-g[0]) / b[0];
    
        for (int i = 1; i < N; i++)
            p[i] = (F0*weight[i] - g[i] - a[i] * p[i - 1]) * k[i];
    
//        p[N-1] = - asb * p[N - 2] * k[N-1];
    
        ///* loop from N - 2 to 0 inclusive */
        for (int i = N - 2; i >=0; i--)
            p[i] = p[i] - c[i] * p[i + 1];
    }
    //COMPUTE derivates of vector field
    
    dXdt[0] = X[0+N];
    dXdt[0+N]= (p[0] - g[0])*m0;

    for(int i=1; i<N; i++)
    {
        dXdt[i] = X[i+N];
        dXdt[i+N]= p[i] - g[i];
    }

}

void Cochlea_t::set_arrays()
{
    g.resize(N,0);
    a.resize(N,0);
    b.resize(N,0);
    c.resize(N,0);
    p.resize(N,0);
    k.resize(N,0);  
    X.resize(N*2,0);
}
