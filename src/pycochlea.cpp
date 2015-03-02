#include "stdio.h"
#include <boost/numeric/odeint.hpp>
#include "cochlea.hpp"


using namespace std;
using namespace boost::numeric::odeint;

//extern void _main();
//------------------------------------------------------------------

extern "C" void cochlea(flotante *X_t, flotante *tt, flotante *stimulus, flotante* weight, flotante* sdivm, flotante* ddivm, flotante* params, int* dimension, flotante* solver_opts ) 
{

//params dv, dvy2, alpha, beta
//flotante dv, dvy2, alpha, beta
//dimensions
//int n_t, n_osc, fs
//int solver, flotante abs_tol = 1e-4,flotante rel_tol = 1e-3

    int solver = (int)solver_opts[0];
    flotante abs_tol = solver_opts[1];
    flotante rel_tol = solver_opts[1];   

    Cochlea_t C;
 
    C.alpha = params[0];
    C.beta = params[1];
    C.gamma = params[2];
    C.delta = params[3];
    C.eps = params[4];
    C.zeta = params[5];
    C.eta = params[6];
    C.nu = params[7];
    C.m0 = params[8];


//    cout << C.alpha << "\t" << C.beta << "\t"<< C.gamma << "\t"<< C.delta << "\t"<< C.eps << "\t"<< C.zet << "\t"<< C.eta << "\t" << C.fluid << endl;

    C.n_t = dimension[0];
    C.n_osc = dimension[1];
    int fs = dimension[2];
    int dec = dimension[3];
     
    C.stimulus = stimulus;
    C.weight = weight;
    C.ddivm = ddivm;
    C.sdivm = sdivm;
    
    flotante dt = 1.0/(flotante) fs;  
    flotante tmax = ( C.n_t-1)*dt;

    C.dt = dt;
    C.setup(); 
    
    vector< flotante > X;

    X.resize(C.n_osc*2,0);

    //Runge Kutta constant dt 
    if(solver==0)
    {
        Observer_t obs(X_t,tt,C.n_osc*2,dec);

        runge_kutta4< flotante_vector > stepper;
        integrate_const( stepper , C , X , 0.0, tmax , dt ,obs );
    }
    //Runge Kutta Dopri dense ouput
    if(solver==1)
    {
        Observer_t obs(X_t,tt,C.n_osc*2);

        size_t steps = integrate_const( make_dense_output< runge_kutta_dopri5< flotante_vector > >( abs_tol, rel_tol ),
                        C , X , 0.0 , tmax , dt , obs );

        cout << steps << endl;
    }
    //Runge Kutta Cash Karp54 dense ouput
    if(solver==2)
    {
        Observer_t obs(X_t,tt,C.n_osc*2);
    
        size_t steps = integrate_const( make_controlled< runge_kutta_cash_karp54< flotante_vector >  >( abs_tol , rel_tol) , 
                        C , X , 0.0 , tmax , dt, obs );

        cout << steps << endl;
    }
    
}
