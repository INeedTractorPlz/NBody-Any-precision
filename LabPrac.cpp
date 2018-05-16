#include <gmpxx.h>
#include<gmp.h>
#include<mpfr.h>
#include <variant>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include<boost/lexical_cast.hpp>

#include <unistd.h>
#include <getopt.h>

#include<ios>
#include<fstream>
#include <functional>
#include<utility>
#include<string>
#include<iostream>


using namespace boost::numeric::ublas;
using namespace boost::numeric::odeint;

#define NumberOfIntegrators 3


typedef mpf_class data_type;
typedef matrix<data_type> state_type;
typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
typedef std::pair<state_type,state_type> symplectic_type;

template<typename type, typename Type>
struct RungeKutta4;
template<typename type, typename Type>
struct RungeKutta5_Fehlberg;
template<typename type, typename Type>
struct LeapFrog;


typedef std::variant<RungeKutta4<data_type, state_type> , RungeKutta5_Fehlberg<data_type, state_type>,
    LeapFrog<data_type, symplectic_type> > choosen_integrator_t;
mpf_class atan(mpf_class in,mpfr_rnd_t rnd=MPFR_RNDN){
    mpfr_t out_t;
    mpf_class out;
    mpfr_init_set_f(out_t,in.get_mpf_t(),rnd);
    mpfr_atan(out_t,out_t,rnd);
    mpfr_get_f(out.get_mpf_t(),out_t,rnd);
    return out;
}
mpf_class mpf_class_set_str(const char *s, mpfr_prec_t prec = 64, mpfr_rnd_t rnd=MPFR_RNDN){
    mpfr_t out_t;
    mpfr_init2(out_t, prec);
    mpf_class out;
    mpfr_set_str(out_t, s, 10, rnd);
    mpfr_get_f(out.get_mpf_t(),out_t,rnd);
    return out;
}

struct NBody{
    std::vector< state_type >& m_states;
    std::vector< data_type >& m_times, &energy_result;
    const std::vector<data_type>& m;
    std::vector<state_type >& norm_R;
    data_type G;
    const data_type precision_energy;
    unsigned &current;
    unsigned p,dim;
    NBody(std::vector< state_type >& m_states, std::vector< data_type >& m_times, 
    std::vector<data_type>& energy_result, const data_type &precision_energy,
    const std::vector<data_type>& m,unsigned &current,const unsigned &dim,
    std::vector<state_type >& norm_R):  m_states(m_states), m_times(m_times),
    m(m), energy_result(energy_result), precision_energy(precision_energy), 
    current(current), dim(dim), norm_R(norm_R),p(m.size()){
        data_type pi=atan(data_type(1.))*4;
        G=4*pi*pi;
    }
    void operator()(const state_type& R, state_type& A,  data_type t){
        A=zero_matrix<data_type>(p,dim);
        subrange(A,0,p,0,dim/2)=subrange(R,0,p,dim/2,dim);
        state_type norm_r(p,p);

        norm(R,norm_r);        
        for(unsigned j=0;j<p;++j){
            for(unsigned k=0;k<j;++k){
                subrange(A,j,j+1,dim/2,dim)+=m[k]*(subrange(R,k,k+1,0,dim/2)-subrange(R,j,j+1,0,dim/2))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
            for(unsigned k=j+1;k<p;++k){
                subrange(A,j,j+1,dim/2,dim)+=m[k]*(subrange(R,k,k+1,0,dim/2)-subrange(R,j,j+1,0,dim/2))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
        }
        subrange(A,0,p,dim/2,dim)*=G;
    }    
    void operator()( const state_type &R , data_type t ){
        if(t!=0) ++current;
        m_states.push_back(R); 
        m_times.push_back(t);
        norm_R.push_back(std::move(state_type(p,p))); 
        norm(R,norm_R[current]);
        energy_result.push_back(EnergyIntegral(R,norm_R[current]));
    }
    bool operator()(const state_type& X, const state_type& Y){
        state_type norm_Y=std::move(state_type(p,p));
        norm(Y,norm_Y);
        data_type energy_Y(EnergyIntegral(Y,norm_Y));
        data_type delta=abs((energy_result[0]-energy_Y)/energy_result[0]);
        //std::cout << "delta= " << delta << std::endl;
        if(delta > precision_energy) return 1;
        else return 0;
    }
    void norm(const state_type& R, state_type& norm_r);
    data_type EnergyIntegral(const state_type& R, const state_type& norm_r);
};

struct SymplecticNBody : NBody{
    SymplecticNBody(std::vector< state_type >& m_states, std::vector< data_type >& m_times, 
    std::vector<data_type>& energy_result, const data_type &precision_energy,
    const std::vector<data_type>& m,unsigned &current,const unsigned &dim,
    std::vector<state_type >& norm_R):  NBody(m_states, m_times, energy_result, precision_energy,
    m, current, dim, norm_R){  }
    
    void operator()(const state_type& R, state_type& A,  data_type t){
        A=zero_matrix<data_type>(p,dim/2);
        state_type norm_r(p,p);

        norm(R,norm_r);
        for(unsigned j=0;j<p;++j){
            for(unsigned k=0;k<j;++k){
                row(A,j)+=m[k]*(row(R,k)-row(R,j))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
            for(unsigned k=j+1;k<p;++k){
                row(A,j)+=m[k]*(row(R,k)-row(R,j))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
        }
        A*=G;
    }    
    void operator()( const symplectic_type &R , data_type t ){
        if(t!=0) ++current;
        state_type R_matrix(p,dim);
        subrange(R_matrix,0,p,0,dim/2)=R.first;
        subrange(R_matrix,0,p,dim/2,dim)=R.second;
        
        m_states.push_back(R_matrix); 
        m_times.push_back(t);
        norm_R.push_back(std::move(state_type(p,p))); 
        norm(R.first,norm_R[current]);
        energy_result.push_back(EnergyIntegral(R_matrix,norm_R[current]));
    }
    bool operator()(const symplectic_type& X, const symplectic_type& Y){
        state_type norm_Y=std::move(state_type(p,p));
        state_type Y_matrix(p,dim);
        subrange(Y_matrix,0,p,0,dim/2)=Y.first;
        subrange(Y_matrix,0,p,dim/2,dim)=Y.second;
    
        norm(Y_matrix,norm_Y);
        data_type energy_Y(EnergyIntegral(Y_matrix,norm_Y));
        data_type delta=abs((energy_result[0]-energy_Y)/energy_result[0]);
        //std::cout << "delta= " << delta << std::endl;
        if(delta > precision_energy) return 1;
        else return 0;
    }
};

data_type NBody::EnergyIntegral(const state_type& R, const state_type& norm_r){
    data_type E=0,U;
    for(unsigned i=0;i<m.size();++i){
        U=0;
        for(unsigned j=i+1;j<m.size();++j)
            U-=m[j]*m[i]/(norm_r(i,j));
        E+=U*G;
        E+=m[i]*inner_prod(row(subrange(R,0,p,dim/2,dim),i),row(subrange(R,0,p,dim/2,dim),i))/2;
        //E+=m[i]*(R(i,3)*R(i,3)+R(i,4)*R(i,4)+R(i,5)*R(i,5))/2;
    }    
    return E;
}

template<typename type>
struct Arguments_t{
    unsigned &bits, &number_bodies, &number_steps, &dim;
    type &time_end, &precision_energy, &distance_encounter;
    std::string &name_integrator;
    bool &without_controlled_step;
    Arguments_t(unsigned &number_bodies,type &time_end,unsigned &number_steps,unsigned &bits,
    type &precision_energy, std::string &name_integrator, unsigned &dim,
    type &distance_encounter, bool &without_controlled_step) : 
    bits(bits), number_bodies(number_bodies),time_end(time_end),
    number_steps(number_steps),precision_energy(precision_energy),
    name_integrator(name_integrator), dim(dim), distance_encounter(distance_encounter),
    without_controlled_step(without_controlled_step){ }
    void Parser(int argc, char *const argv[]){
        int i, longIndex=0;
        std::string ss;
        static const char* shortopts = "b:M:t:N:p:";
        static const struct option longOpts[] ={
            {"integrator",required_argument, NULL,0},
            {"without_controlled_step",no_argument, NULL, 0}  
        };
        auto opt = getopt_long(argc,argv,shortopts,longOpts,&longIndex);
        while(opt != -1){
            switch(opt){
                case 'b':
                    ss=optarg;
                    i = boost::lexical_cast<int>(ss);
                    if(i < 64){
                        std::cout << "Don't possible use less 64 bits";
                        exit(-1);
                    }else
                        bits = i;
                    break;
                case 'M':
                    ss=optarg;
                    i = boost::lexical_cast<int>(ss);
                    if(i < 2){
                        std::cout << "Don't possible use less two bodies";
                        exit(-1);
                    }else
                        number_bodies = i;
                    break;
                case 't':
                    if(typeid(type) == typeid(mpf_class))
                        time_end = mpf_class_set_str(optarg,bits);
                    else
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            time_end = boost::lexical_cast<type>(ss);
                        }
                        catch(...){
                            std::cout << "String can't be converted to this type(" << 
                            typeid(type).name() << ")" << std::endl;
                            exit(-1);
                        }
                    break;
                case 'N':
                    ss=optarg;
                    i = boost::lexical_cast<int>(ss);
                    if(i < 1){
                        std::cout << "Are you sure, that want divide a segment on a " << i << "parts?";
                        std::cout << std::endl;
                        exit(-1); 
                    }else
                        number_steps = i;
                    break;
                case 'p':
                    if(typeid(type) == typeid(mpf_class))
                        precision_energy = mpf_class_set_str(optarg,bits);
                    else
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            precision_energy = boost::lexical_cast<type>(ss);
                        }
                        catch(...){
                            std::cout << "String can't be converted to this type(" << 
                            typeid(type).name() << ")" << std::endl;
                            exit(-1);
                        }
                    break;
                case 'd':
                    ss=optarg;
                    i = boost::lexical_cast<int>(ss);
                    if(i < 1){
                        std::cout << "Dimension" << i << "isn't exists";
                        std::cout << std::endl;
                        exit(-1); 
                    }else
                        dim = 2*i;
                    break;
                case 'e':
                    if(typeid(type) == typeid(mpf_class))
                        precision_energy = mpf_class_set_str(optarg,bits);
                    else
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            distance_encounter = boost::lexical_cast<type>(ss);
                        }
                        catch(...){
                            std::cout << "String can't be converted to this type(" << 
                            typeid(type).name() << ")" << std::endl;
                            exit(-1);
                        }
                    break;
                case 0:
                    if("integrator" == longOpts[longIndex].name){
                        name_integrator = optarg;
                    }
                    if("without_controlled_step" == longOpts[longIndex].name){
                        without_controlled_step = 1;
                    }
                    break;
                default :
                    {
                        std::cout << "This key isn't exists" << std::endl;
                        std::cout << std::endl;
                        exit(-1);
                    }
            }  
            opt = getopt_long(argc,argv,shortopts,longOpts,&longIndex);     
        }
    }
};



void NBody::norm(const state_type& R, state_type& norm_r){
   //data_type r1, r2, r3;
   for(unsigned i=0;i<R.size1();++i){
        for(unsigned j=0;j<i;++j)
            norm_r(i,j)=norm_r(j,i);
         for(unsigned j=i+1;j<R.size1();++j){
            norm_r(i,j)=norm_2(row(subrange(R,0,p,0,dim/2),i)-row(subrange(R,0,p,0,dim/2),j));
            //r1=R(i,0)-R(j,0); r2=R(i,1)-R(j,1); 
            //r3=R(i,2)-R(j,2);
            //norm_r(i,j)=sqrt(r1*r1+r2*r2+r3*r3);
        }
    }
}

const char* filenamestr(const char* s1, unsigned j, const char* s2){
    static std::string ss;
    ss="";
    ss+= s1; ss+=boost::lexical_cast<std::string>(j); ss+=s2;
    return ss.c_str();
}


template<typename type, typename Type>
struct RungeKutta4{
    Type k1,k2,k3,k4;
    RungeKutta4(){}
    template<typename sysf>
    void do_step(sysf sysF, const Type &in, Type &out, type time, type step){
        sysF(in,k1,time);
        sysF(in+step*k1/2.,k2,time+step/2.);
        sysF(in+step*k2/2.,k3,time+step/2.);
        sysF(in+step*k3,k4,time+step);
        out=std::move(in + step*(k1+2.*k2+2.*k3+k4)/6.);
    }
};

template<typename type, typename Type>
struct RungeKutta5_Fehlberg{
    Type k1,k2,k3,k4,k5,k6;
    RungeKutta5_Fehlberg(){}
    template<typename sysf>
    void do_step(sysf sysF, const Type &in, Type &out, type time, type step){
        sysF(in,k1,time);
        sysF(in+step*k1/4.,k2,time+step/4.);
        sysF(in+step*(3*k1/32.+9*k2/32.),k3,time+3*step/8);
        sysF(in+step*(1932*k1/2197-7200*k2/2197+7296*k3/2197),k4,time+12*step/13);
        sysF(in+step*(439*k1/216-8*k2+3680*k3/513-845*k4/4104),k5,time+step);
        sysF(in+step*(-8*k1/27+ 2*k2-3544*k3/2565+1859*k4/4104-11*k5/40),k6,time+step/2);
        
        out=std::move(in + step*(16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55));
    }
};

template<typename type, typename Type>
struct LeapFrog{
    Type intermediate;
    bool first_in=1;
    LeapFrog(){}
    template<typename sysf>
    void do_step(sysf sysF, const Type &in, Type &out, type time, type step){
        if(first_in){
            sysF(in.first,intermediate.second,time);
            first_in=0;
        }
        intermediate.first=in.second+intermediate.second*step/2;
        out.first=in.first+intermediate.first*step;
        sysF(out.first,intermediate.second,time+step);
        out.second=intermediate.first+intermediate.second*step/2;
    }
};

template<typename Type, typename type>
struct Encounter{
    type distance_encounter;
    Encounter(const type& in) : distance_encounter(in){}
    bool operator()(const Type& X, const type& t){
        if(norm_2(subrange(row(X,0) - row(X,1),0,3)) <= distance_encounter){
            std::cout << "ENCOUNTER!! in time " << t << std::endl;
            return 0;
        }else{
            std::cout << norm_2(subrange(row(X,0) - row(X,1),0,3)) << "   " << t << std::endl;
            return 1;
        }
    } 
};

template<typename Type, typename type>
struct SymplecticEncounter{
    type distance_encounter;
    SymplecticEncounter(const type& in) : distance_encounter(in){}
    bool operator()(const Type& X, const type& t){
        if(norm_2(row(X.first,0) - row(X.first,1)) <= distance_encounter){
            std::cout << "ENCOUNTER!! in time " << t << std::endl;
            return 0;
        }else{
            std::cout << norm_2(row(X.first,0) - row(X.first,1)) << "   " << t << std::endl;
            return 1;
        }
    } 
};

template<typename type, typename Integrator, typename sysf, typename observer,typename Type, 
typename controller=std::function<int(const Type& X, const Type& Y)>, 
typename exit=std::function<bool(const Type&X, const type &t)> >
void integrate(Integrator rk, sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){
    int state;
    Type Y;
    type T=t+h*n;
    Obs(X0,t);
    while(t<T && Bad(X0,t)){
        state=-1;
        rk.do_step(sysF,X0,Y,t,h);
        if(Err(X0,Y)==-1){
            goto obs_point;
        } 
        while(Err(X0,Y)){
            h=h/2.;
            rk.do_step(sysF, X0,Y,t,h);
            std::cout << "h= " << h << std::endl;
            //std::cout << "Y :" << std::endl << Y << std::endl;
            state=1;
        }
        if(state==1){
            std::cout << "STOP1" << std::endl;
            h=h*2;
            rk.do_step(sysF,X0,Y,t,h);
            goto obs_point;
        }
        while(!Err(X0,Y)){
            std::cout << "STOP" << std::endl;
            h=h*2;
            rk.do_step(sysF, X0,Y,t,h);
            state=0;
        }
        if(state==0){
            h=h/2;
            rk.do_step(sysF,X0,Y,t,h);
        }
        obs_point:
        X0=Y;
        t+=h;
        Obs(X0,t);
    }
}
/*template<int N, int... Rest>
struct SetOfIndex_t{
    static constexpr auto& value = SetOfIndex_t<N - 1, N, Rest...>::value;
};
template<int... Rest>
struct SetOfIndex_t<0, Rest...>{
    static constexpr int value[] = {0, Rest...};
};
*/
template<typename data_type>
struct Integrate_based_args_t{
    Arguments_t<data_type> args;
    std::string integrator_type;
    choosen_integrator_t choosen_integrator;
    Integrate_based_args_t(const Arguments_t<data_type> &args) : args(args) { }
    std::string choise_integrator(){
        if (args.name_integrator == "RK4"){
            choosen_integrator = RungeKutta4<data_type, state_type>();
            integrator_type = "standart";
            return integrator_type;
        }
        if (args.name_integrator == "RK5"){
            choosen_integrator = RungeKutta5_Fehlberg<data_type, state_type>();
            integrator_type = "standart";
            return integrator_type;
        }
        if (args.name_integrator == "LP"){
            choosen_integrator = LeapFrog<data_type, symplectic_type>();
            integrator_type = "symplectic";
            return integrator_type;
        }
        throw;
    } 
    template<typename type, typename Integrator, typename sysf, typename observer,typename Type, 
    typename controller=std::function<int(const Type& X, const Type& Y)>, 
    typename exit=std::function<bool(const Type&X, const type &t)> >
    void integrator_call(Integrator rk, sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){
        if(args.distance_encounter == 0){
            if(args.without_controlled_step)
                integrate(rk, sysF, X0, t, h, n, Obs);
            else
                integrate(rk, sysF, X0, t, h, n, Obs,[](const Type& X, const type& t){return 1;},Err);
        }else{
            if(args.without_controlled_step)
                integrate(rk, sysF, X0, t, h, n, Obs,Bad);
            else
                integrate(rk, sysF, X0, t, h, n, Obs,Bad, Obs);
        }
    }
    template<typename type, typename sysf, typename observer,typename Type, 
    typename controller=std::function<int(const Type& X, const Type& Y)>, 
    typename exit=std::function<bool(const Type&X, const type &t)>>
    void integrator_call_call(sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){
                    //static constexpr auto& SetOfIndex = SetOfIndex_t<NumberOfIntegrators>::value;
                    if(auto rk = std::get_if<0>(&choosen_integrator))
                        integrator_call(std::get<0>(choosen_integrator),sysF,X0,t,h,n,Obs,Bad,Err);
                    else{
                            if(auto rk = std::get_if<1>(&choosen_integrator))
                                integrator_call(std::get<1>(choosen_integrator),sysF,X0,t,h,n,Obs,Bad,Err);
                        else{
                            if(auto rk = std::get_if<2>(&choosen_integrator)){}
                                //integrator_call(std::get<2>(choosen_integrator),sysF,X0,t,h,n,Obs,Bad,Err);
                            else
                                std::cout << "INTEGRATOR CALL ERROR";
                            }
                    }
                    /*try{
                        integrator_call(std::get<0>(choosen_integrator),sysF,X0,t,h,n,Obs,Bad,Err);
                    }catch(...){
                        try{
                            integrator_call(std::get<1>(choosen_integrator),sysF,X0,t,h,n,Obs,Bad,Err);
                        }catch(...){
                            //integrator_call(std::get<2>(choosen_integrator),sysF,X0,t,h,n,Obs,Bad,Err);
                        }
                    }*/   
    }
};


int main(int argc, char *const argv[]){
    std::ios_base::sync_with_stdio(false);
    
    std::ofstream out,energy;
    std::ifstream file_in,file_size,file_mass,file_initial,file_RBinitial;
    std::vector<data_type> m;
    std::vector<state_type> state_result;
    std::vector<data_type> time_result, energy_result;
    std::vector<state_type >norm_R;
    unsigned current=0;
    unsigned N,M, bits = 64, dim = 6;
    data_type T,t0(0.), max_delta=0, delta, distance_encounter = 0, precision_energy = 0;
    state_type X0,Y;
    symplectic_type Z;
    std::string name_integrator = "RK4";
    std::string type_integrator = "standart";
    bool without_controlled_step = 0;
    //runge_kutta4<state_type> rk;
    //controlled_stepper_type controlled_stepper;
    
    
    file_size.open("file_size.dat",std::ios_base::in);
    file_size >> T >> N >> M >> dim >> bits >> distance_encounter;
    file_size.close();

    Arguments_t<data_type> keys(N,T,M,bits,precision_energy,name_integrator,dim,
    distance_encounter,without_controlled_step);
    keys.Parser(argc,argv);

    std::cout.precision(bits);
    out.precision(bits);
    energy.precision(bits);
    mpf_set_default_prec(bits);
    std::cout << mpf_get_default_prec() << std::endl;
    

    
    m.resize(M);
    X0.resize(M,dim);
    
    file_mass.open("file_mass.dat",std::ios_base::in);
    for(unsigned i=0;i<M;++i){
        file_mass >> m[i];
    }    
    file_mass.close();
    
    file_initial.open("file_initial.dat",std::ios_base::in);
    for(unsigned i=0;i<M;++i){
        for(int j=0;j<dim;++j)
            file_initial >> X0(i,j);
    }        
    file_initial.close();

    
    Y=X0;
    std::cout << Y << std::endl;
    std::cout << name_integrator << std::endl;
    Integrate_based_args_t Integrate_based_args(keys);
    try{
        type_integrator = Integrate_based_args.choise_integrator();
    }
    catch(...){
        std::cout << "This integrator isn't exists" << std::endl;
        exit(-1);
    }
    NBody nb(state_result, time_result, energy_result, precision_energy, m,current,dim,norm_R);
    SymplecticNBody snb(state_result, time_result, energy_result,  precision_energy, m,current,dim,norm_R);
        
    if(type_integrator == "standart"){
        Integrate_based_args.integrator_call_call(nb , Y , t0, data_type(T/N), N, nb,
                Encounter<state_type, data_type>(distance_encounter), nb);
    }
    /*if(type_integrator == "symplectic"){
        Z.first=subrange(Y,0,M,0,dim/2);
        Z.second=subrange(Y,0,M,dim/2,dim);
        Integrate_based_args.integrator_call_call(snb , Z , t0, data_type(T/N), N, snb,
                SymplecticEncounter<symplectic_type, data_type>(distance_encounter), snb);
    }*/
    auto index_integrator = Integrate_based_args.choosen_integrator.index();
    std::cout << index_integrator << std::endl;
    std::cout << typeid(Integrate_based_args.choosen_integrator).name() << std::endl;
    
    //size_t steps=integrate_n_steps( rk , nb , Y , 0., T/N, N, nb);
    //size_t steps=integrate_adaptive( make_controlled< error_stepper_type >( 1.0e-12 , 1.0e-8 ),
    //nb, Y, 0., T, T/N, nb);
    //std::cout << "steps=" << steps << std::endl;
    //std::cout << Y << std::endl;

    /*
    out.open("res.dat",std::ios_base::trunc);
    for(unsigned i=0;i<=current;++i){
        for(unsigned j=0;j<M;++j){
            for(unsigned k=0;k<dim;++k)
                out << state_result[i](j,k) << " ";
            out << std::endl;
        }
        out << std::endl << std::endl;
    }
    out.close();
      for(unsigned j=0;j<M;++j){
        out.open(filenamestr("NBodyRK_",j,".dat"),std::ios_base::trunc);
        std::cout << filenamestr("NBodyRK_",j,".dat") << std::endl;
        for(unsigned i=0;i<=current;++i){
            out << time_result[i] << " ";
            out << row(state_result[i],j) << std::endl;
        }        
        out.close();
    }
    
    energy.open("NBodyRK_energy.dat",std::ios_base::trunc);
    for(unsigned i=0;i<=current;++i){
        energy <<  "h(" << i << ")= " << energy_result[i] << std::endl;
        delta=abs((energy_result[0]-energy_result[i])/energy_result[i]);
        if(delta > max_delta) max_delta=delta;
    }
    std::cout << max_delta << std::endl;
    */
   return 0;
}