#include <gmpxx.h>
#include<gmp.h>
#include<mpfr.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include<boost/lexical_cast.hpp>

#include <unistd.h>
#include <getopt.h>

#include <algorithm>
#include<ios>
#include<fstream>
#include <functional>
#include<utility>
#include<string>
#include<iostream>

#define DOUBLE_PRESICION

#ifdef ANY_PRESICION
#define ABS_ abs
typedef mpf_class data_type;
#endif

#ifdef DOUBLE_PRESICION
#define ABS_ fabs
typedef double data_type;
#endif


#ifdef FLOAT_PRESICION
#define ABS_ fabs
typedef float data_type;
#endif

using namespace boost::numeric::ublas;
//using namespace boost::numeric::odeint;

typedef matrix<data_type> state_type;
typedef vector<data_type> state_vector;
//typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
//typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
typedef std::pair<state_type,state_type> symplectic_type;

template<typename type, typename Type>
struct RungeKutta4;
template<typename type, typename Type>
struct RungeKutta5_Fehlberg;
template<typename type, typename Type>
struct LeapFrog;

template<typename type>
std::ostream& operator<<(std::ostream& os, const std::pair<type,type> in){
    os << in.first << std::endl;
    os << in.second << std::endl;
    return os;
} 

mpf_class atan(mpf_class in,mpfr_rnd_t rnd=MPFR_RNDN){
    mpfr_t out_t;
    mpf_class out;
    mpfr_init_set_f(out_t,in.get_mpf_t(),rnd);
    mpfr_atan(out_t,out_t,rnd);
    mpfr_get_f(out.get_mpf_t(),out_t,rnd);
    return out;
}

mpf_class acos(mpf_class in,mpfr_rnd_t rnd=MPFR_RNDN){
    mpfr_t out_t;
    mpf_class out;
    mpfr_init_set_f(out_t,in.get_mpf_t(),rnd);
    mpfr_acos(out_t,out_t,rnd);
    mpfr_get_f(out.get_mpf_t(),out_t,rnd);
    return out;
}

mpf_class cos(mpf_class in,mpfr_rnd_t rnd=MPFR_RNDN){
    mpfr_t out_t;
    mpf_class out;
    mpfr_init_set_f(out_t,in.get_mpf_t(),rnd);
    mpfr_cos(out_t,out_t,rnd);
    mpfr_get_f(out.get_mpf_t(),out_t,rnd);
    return out;
}

mpf_class atan2(mpf_class cos_x, mpf_class sin_x, mpfr_rnd_t rnd=MPFR_RNDN){
    mpfr_t cos_t,sin_t, out_t;
    mpf_class out;
    mpfr_init_set_f(cos_t,cos_x.get_mpf_t(),rnd);
    mpfr_init_set_f(sin_t,sin_x.get_mpf_t(),rnd);
    mpfr_init2(out_t,cos_x.get_prec());
    mpfr_atan2(out_t,sin_t,cos_t,rnd);
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

state_vector cross_product(const state_vector &u, const state_vector &v){
    state_vector result(u.size());
    result(0) = u(1)*v(2) - u(2)*v(1);
    result(1) = u(2)*v(0) - u(0)*v(2);
    result(2) = u(0)*v(1) - u(1)*v(0);
    return result;
}

struct NBody{
    std::ofstream res_file,energy_file,orbital_file;
    std::string name_res,name_energy,name_orbital;
    std::vector< state_type >& m_states;
    std::vector< data_type >& m_times, &energy_result;
    std::vector<std::vector<state_vector> > &orbital_elements_result; 
    const std::vector<data_type>& m;
    std::vector<state_type >& norm_R;
    std::vector<std::pair<unsigned, unsigned> > pairs_of_bodies;
    unsigned number_orbital_elements;
    data_type &G,pi = atan(data_type(1.0))*4;
    const data_type precision_energy;
    unsigned &current;
    unsigned number_bodies,dim;
    NBody(std::vector< state_type >& m_states, std::vector< data_type >& m_times, 
    std::vector<data_type>& energy_result, const data_type &precision_energy,
    const std::vector<data_type>& m,unsigned &current,const unsigned &dim,
    std::vector<state_type >& norm_R, data_type &G,
    std::vector<std::vector<state_vector> > &orbital_elements_result,
    std::vector<std::pair<unsigned, unsigned> > pairs_of_bodies,
    unsigned number_orbital_elements, const char *name_res, const char *name_energy,
    const char *name_orbital):m_states(m_states),
    m_times(m_times),m(m), energy_result(energy_result),precision_energy(precision_energy),
    current(current),dim(dim), norm_R(norm_R),number_bodies(m.size()), G(G),
    orbital_elements_result(orbital_elements_result),
    pairs_of_bodies(pairs_of_bodies),number_orbital_elements(number_orbital_elements),
    name_res(name_res),name_energy(name_energy),name_orbital(name_orbital){    }

    void operator()(const state_type& R, state_type& A,  data_type t){
        A=zero_matrix<data_type>(number_bodies,dim);
        subrange(A,0,number_bodies,0,dim/2)=subrange(R,0,number_bodies,dim/2,dim);
        state_type norm_r(number_bodies,number_bodies);

        norm(R,norm_r);        
        for(unsigned j=0;j<number_bodies;++j){
            for(unsigned k=0;k<j;++k){
                subrange(A,j,j+1,dim/2,dim)+=m[k]*(subrange(R,k,k+1,0,dim/2)-subrange(R,j,j+1,0,dim/2))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
            for(unsigned k=j+1;k<number_bodies;++k){
                subrange(A,j,j+1,dim/2,dim)+=m[k]*(subrange(R,k,k+1,0,dim/2)-subrange(R,j,j+1,0,dim/2))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
        }
        subrange(A,0,number_bodies,dim/2,dim)*=G;
    }    
    void operator()( const state_type &R , data_type t ){
        if(t!=0) ++current;
        result_record(R,t);
    }

    bool operator()(const state_type& X, const state_type& Y){
        state_type norm_Y=std::move(state_type(number_bodies,number_bodies));
        norm(Y,norm_Y);
        data_type energy_Y(EnergyIntegral(Y,norm_Y));
        data_type delta=ABS_((energy_result[current]-energy_Y)/energy_result[current]);
        //std::cout << "delta= " << delta << std::endl;
        if(delta > precision_energy) return 1;
        else return 0;
    }
    void norm(const state_type& R, state_type& norm_r);
    data_type EnergyIntegral(const state_type& R, const state_type& norm_r);
    void get_orbital_elements(const state_vector &R, unsigned i, unsigned j,
    state_vector& orbital_elements);
    void result_record(const state_type &R, data_type t);
    void output_to_file(const state_type &R, data_type t);
};

struct SymplecticNBody : NBody{
    SymplecticNBody(std::vector< state_type >& m_states, std::vector< data_type >& m_times, 
    std::vector<data_type>& energy_result, const data_type &precision_energy,
    const std::vector<data_type>& m,unsigned &current,const unsigned &dim,
    std::vector<state_type >& norm_R, data_type &G,
    std::vector<std::vector<state_vector> > &orbital_elements_result,
    std::vector<std::pair<unsigned,unsigned> > pairs_of_bodies,
    unsigned number_orbital_elements,const char *name_res, const char *name_energy,
    const char *name_orbital):
    NBody(m_states, m_times, energy_result, precision_energy,m,
    current, dim, norm_R, G,orbital_elements_result,pairs_of_bodies,
    number_orbital_elements,name_res,name_energy,name_orbital){  }
    
    void operator()(const state_type& R, state_type& A,  data_type t){
        A=zero_matrix<data_type>(number_bodies,dim/2);
        state_type norm_r(number_bodies,number_bodies);

        norm(R,norm_r);
        for(unsigned j=0;j<number_bodies;++j){
            for(unsigned k=0;k<j;++k){
                row(A,j)+=m[k]*(row(R,k)-row(R,j))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
            for(unsigned k=j+1;k<number_bodies;++k){
                row(A,j)+=m[k]*(row(R,k)-row(R,j))/
                (norm_r(k,j)*norm_r(k,j)*norm_r(k,j));
            }
        }
        A*=G;
    }    
    void operator()( const symplectic_type &R , data_type t ){
        if(t!=0) ++current;
        state_type R_matrix(number_bodies,dim);
        subrange(R_matrix,0,number_bodies,0,dim/2)=R.first;
        subrange(R_matrix,0,number_bodies,dim/2,dim)=R.second;
        
        m_states.push_back(R_matrix); 
        m_times.push_back(t);
        norm_R.push_back(std::move(state_type(number_bodies,number_bodies))); 
        norm(R.first,norm_R[current]);
        energy_result.push_back(EnergyIntegral(R_matrix,norm_R[current]));
        for(unsigned i=0;i<number_orbital_elements;++i){
            state_vector distance_i = row(R_matrix,pairs_of_bodies[i].first)
            -row(R_matrix,pairs_of_bodies[i].second);
            std::for_each(distance_i.begin(),distance_i.end(),[](data_type &x){x=ABS_(x);});
            orbital_elements_result[i].push_back(std::move(state_vector(4)));
            get_orbital_elements(distance_i,pairs_of_bodies[i].first,pairs_of_bodies[i].second,
            orbital_elements_result[i][current]);
        }
    }
    bool operator()(const symplectic_type& X, const symplectic_type& Y){
        state_type norm_Y=std::move(state_type(number_bodies,number_bodies));
        state_type Y_matrix(number_bodies,dim);
        subrange(Y_matrix,0,number_bodies,0,dim/2)=Y.first;
        subrange(Y_matrix,0,number_bodies,dim/2,dim)=Y.second;
        
        norm(Y_matrix,norm_Y);
        data_type energy_Y(EnergyIntegral(Y_matrix,norm_Y));
        data_type delta=ABS_((energy_result[current]-energy_Y)/energy_result[current]);
        //std::cout << "ABS_(energy_result[0]-energy_Y)= " << ABS_(energy_result[0]-energy_Y) << std::endl;
        //std::cout << "energy_result[0]= " << energy_result[0] << std::endl;
        //std::cout << "energy_Y= " << energy_Y << std::endl;
        //std::cout << "delta= " << delta << std::endl;
        if(delta > precision_energy) return 1;
        else return 0;
    }
};

void NBody::result_record(const state_type &R, data_type t){
    m_states.push_back(R); 
    m_times.push_back(t);
    norm_R.push_back(std::move(state_type(number_bodies,number_bodies))); 
    norm(R,norm_R[current]);
    energy_result.push_back(EnergyIntegral(R,norm_R[current]));
    for(unsigned i=0;i<number_orbital_elements;++i){
        state_vector distance_i = row(R,pairs_of_bodies[i].first)-row(R,pairs_of_bodies[i].second);
        std::for_each(distance_i.begin(),distance_i.end(),[](data_type &x){x=ABS_(x);});
        orbital_elements_result[i].push_back(std::move(state_vector(4)));
        get_orbital_elements(distance_i,pairs_of_bodies[i].first,pairs_of_bodies[i].second,
        orbital_elements_result[i][current]);
    }
}

void NBody::get_orbital_elements(const state_vector &R, unsigned i, unsigned j,
state_vector& orbital_elements){
    data_type kappa_quad,r,v;
    
    kappa_quad=G*(m[i]+m[j]);
    r=norm_R[current](i,j);
    v=norm_2(subrange(R,dim/2,dim));
    
    data_type h = v*v/2 - kappa_quad/r;
    state_vector c = cross_product(subrange(R,0,dim/2),subrange(R,dim/2,dim));
    
    data_type a,e,inclination;
    e = sqrt(1.+2*h*norm_2(c)*norm_2(c)/(kappa_quad*kappa_quad));
    a = -kappa_quad/2/h;
    
    inclination = acos(c(2)/norm_2(c));
    
    orbital_elements(0) = a;
    orbital_elements(1) = e;
    orbital_elements(2) = inclination*180/pi;
    orbital_elements(3) = (1-e)*a;
}

data_type NBody::EnergyIntegral(const state_type& R, const state_type& norm_r){
    data_type E=0,U;
    for(unsigned i=0;i<m.size();++i){
        U=0;
        for(unsigned j=i+1;j<m.size();++j)
            U-=m[j]*m[i]/(norm_r(i,j));
        E+=U*G;
        E+=m[i]*inner_prod(row(subrange(R,0,number_bodies,dim/2,dim),i),
        row(subrange(R,0,number_bodies,dim/2,dim),i))/2;
    }    
    return E;
}


void NBody::norm(const state_type& R, state_type& norm_r){
   for(unsigned i=0;i<R.size1();++i){
        for(unsigned j=0;j<i;++j)
            norm_r(i,j)=norm_r(j,i);
         for(unsigned j=i+1;j<R.size1();++j){
            norm_r(i,j)=norm_2(row(subrange(R,0,number_bodies,0,dim/2),i)
            -row(subrange(R,0,number_bodies,0,dim/2),j));
        }
    }
}

template<typename type>
struct Arguments_t{
    unsigned &bits, &number_bodies, &number_steps, &dim,&number_orbital_elements;
    type &time_end, &precision_energy, &distance_encounter, &G;
    std::string &name_integrator;
    bool &without_controlled_step;
    Arguments_t(unsigned &number_bodies,type &time_end,unsigned &number_steps,unsigned &bits,
    type &precision_energy, std::string &name_integrator, unsigned &dim,
    type &distance_encounter, bool &without_controlled_step, type &G,
    unsigned &number_orbital_elements) : bits(bits), number_bodies(number_bodies),time_end(time_end),
    number_steps(number_steps),precision_energy(precision_energy),
    name_integrator(name_integrator), dim(dim), distance_encounter(distance_encounter),
    without_controlled_step(without_controlled_step), G(G),
    number_orbital_elements(number_orbital_elements){ }
    void Parser(int argc, char *const argv[]){
        int i, longIndex=0;
        std::string ss;
        static const char* shortopts = "b:M:t:N:p:d:e:G:";
        static const struct option longOpts[] ={
            {"integrator",required_argument, NULL,0},
            {"without_controlled_step",no_argument, NULL, 0} ,
            {"without_orbital_elements",required_argument, NULL, 0}
        };
        auto opt = getopt_long(argc,argv,shortopts,longOpts,&longIndex);
        while(opt != -1){
            switch(opt){
                case 'b':
                    #ifdef ANY_PRESICION
                        ss=optarg;
                        i = boost::lexical_cast<int>(ss);
                        if(i < 64){
                            std::cout << "Don't possible use less 64 bits";
                            exit(-1);
                        }else
                            bits = i;
                    #endif
                    
                    #ifdef DOUBLE_PRESICION
                        std::cout << "Precision is double, you can't change number of bits";
                        exit(-1);
                    #endif
                    
                    #ifdef FLOAT_PRESICION
                        std::cout << "Precision is float, you can't change number of bits";
                        exit(-1);
                    #endif
                    
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
                case 'G':
                    #ifdef ANY_PRESICION
                        G=mpf_class_set_str(optarg,bits);
                    #else
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            G = boost::lexical_cast<type>(ss);
                        }
                        catch(...){
                            std::cout << "String can't be converted to this type(" << 
                            typeid(type).name() << ")" << std::endl;
                            exit(-1);
                        }
                    #endif
                    break;
                case 't':
                    #ifdef ANY_PRESICION
                        time_end = mpf_class_set_str(optarg,bits);
                    #else
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            time_end = boost::lexical_cast<type>(ss);
                        }
                        catch(...){
                            std::cout << "String can't be converted to this type(" << 
                            typeid(type).name() << ")" << std::endl;
                            exit(-1);
                        }
                    #endif
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
                    #ifdef ANY_PRESICION
                        precision_energy = mpf_class_set_str(optarg,bits);
                    #else
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            precision_energy = boost::lexical_cast<type>(ss);
                        }
                        catch(...){
                            std::cout << "String can't be converted to this type(" << 
                            typeid(type).name() << ")" << std::endl;
                            exit(-1);
                        }
                    #endif
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
                    #ifdef ANY_PRESICION
                        distance_encounter = mpf_class_set_str(optarg,bits);
                    #else
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            distance_encounter = boost::lexical_cast<type>(ss);
                        }
                        catch(...){
                            std::cout << "String can't be converted to this type(" << 
                            typeid(type).name() << ")" << std::endl;
                            exit(-1);
                        }
                    #endif
                    break;
                case 0:
                    if("integrator" == longOpts[longIndex].name){
                        name_integrator = optarg;
                    }
                    if("without_controlled_step" == longOpts[longIndex].name){
                        without_controlled_step = 1;
                    }
                    if("without_orbital_elements" == longOpts[longIndex].name){
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            number_orbital_elements = boost::lexical_cast<unsigned>(ss);
                        }catch(...){
                            number_orbital_elements = 0;
                        }
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
    type min_distance=1.0;
    Encounter(const type& in) : distance_encounter(in){}
    bool operator()(const Type& X, const type& t){
        if(norm_2(subrange(row(X,0) - row(X,1),0,3)) <= distance_encounter){
            std::cout << "ENCOUNTER!! in time " << t << std::endl;
            return 0;
        }else{
            //std::cout << "min_distance = " << min_distance << std::endl;
            //std::cout << "current_distance= " << norm_2(subrange(row(X,0) - row(X,1),0,3)) << "   " << t << std::endl;
            if(min_distance > norm_2(subrange(row(X,0) - row(X,1),0,3)))
                min_distance = norm_2(subrange(row(X,0) - row(X,1),0,3));
            //std::cout << norm_2(subrange(row(X,0) - row(X,1),0,3)) << "   " << t << std::endl;
            return 1;
        }
    } 
};

template<typename Type, typename type>
struct SymplecticEncounter{
    type distance_encounter;
    type min_distance=1.0;
    SymplecticEncounter(const type& in) : distance_encounter(in){}
    bool operator()(const Type& X, const type& t){
        if(norm_2(row(X.first,0) - row(X.first,1)) <= distance_encounter){
            std::cout << "ENCOUNTER!! in time " << t << std::endl;
            return 0;
        }else{
            //std::cout << "min_distance = " << min_distance << std::endl;
            //std::cout << "current_distance= " << norm_2(row(X.first,0) - row(X.first,1)) << "   " << t << std::endl;
            if(min_distance > norm_2(row(X.first,0) - row(X.first,1)))
                min_distance = norm_2(row(X.first,0) - row(X.first,1));
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
    type h_begin=h;
    type T=t+h*n;
    Obs(X0,t);
    while(t<T && Bad(X0,t)){
        state=-1;
        std::cout << t << '\r';
        rk.do_step(sysF,X0,Y,t,h);
        if(Err(X0,Y)==-1){
            goto obs_point;
        } 
        while(Err(X0,Y)){
            h=h/2.;
            rk.do_step(sysF, X0,Y,t,h);
            //std::cout << "h_decrease= " << h << std::endl;
            //std::cout << "Y :" << std::endl << Y << std::endl;
            state=1;
        }
        if(state==1){
            //std::cout << "STOP1" << std::endl;
            h=h*2;
            //std::cout << "h_increase= " << h << std::endl;
            rk.do_step(sysF,X0,Y,t,h);
            goto obs_point;
        }
        while(!Err(X0,Y)){
            if(h > h_begin){
                state=0;
                break;
            }
            //std::cout << "STOP" << std::endl;
            h=h*2;
            //std::cout << "h_increase= " << h << std::endl;
            rk.do_step(sysF, X0,Y,t,h);
            state=0;
        }
        if(state==0){
            h=h/2;
            //std::cout << "h_decrease= " << h << std::endl; 
            rk.do_step(sysF,X0,Y,t,h);
        }
        obs_point:
        X0=Y;
        t+=h;
        Obs(X0,t);
    }
    std::cout << std::endl;
}

template<typename type>
struct Integrate_based_args_t{
    Arguments_t<type> args;
    std::string integrator_type;
    Integrate_based_args_t(const Arguments_t<type> &args) : args(args) { }
    std::string choice_integrator(){
        if (args.name_integrator == "RK4"){
            integrator_type = "standard";
            return integrator_type;
        }
        if (args.name_integrator == "RK5"){
            integrator_type = "standard";
            return integrator_type;
        }
        if (args.name_integrator == "LP"){
            integrator_type = "symplectic";
            return integrator_type;
        }
        throw;
    } 
    template<typename Integrator, typename sysf, typename observer,typename Type, 
    typename controller=std::function<int(const Type& X, const Type& Y)>, 
    typename exit=std::function<bool(const Type&X, const type &t)> >
    void integrator_call(Integrator rk, sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){
        if(args.distance_encounter == 0){
            if(args.without_controlled_step)
                integrate(rk, sysF, X0, t, h, n, Obs,[](const Type& X, const type& t){return 1;},
                [](const Type& X, const Type& Y){return -1;});
            else
                integrate(rk, sysF, X0, t, h, n, Obs,[](const Type& X, const type& t){return 1;},Err);
        }else{
            if(args.without_controlled_step)
                integrate(rk, sysF, X0, t, h, n, Obs,Bad,[](const Type& X, const Type& Y){return -1;});
            else
                integrate(rk, sysF, X0, t, h, n, Obs,Bad, Err);
        }
    }
    template<typename sysf, typename observer,typename Type, 
    typename controller=std::function<int(const Type& X, const Type& Y)>, 
    typename exit=std::function<bool(const Type&X, const type &t)>>
    void integrator_call_standard(sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){

        if (args.name_integrator == "RK4"){
            integrator_call(RungeKutta5_Fehlberg<type, Type>(),sysF,X0,t,h,n,Obs,Bad,Err);
        }
        if (args.name_integrator == "RK5"){
            integrator_call(RungeKutta4<type, Type>(),sysF,X0,t,h,n,Obs,Bad,Err);
        }
    }
    template<typename sysf, typename observer,typename Type, 
    typename controller=std::function<int(const Type& X, const Type& Y)>, 
    typename exit=std::function<bool(const Type&X, const type &t)>>
    void integrator_call_symplectic(sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){
        if (args.name_integrator == "LP"){
            integrator_call(LeapFrog<type, Type>(),sysF,X0,t,h,n,Obs,Bad,Err);
        }
    }
};


int main(int argc, char *const argv[]){
    std::ios_base::sync_with_stdio(false);
    
    std::ofstream out,energy,moon;
    std::ifstream file_in,file_size,file_mass,file_initial;
    std::vector<data_type> m;
    std::vector<state_type> state_result;
    std::vector<data_type> time_result, energy_result;
    std::vector<state_type >norm_R;
    std::vector<std::vector<state_vector> > orbital_elements_result;
    unsigned current=0;
    unsigned N,M, bits = 64, dim = 6;
    data_type pi=atan(data_type(1.))*4;
    data_type G=4*pi*pi;
    //G=1.0;
    data_type T,t0(0.), max_delta=0, delta, distance_encounter = 0, precision_energy = 0;
    state_type X0,Y;
    symplectic_type Z;
    std::string name_integrator = "RK4";
    std::string type_integrator = "standard";
    bool without_controlled_step = 0;
    unsigned number_orbital_elements;
    std::pair<unsigned,unsigned> two_bodies;
    std::vector<std::pair<unsigned,unsigned> > pairs_of_bodies;
    //runge_kutta4<state_type> rk;
    //controlled_stepper_type controlled_stepper;
    
    
    file_size.open("file_size.dat",std::ios_base::in);
    file_size >> T >> N >> M >> dim >> bits >> distance_encounter;
    while(file_size >> two_bodies.first >> two_bodies.second){
        pairs_of_bodies.push_back(two_bodies);
    }
    file_size.close();

    number_orbital_elements = pairs_of_bodies.size();
    
    std::cout << "distance_encounter= " << distance_encounter << std::endl;
    Arguments_t<data_type> keys(M,T,N,bits,precision_energy,name_integrator,dim,
    distance_encounter,without_controlled_step,G,number_orbital_elements);
    keys.Parser(argc,argv);

    #ifdef ANY_PRESICION
    std::cout.precision(bits);
    out.precision(bits);
    energy.precision(bits);
    mpf_set_default_prec(bits);
    std::cout << "PRESICION: " << mpf_get_default_prec() << std::endl;
    #endif

    #ifdef DOUBLE_PRESICION
    std::cout << "DOUBLE PRESICION" << std::endl;
    #endif

    #ifdef FLOAT_PRESICION
    std::cout << "FLOAT PRESICION" << std::endl;
    #endif
    
    std::cout << sqrt(1.0009534/3)*atan(data_type(1.0))*8 << std::endl;

    m.resize(M);
    X0.resize(M,dim);
    orbital_elements_result.resize(number_orbital_elements);

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
        type_integrator = Integrate_based_args.choice_integrator();
    }
    catch(...){
        std::cout << "This integrator isn't exists" << std::endl;
        exit(-1);
    }

    NBody nb(state_result, time_result, energy_result, precision_energy,
    m,current,dim,norm_R,G,orbital_elements_result, pairs_of_bodies,number_orbital_elements,
    "res.dat","orbital_elements.dat","energy.dat");
    SymplecticNBody snb(state_result, time_result, energy_result,  precision_energy, m,current,
    dim,norm_R,G, orbital_elements_result,pairs_of_bodies,number_orbital_elements,
    "res.dat","orbital_elements.dat","energy.dat");
    std::cout << type_integrator << std::endl;
    if(type_integrator == "standard"){
        Integrate_based_args.integrator_call_standard(nb , Y , t0, data_type(T/N), N, nb,
                Encounter<state_type, data_type>(distance_encounter), nb);
        std::cout << "steps=" << current << std::endl;
        std::cout << Y << std::endl;
    }
    if(type_integrator == "symplectic"){
        Z.first=subrange(Y,0,M,0,dim/2);
        Z.second=subrange(Y,0,M,dim/2,dim);
        Integrate_based_args.integrator_call_symplectic(snb , Z , t0, data_type(T/N), N, snb,
                SymplecticEncounter<symplectic_type, data_type>(distance_encounter), snb);
        std::cout << "steps=" << current << std::endl;
        std::cout << Z.first << std::endl;
        std::cout << Z.second << std::endl;
    }
    
    //size_t steps=integrate_n_steps( rk , nb , Y , 0., T/N, N, nb);
    //size_t steps=integrate_adaptive( make_controlled< error_stepper_type >( 1.0e-12 , 1.0e-8 ),
    //nb, Y, 0., T, T/N, nb);
    
    
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
    if(number_orbital_elements>0){
        moon.open("orbital_elements.dat",std::ios_base::trunc);
        for(unsigned i=0;i<=current;++i){
            for(unsigned l=0;l<number_orbital_elements;++l){
                for(unsigned k=0;k<4;++k)
                    moon << orbital_elements_result[l][i](k) << " ";
                moon << std::endl;
            }
            moon << std::endl << std::endl;
        }
        moon.close();
    }
    /*
      for(unsigned j=0;j<M;++j){
        out.open(filenamestr("NBodyRK_",j,".dat"),std::ios_base::trunc);
        std::cout << filenamestr("NBodyRK_",j,".dat") << std::endl;
        for(unsigned i=0;i<=current;++i){
            out << time_result[i] << " ";
            out << row(state_result[i],j) << std::endl;
        }        
        out.close();
    }
    */
    energy.open("energy.dat",std::ios_base::trunc);
    for(unsigned i=0;i<=current;++i){
        energy <<  "h(" << i << ")= " << energy_result[i] << std::endl;
        delta=ABS_((energy_result[0]-energy_result[i])/energy_result[i]);
        if(delta > max_delta) max_delta=delta;
    }
    std::cout << max_delta << std::endl;
    
   return 0;
}