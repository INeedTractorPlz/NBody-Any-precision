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

//#define INSTANT_OUTPUT_TO_FILE
#define MEMORY_FREEING_BY_CURRENT


#define DOUBLE_PRESICION
//#define ANY_PRESICION

#ifdef ANY_PRESICION
#define ABS_ abs
typedef mpf_class data_type;
#endif

#ifdef LONG_DOUBLE_PRESICION
#define ABS_ fabs
typedef long double data_type;
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



const char* filenamestr(const char* s1, unsigned j, const char* s2){
    static std::string ss;
    ss="";
    ss+= s1; ss+=boost::lexical_cast<std::string>(j); ss+=s2;
    return ss.c_str();
}

template<typename type>
std::ostream& operator<<(std::ostream& os, const std::pair<type,type> in){
    os << in.first << std::endl;
    os << in.second << std::endl;
    return os;
} 

unsigned mpf_to_unsigned(mpf_class in,mpfr_rnd_t rnd=MPFR_RNDN){
    mpfr_t out_t;
    unsigned result;
    mpf_class out;
    mpfr_init_set_f(out_t,in.get_mpf_t(),rnd);
    return mpfr_get_ui(out_t,rnd);
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
    //static std::ofstream res_file,energy_file,orbital_file;
    std::string res_file_name,energy_file_name,orbital_file_name,distance_file_name;
    data_type &max_delta,&energy0,&energy;
    #ifdef MEMORY_FREEING_BY_CURRENT
    std::vector< state_type >& m_states;
    std::vector< data_type >& m_times, &energy_result;
    std::vector<std::vector<state_vector> > &orbital_elements_result; 
    std::vector<state_type >& norm_R;
    #endif
    const std::vector<data_type>& m;
    std::vector<std::pair<unsigned, unsigned> > pairs_of_bodies;
    unsigned number_orbital_elements;
    data_type &G,pi = atan(data_type(1.0))*4;
    const data_type precision_energy;
    unsigned &current,number_steps;
    #ifdef MEMORY_FREEING_BY_CURRENT
    int &freeing_memory_current;
    #endif
    int &recording_frequency, &focal_body;
    unsigned number_bodies,dim;
    NBody(
    #ifdef MEMORY_FREEING_BY_CURRENT
    std::vector< state_type >& m_states, std::vector< data_type >& m_times, 
    std::vector<data_type>& energy_result, std::vector<state_type >& norm_R,
    std::vector<std::vector<state_vector> > &orbital_elements_result,
    #endif
    const std::vector<data_type>& m,unsigned &current,const unsigned &dim,
    const data_type &precision_energy, data_type &G,
    std::vector<std::pair<unsigned, unsigned> > pairs_of_bodies,
    unsigned number_orbital_elements
    ,const char *res_file_name,const char *energy_file_name, const char *orbital_file_name,
    const char *distance_file_name, data_type &max_delta,data_type &energy0,data_type &energy,
    #ifdef MEMORY_FREEING_BY_CURRENT
    int &freeing_memory_current,
    #endif
    int &recording_frequency, int &focal_body
    ):
    #ifdef MEMORY_FREEING_BY_CURRENT
    m_states(m_states),m_times(m_times), energy_result(energy_result),norm_R(norm_R),
    orbital_elements_result(orbital_elements_result),
    #endif
    m(m),current(current),dim(dim), number_bodies(m.size()), precision_energy(precision_energy), G(G),
    pairs_of_bodies(pairs_of_bodies),number_orbital_elements(number_orbital_elements)
    ,res_file_name(res_file_name),energy_file_name(energy_file_name),orbital_file_name(orbital_file_name),
    distance_file_name(distance_file_name),max_delta(max_delta),energy0(energy0), energy(energy),
    #ifdef MEMORY_FREEING_BY_CURRENT
    freeing_memory_current(freeing_memory_current),
    #endif
    recording_frequency(recording_frequency),focal_body(focal_body)
    {  
        max_delta = 0;
    }

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
        if(t!=0) ++number_steps;
        
        if(number_steps%recording_frequency == 0){
            if(t!=0) ++current;
            #if defined(MEMORY_FREEING_BY_CURRENT)
            if(freeing_memory_current != -1){
                if(current == freeing_memory_current){
                    output_to_file();
                    std::cout << "Clear memory current = " << current << std::endl;
                    energy_result.clear();
                    m_states.clear();
                    m_times.clear();
                    for(unsigned i = 0; i < number_orbital_elements; ++i)
                        orbital_elements_result[i].clear();
                    norm_R.clear();
                    current = 0;
                }
            }
            result_record(R,t);
            energy = energy_result[current];
            #endif
        
            #ifdef MEMORY_FREEING_BY_CURRENT
            if(freeing_memory_current != -1){
            #endif
                if(t == 0){
                    state_type norm_r(number_bodies,number_bodies);
                    norm(R,norm_r);
                    energy0 = EnergyIntegral(R,norm_r);    
                    energy = energy0;
                }
            #ifdef MEMORY_FREEING_BY_CURRENT
            }
            #endif

            #ifdef INSTANT_OUTPUT_TO_FILE
            output_to_file(R,t);
            #endif
        }else{
            state_type norm_r(number_bodies,number_bodies);
            norm(R,norm_r);
            energy = EnergyIntegral(R,norm_r);
        }
    }

    bool operator()(const state_type& X, const state_type& Y){
        state_type norm_Y=std::move(state_type(number_bodies,number_bodies));
        norm(Y,norm_Y);
        data_type energy_Y(EnergyIntegral(Y,norm_Y));
        data_type delta=ABS_((energy-energy_Y)/energy);
        //std::cout << "delta= " << delta << std::endl;
        if(delta > precision_energy) return 1;
        else return 0;
    }
    void norm(const state_type& R, state_type& norm_r);
    data_type EnergyIntegral(const state_type& R, const state_type& norm_r);
    void get_orbital_elements(const state_vector &R, unsigned i, unsigned j,
    state_vector& orbital_elements,const state_type &norm_r);
    #ifdef MEMORY_FREEING_BY_CURRENT
    void result_record(const state_type &R, data_type t);
    #endif
    #ifdef INSTANT_OUTPUT_TO_FILE
    void output_to_file(const state_type &R, data_type t);
    #endif
    #ifdef MEMORY_FREEING_BY_CURRENT
    void output_to_file();
    #endif
};

struct SymplecticNBody : NBody{
    SymplecticNBody(
    #ifdef MEMORY_FREEING_BY_CURRENT
    std::vector< state_type >& m_states, std::vector< data_type >& m_times, 
    std::vector<data_type>& energy_result,std::vector<state_type >& norm_R,
    std::vector<std::vector<state_vector> > &orbital_elements_result, 
    #endif
    const std::vector<data_type>& m,unsigned &current,const unsigned &dim,
    const data_type &precision_energy, data_type &G,
    std::vector<std::pair<unsigned,unsigned> > pairs_of_bodies,
    unsigned number_orbital_elements
    ,const char *res_file_name,const char *energy_file_name, const char *orbital_file_name,
    const char *distance_file_name, data_type &max_delta,data_type &energy0,data_type &energy,
    #ifdef MEMORY_FREEING_BY_CURRENT
    int &freeing_memory_current,
    #endif
    int &recording_frequency, int &focal_body
    ):
    NBody(
    #ifdef MEMORY_FREEING_BY_CURRENT
    m_states, m_times, energy_result,norm_R,orbital_elements_result, 
    #endif
    m,current, dim, precision_energy, G,pairs_of_bodies,
    number_orbital_elements,res_file_name,energy_file_name,orbital_file_name,
    distance_file_name,max_delta,energy0,energy,
    #ifdef MEMORY_FREEING_BY_CURRENT
    freeing_memory_current,
    #endif
    recording_frequency,focal_body
    ){  }
    
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
        if(t!=0) ++number_steps;
        
        state_type R_matrix(number_bodies,dim);
        subrange(R_matrix,0,number_bodies,0,dim/2)=R.first;
        subrange(R_matrix,0,number_bodies,dim/2,dim)=R.second;
            
        if(number_steps%recording_frequency == 0){
            if(t!=0){
                ++current;
            }
                    
            #if defined(MEMORY_FREEING_BY_CURRENT)
            if(freeing_memory_current != -1){
                if(current == freeing_memory_current){
                    output_to_file();
                    energy_result.clear();
                    m_states.clear();
                    m_times.clear();
                    for(unsigned i = 0; i < number_orbital_elements; ++i)
                        orbital_elements_result[i].clear();
                    norm_R.clear();
                    current = 0;
                }
            }
            result_record(R_matrix,t);
            energy = energy_result[current];
            #endif

            #ifdef MEMORY_FREEING_BY_CURRENT
            if(freeing_memory_current != -1){
            #endif
                if(t==0){
                    state_type norm_r(number_bodies,number_bodies);
                    norm(R_matrix,norm_r);
                    energy0 = EnergyIntegral(R_matrix,norm_r);    
                    energy=energy0;
                }
            #ifdef MEMORY_FREEING_BY_CURRENT
            }
            #endif

            #ifdef INSTANT_OUTPUT_TO_FILE
            output_to_file(R_matrix,t);
            #endif
        }else{
            state_type norm_r(number_bodies,number_bodies);
            norm(R_matrix,norm_r);
            energy = EnergyIntegral(R_matrix,norm_r);
        }
    }
    bool operator()(const symplectic_type& X, const symplectic_type& Y){
        state_type norm_Y=std::move(state_type(number_bodies,number_bodies));
        state_type Y_matrix(number_bodies,dim);
        subrange(Y_matrix,0,number_bodies,0,dim/2)=Y.first;
        subrange(Y_matrix,0,number_bodies,dim/2,dim)=Y.second;
        
        norm(Y_matrix,norm_Y);
        data_type energy_Y(EnergyIntegral(Y_matrix,norm_Y));
        data_type delta=ABS_((energy-energy_Y)/energy);
        if(delta > precision_energy) return 1;
        else return 0;
    }
};



#ifdef MEMORY_FREEING_BY_CURRENT    
void NBody::result_record(const state_type &R, data_type t){
    m_states.push_back(R); 
    m_times.push_back(t);
    norm_R.push_back(std::move(state_type(number_bodies,number_bodies))); 
    norm(R,norm_R[current]);
    energy_result.push_back(EnergyIntegral(R,norm_R[current]));
    for(unsigned i=0;i<number_orbital_elements;++i){
        state_vector distance_i = row(R,pairs_of_bodies[i].first) - row(R,pairs_of_bodies[i].second);
        //std::for_each(distance_i.begin(),distance_i.end(),[](data_type &x){x=ABS_(x);});
        orbital_elements_result[i].push_back(std::move(state_vector(4)));
        get_orbital_elements(distance_i,pairs_of_bodies[i].first,pairs_of_bodies[i].second,
        orbital_elements_result[i][current],norm_R[current]);
    }
}
#endif

#ifdef INSTANT_OUTPUT_TO_FILE
void NBody::output_to_file(const state_type &R, data_type t){
    std::ofstream res_file,energy_file,orbital_file,distance_between_pairs;
    
    if(focal_body != -1){
        res_file.open(res_file_name,std::ios_base::app);
        for(unsigned j=0;j<number_bodies;++j){
            for(unsigned k=0;k<dim;++k)
                res_file << R(j,k) - R(focal_body,k) << " ";
            res_file << std::endl;
        }
        res_file << std::endl << std::endl;
        res_file.close();
    }else{
        res_file.open(res_file_name,std::ios_base::app);
        for(unsigned j=0;j<number_bodies;++j){
            for(unsigned k=0;k<dim;++k)
                res_file << R(j,k) << " ";
            res_file << std::endl;
        }
        res_file << std::endl << std::endl;
        res_file.close();
    }
    state_type norm_r(number_bodies,number_bodies);
    norm(R,norm_r);
    
    distance_between_pairs.open(distance_file_name,std::ios_base::app);
    distance_between_pairs << t << " ";
    for(unsigned l = 0; l < number_orbital_elements; ++l)
            distance_between_pairs << 
            norm_r(pairs_of_bodies[l].first,pairs_of_bodies[l].second) << " ";
    distance_between_pairs << std::endl;
    distance_between_pairs << std::endl << std::endl;
    distance_between_pairs.close();


    if(number_orbital_elements>0){
        orbital_file.open(orbital_file_name,std::ios_base::app);
        for(unsigned i=0;i<number_orbital_elements;++i){
            state_vector distance_i = row(R,pairs_of_bodies[i].first)-row(R,pairs_of_bodies[i].second);
            std::for_each(distance_i.begin(),distance_i.end(),[](data_type &x){x=ABS_(x);});
            state_vector oe_out(4); 
            get_orbital_elements(distance_i,pairs_of_bodies[i].first,pairs_of_bodies[i].second,
            oe_out,norm_r);
            for(unsigned k=0;k<4;++k)
                orbital_file << oe_out(k) << " ";
            orbital_file << std::endl;
        }
        orbital_file << std::endl << std::endl;
        orbital_file.close();
    }
    
    energy = EnergyIntegral(R,norm_r);
    energy_file.open(energy_file_name,std::ios_base::app);
        energy_file <<  "h(" << current << ")= " << energy << std::endl;
        data_type delta=ABS_((energy0-energy)/energy0);
            if(delta > max_delta) max_delta=delta;
    energy_file.close();
}
#endif
#ifdef MEMORY_FREEING_BY_CURRENT
void NBody::output_to_file(){
    std::ofstream res_file,energy_file,orbital_file, distance_between_pairs;
    
    if(focal_body != -1){
        res_file.open(res_file_name,std::ios_base::app);
        for(unsigned i=0;i<current;++i){
            for(unsigned j=0;j<number_bodies;++j){
                for(unsigned k=0;k<dim;++k)
                    res_file << m_states[i](j,k) - m_states[i](focal_body,k) << " ";
                res_file << std::endl;
            }
            res_file << std::endl << std::endl;
        }
        res_file.close();
    }else{
        res_file.open(res_file_name,std::ios_base::app);
        for(unsigned i=0;i<current;++i){
            for(unsigned j=0;j<number_bodies;++j){
                for(unsigned k=0;k<dim;++k)
                    res_file << m_states[i](j,k) << " ";
                res_file << std::endl;
            }
            res_file << std::endl << std::endl;
        }
        res_file.close();
    }

    distance_between_pairs.open(distance_file_name,std::ios_base::app);
    for(unsigned i=0;i<current;++i){
        distance_between_pairs << m_times[i] << " ";
        for(unsigned l = 0; l < number_orbital_elements; ++l)
            distance_between_pairs << 
            norm_R[i](pairs_of_bodies[l].first,pairs_of_bodies[l].second) << " ";
        distance_between_pairs << std::endl;
        distance_between_pairs << std::endl << std::endl;
    }
    distance_between_pairs.close();

    if(number_orbital_elements>0){
        orbital_file.open(orbital_file_name,std::ios_base::app);
        for(unsigned i=0;i<current;++i){
            for(unsigned l=0;l<number_orbital_elements;++l){
                orbital_file << m_times[i] << " ";
                for(unsigned k=0;k<4;++k)
                    orbital_file << orbital_elements_result[l][i](k) << " ";
                orbital_file << std::endl << std::endl;
            }
            orbital_file << std::endl << std::endl;
        }
        orbital_file.close();
    }
    energy_file.open(energy_file_name,std::ios_base::app);
    for(unsigned i=0;i<current;++i){
        energy_file <<  "h(" << i << ")= " << energy_result[i] << std::endl;
        data_type delta=ABS_((energy0-energy_result[i])/energy0);
        if(delta > max_delta) max_delta=delta;
    }

}
#endif
void NBody::get_orbital_elements(const state_vector &R, unsigned i, unsigned j,
state_vector& orbital_elements,const state_type &norm_r){
    data_type kappa_quad,r,v;
    data_type a,e,inclination;
    
    kappa_quad=G*(m[i]+m[j]);
    r=norm_r(i,j);
    v=norm_2(subrange(R,dim/2,dim));
    
    data_type h = v*v/2 - kappa_quad/r;
    state_vector c = cross_product(subrange(R,0,dim/2),subrange(R,dim/2,dim));
    
    e = sqrt(1.+2.*h*inner_prod(c,c)/(kappa_quad*kappa_quad));
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
    unsigned &bits, &number_bodies, &number_steps, &dim,&number_orbital_elements,
    &number_encounters;
    type &time_end, &precision_energy, &G, &step;
    std::string &name_integrator;
    bool &without_controlled_step;
    #ifdef MEMORY_FREEING_BY_CURRENT
    int &freeing_memory_current;
    #endif
    int &recording_frequency, &focal_body;
    Arguments_t(unsigned &number_bodies,type &time_end,unsigned &number_steps,unsigned &bits,
    type &precision_energy, std::string &name_integrator, unsigned &dim,
    unsigned &number_encounters, bool &without_controlled_step, type &G, type &step,
    unsigned &number_orbital_elements,
    #ifdef MEMORY_FREEING_BY_CURRENT
    int &freeing_memory_current,
    #endif
    int &recording_frequency, int &focal_body
    ) : bits(bits), number_bodies(number_bodies),time_end(time_end),
    number_steps(number_steps),precision_energy(precision_energy),
    name_integrator(name_integrator), dim(dim), number_encounters(number_encounters),
    without_controlled_step(without_controlled_step), G(G),step(step),
    number_orbital_elements(number_orbital_elements),
    #ifdef MEMORY_FREEING_BY_CURRENT
    freeing_memory_current(freeing_memory_current),
    #endif
    recording_frequency(recording_frequency),focal_body(focal_body)
    { }
    void Parser(int argc, char *const argv[]){
        int i, longIndex=0;
        std::string ss;
        static const char* shortopts = "b:M:t:N:p:d:G:h:";
        static const struct option longOpts[] ={
            {"integrator",required_argument, NULL,0},
            {"without_controlled_step",no_argument, NULL, 0} ,
            {"number_orbital_elements",required_argument, NULL, 0},
            #ifdef MEMORY_FREEING_BY_CURRENT
            {"freeing_memory_current",required_argument, NULL, 0},
            #endif
            {"recording_frequency",required_argument, NULL, 0},
            {"focal_body",required_argument, NULL, 0},
            {"number_encounters",required_argument, NULL, 0}
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
                    try {
                        i = boost::lexical_cast<int>(ss);
                    }catch(...){
                        #ifdef ANY_PRESICION
                        try{
                            i = mpf_to_unsigned(mpf_class_set_str(optarg,bits));
                        }
                        #else
                        try{
                            i = static_cast<int>(boost::lexical_cast<type>(ss));
                        }
                        #endif
                        catch(...){
                            std::cout << "This string can't be converted to int or type " << 
                            typeid(type).name() << std::endl;
                        }
                    }
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
                case 'h':
                    #ifdef ANY_PRESICION
                        step = mpf_class_set_str(optarg,bits);
                    #else
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            step = boost::lexical_cast<type>(ss);
                        }
                        catch(...){
                            std::cout << "This string can't be converted to this type(" << 
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
                    if("number_orbital_elements" == longOpts[longIndex].name){
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            number_orbital_elements = boost::lexical_cast<unsigned>(ss);
                        }catch(...){
                            number_orbital_elements = 0;
                        }
                    }
                    if("number_encounters" == longOpts[longIndex].name){
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            number_encounters = boost::lexical_cast<unsigned>(ss);
                        }catch(...){
                            number_encounters = 0;
                        }
                    }
                    if("focal_body" == longOpts[longIndex].name){
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            focal_body = boost::lexical_cast<int>(ss);
                        }catch(...){
                            std::cout << "Unsuitable type for number of body" << std::endl;
                        }
                    }
                    if("recording_frequency"== longOpts[longIndex].name){
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            if(ss == "default"){
                                recording_frequency = 100;
                            }else{                        
                                recording_frequency = boost::lexical_cast<int>(ss);
                            }
                        }catch(...){
                        }
                    }
                    
                    #ifdef MEMORY_FREEING_BY_CURRENT
                    if("freeing_memory_current"== longOpts[longIndex].name){
                        try{
                            std::string ss=boost::lexical_cast<std::string>(optarg);
                            if(ss == "default"){
                                freeing_memory_current = 1000000;
                            }else{                        
                                freeing_memory_current = boost::lexical_cast<int>(ss);
                            }
                        }catch(...){
                        }
                    }
                    #endif
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
    std::vector<std::pair<unsigned,unsigned> > pairs_of_encounters;
    std::vector<type> encounter_distances;
    unsigned number_encounters;
    Encounter(const std::vector<std::pair<unsigned,unsigned> >& pairs_of_encounters,
    const std::vector<type> &encounter_distances, unsigned number_encounters)
    : pairs_of_encounters(pairs_of_encounters),encounter_distances(encounter_distances),
    number_encounters(number_encounters)
    {   }
    bool operator()(const Type& X, const type& t){
        for(unsigned l = 0; l < number_encounters; ++l){
            if(norm_2(subrange(row(X,pairs_of_encounters[l].first)
            - row(X,pairs_of_encounters[l].second),0,3)) <= encounter_distances[l]){
                std::cout << "ENCOUNTER!! between " << pairs_of_encounters[l].first 
                << " and " << pairs_of_encounters[l].second 
                << " in time " << t << std::endl;
                return 0;
            }
        }
        return 1;
    } 
};

template<typename Type, typename type>
struct SymplecticEncounter{
    std::vector<std::pair<unsigned,unsigned> > pairs_of_encounters;
    std::vector<type> encounter_distances;
    unsigned number_encounters;
    SymplecticEncounter(const std::vector<std::pair<unsigned,unsigned> >& pairs_of_encounters,
    const std::vector<type> &encounter_distances, unsigned number_encounters)
    : pairs_of_encounters(pairs_of_encounters),encounter_distances(encounter_distances),
    number_encounters(number_encounters)
    {   }
    bool operator()(const Type& X, const type& t){
        for(unsigned l = 0; l < number_encounters; ++l){
            if(norm_2(row(X.first,pairs_of_encounters[l].first)
            - row(X.first,pairs_of_encounters[l].second)) <= encounter_distances[l]){
                std::cout << "ENCOUNTER!! between " << pairs_of_encounters[l].first 
                << " and " << pairs_of_encounters[l].second 
                << " in time " << t << std::endl;
                return 0;
            }
        }
        return 1;
    } 
};

template<typename type, typename Integrator, typename sysf, typename observer,typename Type, 
typename controller=std::function<int(const Type& X, const Type& Y)>, 
typename exit=std::function<bool(const Type&X, const type &t)> >
int integrate(Integrator rk, sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){
    int state, number_steps = 0;
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
        ++number_steps;
        Obs(X0,t);
    }
    std::cout << std::endl;
    return number_steps;
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
    int integrator_call(Integrator rk, sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){
        if(args.number_encounters == 0){
            if(args.without_controlled_step)
                return integrate(rk, sysF, X0, t, h, n, Obs,[](const Type& X, const type& t){return 1;},
                [](const Type& X, const Type& Y){return -1;});
            else
                return integrate(rk, sysF, X0, t, h, n, Obs,[](const Type& X, const type& t){return 1;},Err);
        }else{
            if(args.without_controlled_step)
                return integrate(rk, sysF, X0, t, h, n, Obs,Bad,[](const Type& X, const Type& Y){return -1;});
            else
                return integrate(rk, sysF, X0, t, h, n, Obs,Bad, Err);
        }
    }
    template<typename sysf, typename observer,typename Type, 
    typename controller=std::function<int(const Type& X, const Type& Y)>, 
    typename exit=std::function<bool(const Type&X, const type &t)>>
    int integrator_call_standard(sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){

        if (args.name_integrator == "RK4"){
            return integrator_call(RungeKutta5_Fehlberg<type, Type>(),sysF,X0,t,h,n,Obs,Bad,Err);
        }
        if (args.name_integrator == "RK5"){
            return integrator_call(RungeKutta4<type, Type>(),sysF,X0,t,h,n,Obs,Bad,Err);
        }
    }
    template<typename sysf, typename observer,typename Type, 
    typename controller=std::function<int(const Type& X, const Type& Y)>, 
    typename exit=std::function<bool(const Type&X, const type &t)>>
    int integrator_call_symplectic(sysf sysF, Type& X0, type t, type h, int n, observer Obs,
                exit Bad = [](const Type& X, const type& t){return 1;},
                controller Err = [](const Type& X, const Type& Y){return -1;}){
        if (args.name_integrator == "LP"){
            return integrator_call(LeapFrog<type, Type>(),sysF,X0,t,h,n,Obs,Bad,Err);
        }
    }
};


int main(int argc, char *const argv[]){
    std::ios_base::sync_with_stdio(false);
    
    #ifdef MEMORY_FREEING_BY_CURRENT
    std::vector<state_type> state_result;
    std::vector<data_type> time_result, energy_result;
    std::vector<state_type >norm_R;
    std::vector<std::vector<state_vector> > orbital_elements_result;
    #endif
    std::ofstream res_file,energy_file,orbital_file,distance_between_pairs;
    std::ifstream file_in,file_size,file_mass,file_initial,file_encounters;
    std::vector<data_type> m;
    unsigned current = 0;
    #ifdef MEMORY_FREEING_BY_CURRENT
    int freeing_memory_current = -1;
    #endif
    int number_steps,recording_frequency = 1, focal_body = -1;
    unsigned N,M, bits = 64, dim = 6;
    data_type pi=atan(data_type(1.))*4;
    data_type G=4*pi*pi;
    //G=1.0;
    data_type T,t0(0.), max_delta=0, delta, precision_energy = 0,
    step = 0,energy0,energy;
    state_type X0,Y;
    symplectic_type Z;
    std::string name_integrator = "RK4";
    std::string type_integrator = "standard";
    bool without_controlled_step = 0;
    unsigned number_orbital_elements, number_encounters = 0;
    std::pair<unsigned,unsigned> two_bodies;
    std::vector<std::pair<unsigned,unsigned> > pairs_of_bodies, pairs_of_encounters;
    std::vector<data_type> encounter_distances;
    data_type distance_encounter;
    //runge_kutta4<state_type> rk;
    //controlled_stepper_type controlled_stepper;
    
    
    file_size.open("file_size.dat",std::ios_base::in);
    file_size >> T >> N >> M >> dim >> bits;
    while(file_size >> two_bodies.first >> two_bodies.second){
        pairs_of_bodies.push_back(two_bodies);
    }
    file_size.close();

    number_orbital_elements = pairs_of_bodies.size();
    
    file_encounters.open("file_encounters.dat",std::ios_base::in);
    while(file_encounters >> two_bodies.first >> two_bodies.second >> distance_encounter){
        pairs_of_encounters.push_back(two_bodies);
        encounter_distances.push_back(distance_encounter);
    }
    file_encounters.close();

    number_encounters = encounter_distances.size();


    Arguments_t<data_type> keys(M,T,N,bits,precision_energy,name_integrator,dim,
    number_encounters,without_controlled_step,G,step,number_orbital_elements,
    #ifdef MEMORY_FREEING_BY_CURRENT
    freeing_memory_current,
    #endif
    recording_frequency,focal_body
    );
    keys.Parser(argc,argv);


    std::cout << "number_orbital_elements= " << number_orbital_elements << std::endl;
    std::cout << "number_encounters= " << number_encounters << std::endl;

    #ifdef ANY_PRESICION
    std::cout.precision(bits);
    res_file.precision(bits);
    energy_file.precision(bits);
    orbital_file.precision(bits);
    mpf_set_default_prec(bits);
    std::cout << "PRESICION: " << mpf_get_default_prec() << std::endl;
    #endif

    #ifdef LONG_DOUBLE_PRESICION
    std::cout << "LONG DOUBLE PRESICION" << std::endl;
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
    #ifdef MEMORY_FREEING_BY_CURRENT
    orbital_elements_result.resize(number_orbital_elements);
    #endif

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

    NBody nb(
    #ifdef MEMORY_FREEING_BY_CURRENT
    state_result, time_result, energy_result,norm_R,orbital_elements_result, 
    #endif
    m,current,dim,precision_energy,G, pairs_of_bodies,number_orbital_elements
    ,"res.dat","energy.dat","orbital_elements.dat","distance_between_pairs.dat",max_delta,
    energy0,energy,
    #ifdef MEMORY_FREEING_BY_CURRENT
    freeing_memory_current,
    #endif
    recording_frequency,focal_body
    );
    SymplecticNBody snb(
    #ifdef MEMORY_FREEING_BY_CURRENT
    state_result, time_result, energy_result, norm_R, orbital_elements_result,
    #endif
    m,current,dim,precision_energy,G, pairs_of_bodies,number_orbital_elements
    ,"res.dat","energy.dat","orbital_elements.dat","distance_between_pairs.dat",max_delta,
    energy0,energy,
    #ifdef MEMORY_FREEING_BY_CURRENT
    freeing_memory_current,
    #endif
    recording_frequency,focal_body
    );

    //NBody &nb_reference=nb;
    //SymplecticNBody &snb_reference=snb;

    res_file.open("res.dat",std::ios_base::trunc);
    res_file.close();
    energy_file.open("energy.dat",std::ios_base::trunc);
    energy_file.close();
    orbital_file.open("orbital_elements.dat",std::ios_base::trunc);
    orbital_file.close();
    distance_between_pairs.open("distance_between_pairs.dat",std::ios_base::trunc);
    distance_between_pairs.close();

    if(step == 0){
        step = data_type(T/N);
    }else{
        #ifdef ANY_PRESICION
        N=mpf_to_unsigned(T/step + data_type(1.));
        #else
        N=ceil(T/step);
        #endif
    }

    std::cout << type_integrator << std::endl;
    if(type_integrator == "standard"){
        number_steps = Integrate_based_args.integrator_call_standard(nb, Y , t0, step, 
        N, nb, 
        Encounter<state_type, data_type>(pairs_of_encounters,
        encounter_distances,number_encounters),nb); 
        std::cout << "steps=" << number_steps << std::endl;
        std::cout << Y << std::endl;
    }
    if(type_integrator == "symplectic"){
        Z.first=subrange(Y,0,M,0,dim/2);
        Z.second=subrange(Y,0,M,dim/2,dim);
        number_steps = Integrate_based_args.integrator_call_symplectic(snb, Z , t0, step,
        N, snb,
        SymplecticEncounter<symplectic_type, data_type>(pairs_of_encounters,
        encounter_distances,number_encounters),snb);  
        std::cout << "steps=" << number_steps << std::endl;
        std::cout << Z.first << std::endl;
        std::cout << Z.second << std::endl;
    }
    
    //size_t steps=integrate_n_steps( rk , nb , Y , 0., T/N, N, nb);
    //size_t steps=integrate_adaptive( make_controlled< error_stepper_type >( 1.0e-12 , 1.0e-8 ),
    //nb, Y, 0., T, T/N, nb);
    
    #ifdef INSTANT_OUTPUT_TO_FILE
    std::cout << "max_delta= " << max_delta << std::endl;
    #endif

    #ifdef MEMORY_FREEING_BY_CURRENT
    res_file.open("res.dat",std::ios_base::app);
    if(focal_body !=-1){
        for(unsigned i=0;i<=current;++i){
            for(unsigned j=0;j<M;++j){
                for(unsigned k=0;k<dim;++k)
                    res_file << state_result[i](j,k) - state_result[i](focal_body,k) << " ";
                res_file << std::endl;
            }
            res_file << std::endl << std::endl;
        }
        res_file.close();
    }else{
        for(unsigned i=0;i<=current;++i){
            for(unsigned j=0;j<M;++j){
                for(unsigned k=0;k<dim;++k)
                    res_file << state_result[i](j,k) << " ";
                res_file << std::endl;
            }
            res_file << std::endl << std::endl;
        }
        res_file.close();
    }

    distance_between_pairs.open("distance_between_pairs.dat",std::ios_base::app);
    for(unsigned i=0;i<=current;++i){
        distance_between_pairs << time_result[i] << " ";
        for(unsigned l = 0; l < number_orbital_elements; ++l)
            distance_between_pairs << 
            norm_R[i](pairs_of_bodies[l].first,pairs_of_bodies[l].second) << " ";
        distance_between_pairs << std::endl;
        distance_between_pairs << std::endl << std::endl;
    }
    distance_between_pairs.close();

    
    if(number_orbital_elements>0){
        orbital_file.open("orbital_elements.dat",std::ios_base::app);
        for(unsigned i=0;i<=current;++i){
            for(unsigned l=0;l<number_orbital_elements;++l){
                orbital_file << time_result[i] << " ";
                for(unsigned k=0;k<4;++k)
                    orbital_file << orbital_elements_result[l][i](k) << " ";
                orbital_file << std::endl << std::endl << std::endl;
            }
            orbital_file << std::endl << std::endl;
        }
        orbital_file.close();
    }
    /*
      for(unsigned j=0;j<M;++j){
        res_file.open(filenamestr("NBodyRK_",j,".dat"),std::ios_base::trunc);
        std::cout << filenamestr("NBodyRK_",j,".dat") << std::endl;
        for(unsigned i=0;i<=current;++i){
            res_file << time_result[i] << " ";
            res_file << row(state_result[i],j) << std::endl;
        }        
        res_file.close();
    }
    */
    energy_file.open("energy.dat",std::ios_base::app);    
    for(unsigned i=0;i<=current;++i){
        energy_file <<  "h(" << i << ")= " << energy_result[i] << std::endl;
        delta=ABS_((energy_result[0]-energy_result[i])/energy_result[i]);
        if(delta > max_delta) max_delta=delta;
    }
    std::cout << "max_delta= " << max_delta << std::endl;
   #endif

   return 0;
}