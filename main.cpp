#include "biclique.h"
using namespace std;
using namespace std::chrono;

extern long double _cate, _wedge, _btf;
vector<long double> estis, naive_estis;
extern bool one_round, count_cate; 

long double p, Eps, Eps0, Eps1, Eps2;
// private degrees:
vector<int> priv_deg;
int priv_dmax_1, priv_dmax_2, iteration, num_rounds; 
bool input_epsilon =true; 
extern int alpha;
extern bool count_cc; 
extern stats::rand_engine_t engine;  // Declare the engine as extern
long double m3__ = 0, m2__ = 0, m1__ = 0, m0__ = 0, real_stars = 0,
			RR_time, server_side_time, naive_server_side, local_count_time = 0, deg_esti_time = 0,
			communication_cost = 0;


extern bool eva_comm ;


unsigned long long int  real ;

int P___, K___ ; 
// biclique related 
vector<vector<int>> up_options, lo_options; 

long double avg_estimated_variance = 0 ;

extern bool two_noisy_graph_switch; 

extern bool multi_estimator_switch ; 

// I  think we can have another way of counting caterpillars. 
// for each vertex, we consider the wedges. 
int main(int argc, char *argv[]) {

    eva_comm = false;

    if(eva_comm){
        cout<<"evaluating communication costs"<<endl;
    }

    // input parser:
    long double param = stold(argv[1]);
    string dataset = argv[2];
    num_rounds = atoi(argv[3]);
    int algorithm_switch = atoi(argv[4]);

    P___ = atoi(argv[5]);
    K___ = atoi(argv[6]);

    cout<<"P___ = "<<P___ <<endl;
    cout<<"K___ = "<<K___ <<endl;

    // initialize time
    RR_time = 0, server_side_time = 0, naive_server_side = 0;

    std::mt19937 rng(std::random_device{}());  // for seeding

    BiGraph g(dataset);

    
    /*
    long double S1 = binomial(g.num_v1, P___) * g.v1_max_degree;
    long double S2 = binomial(g.num_v2, K___) * g.v2_max_degree;
    cout<<"g.v1_max_degree = "<<g.v1_max_degree <<endl;
    cout<<"g.v2_max_degree = "<<g.v2_max_degree <<endl;

    cout<<"S1 = "<<S1 <<endl;
    cout<<"S2 = "<<S2 <<endl;

    if(g.num_v1 < g.num_v2 && S1 > S2){
        cout<<"found it"<<endl;
    }
    // */



    long double fill_rate = g.num_edges * 1.0 / ((double)g.num_v1 * (double)g.num_v2);

    cout << "fill_rate = " << fill_rate << endl;

    Eps = param;
    cout << "Eps = " << Eps << endl;

    vector<long double> m__;

    // unsigned long long int btf = BFC_EVP(g);
    // cout<<"BTF =  "<<btf<<endl;

    // grab exact biclique coutns from the sqlite database
    fetch_or_compute_biclique_count(P___, K___, dataset, g); 

    estis.resize(num_rounds);
    naive_estis.resize(num_rounds);
    vector<long double> rel_err;

    double t0 = omp_get_wtime();

    for (iteration = 0; iteration < num_rounds; iteration++) {

        // the naive algorithm for pq bcilqiue counting.
        if (algorithm_switch == 0) {
            cout << "Naive algorithm for biclique counting" << endl;
            cout << "EPS = " << Eps << endl;     

            p = 1.0 / (exp(Eps) + 1.0);
            unsigned int seed = rng();
            cout << "random seed = " << seed << endl;


            // printMemoryUsage();
            estis[iteration] = naive_biclique(g, seed, P___, K___);
            // printMemoryUsage();

            cout << "estimate = " << estis[iteration] << endl;

            long double relative_error = abs(estis[iteration] - real) * 1.0 / real;

            cout << "relative error = " << relative_error << endl;
            rel_err.push_back(relative_error);
            cout << endl;
        }
        if (algorithm_switch == 1) {
            // one round algorithms (same communication costs as naive)
            if(P___ == 2 && K___== 2){     
                cout << "Oneround-BTF (existing DBE algorithm)" << endl;
                cout << "EPS = " << Eps << endl;
                // flip probability
                p = 1.0 / (exp(Eps) + 1.0);
                unsigned int seed = rng();
                cout << "seed = " << seed << endl;

                bool avg_btf_switch = false;
                if(avg_btf_switch){
                    // why can we do this? 
                    // in theory, we are able to build two noisy graphs. 
                    // based on each, we can have a BTF estimate. 
                    // however, in experiment evaluations, we dont consider this case 
                    // because the double noisy graph technique
                    // is our contribution to the ADV algorithm
                    long double esti1 = one_round_btf(g, seed);
                    unsigned int seed2 = rng();
                    long double esti2 = one_round_btf(g, seed2);
                    estis[iteration] = (esti1 + esti2)/2;
                }else{
                    // this is a specific algorithm, so it is faster
                    estis[iteration] = one_round_btf(g, seed);

                    // let's see if ths is faster 
                    // estis[iteration] = one_round_biclique_2_K(g, K___, seed); 
                }
                long double rel = abs(estis[iteration] - real) * 1.0 / real;
                cout << "estimate = " << estis[iteration] << endl;
                cout << "relative error = " << rel << endl;
                rel_err.push_back(rel);
                cout << endl;
            }
            if(P___ == 2 && K___ > 2){
                // one-round biclique algorithm for general P, Q values
                cout << "Oneround" << endl;
                cout<<"P___ == " << P___ <<endl;
                cout<<"K___ == " << K___ <<endl;
                cout << "epsilon = " << Eps << endl;
                unsigned int seed = rng();
                cout << "random seed = " << seed << endl;

                p = 1.0 / (exp(Eps) + 1.0);

                // estis[iteration] = one_round_biclique_2_3(g, seed);


                estis[iteration] = one_round_biclique_2_K(g, K___, seed); 

                // estis[iteration] = one_round_biclique(g, seed, P___, K___);

                std::cout << "estimate = " << std::fixed << std::setprecision(10) << estis[iteration] << std::endl;

                cout<<"real = "<< real<<endl;

                long double relative_error = abs(estis[iteration] - real) * 1.0 / real;

                cout << "relative error = " << relative_error << endl;
                rel_err.push_back(relative_error);
                cout << endl;

            }

            if(P___ > 2 ){
                // one-round biclique algorithm for general P, Q values
                cout << "Oneround, P___ > 2" << endl;
                cout << "EPS = " << Eps << endl;
                unsigned int seed = rng();
                cout << "random seed = " << seed << endl;

                // estis[iteration] = wedge_based_btf_avg(g, seed);
                // estis[iteration] = VP_wedge_based_two_round_btf(g, seed);

                p = 1.0 / (exp(Eps) + 1.0);

                estis[iteration] = one_round_biclique(g, seed, P___, K___);
                // cout << "estimate = " << estis[iteration] << endl;
                
                std::cout << "estimate = " << std::fixed << std::setprecision(10) << estis[iteration] << std::endl;

                cout<<"real = "<< real<<endl;

                long double relative_error = abs(estis[iteration] - real) * 1.0 / real;

                cout << "relative error = " << relative_error << endl;
                rel_err.push_back(relative_error);
                cout << endl;
            }
        } 
        if (algorithm_switch >= 2) {

            // manipulate swith
            if(algorithm_switch == 2){
                cout << "\nADV" << endl;
                multi_estimator_switch = false; 
                two_noisy_graph_switch = false ;
            }
            if(algorithm_switch == 3){
                cout << "\nADV+" << endl;
                multi_estimator_switch = true; 
                two_noisy_graph_switch = false ;
            }
            if(algorithm_switch == 4){
                cout << "\nADV++" << endl;
                multi_estimator_switch = true; 
                two_noisy_graph_switch = true ;
            }
            
            // in this algorithm, we only use one vertex as the source of estimation.
            cout << "EPS = " << Eps << endl;
            unsigned int seed = rng();
            cout << "random seed = " << seed << endl;
            if(P___ == 2){
                // this will be the basis of the working example
                estis[iteration] = wedge_based_two_round_2_K_biclique(g, seed);

                // layer_based 
                // estis[iteration] = layer_based_wedge_based_two_round_2_K_biclique(g, seed);
            }
            else if(P___ == 3){
                // this is ready.
                estis[iteration] = wedge_based_two_round_3_K_biclique(g, seed);
            }
            // need to implement two_noisy_graph_switch optimization for P in general
            else{
                cout<<"P = "<<P___ <<endl;
                cout<<"Q = "<<K___ <<endl;
                // we do not test the communication of this
                estis[iteration] = wedge_based_two_round_general_biclique(g, seed, P___, K___);            
                // cout<<"Need to implement baseline for general P values"<<endl;
            }
            
            // cout << "estimate = " << estis[iteration] << endl;
            std::cout << "estimate = " << std::fixed << std::setprecision(10) << estis[iteration] << std::endl;
            cout<<"real = "<<real <<endl;
            long double relative_error = abs(estis[iteration] - real) * 1.0 / real;
            cout << "relative error = " << relative_error << endl;
            rel_err.push_back(relative_error);
            cout << endl;
        }
    }

    double t1 = omp_get_wtime();
    double seconds = t1 - t0;
    printf("time:%f\n", seconds);


    if(eva_comm){
        // byte --> metabyte 
        long double communication_cost_MB = communication_cost / (1024.0 * 1024.0);
        printf("com cost (MB) = %Lf\n", communication_cost_MB);
        return 0;
    }
    
    printf("# Mean = %Lf\n", calculateMean(estis));
    cout << "real count = " << real << endl;

    cout << "adv rel err = " << calculateMean(rel_err) << endl;
    return 0;
}