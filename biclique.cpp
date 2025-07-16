#include "biclique.h"
#include "include/mt19937ar.h"

#include <sys/resource.h>
void printMemoryUsage() {

    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "Memory usage: " << usage.ru_maxrss << " KB" << std::endl;
}

using namespace std;
long double _cate, _wedge, _btf;
vector<int> upper_sample, lower_sample;
double sample_ratio = 1.0;
unordered_map<int, bool> in_sampled_up, in_sampled_lo;
bool samling_one_round = false;
long double verified = 0, not_verified = 0;
extern long double Eps, Eps0, Eps1, Eps2, p, m3__, m2__, m1__, m0__, communication_cost;
extern vector<int> priv_deg;
extern vector<long double> naive_estis;
extern int priv_dmax_1, priv_dmax_2;
extern vector<vector<int>> up_options, lo_options;
extern int iteration;

// double noisy graph optimization?
bool two_noisy_graph_switch ; 
bool multi_estimator_switch ; 

extern int K___ ; 

extern unsigned long long int real; 


extern long double avg_estimated_variance ;

bool one_round = false, 
    // multiR= false, wedge_based = false, 
    edge_clipping = true, count_cate = false, 
	sampling_noisy_graph = false, eva_comm = false;

bool count_cc = false ;  // this is computing the bipartite clustering coefficient. 
double p____ = 0.5;  // sampling ratio.
int alpha = 10; 
stats::rand_engine_t engine(std::time(0));  // used to be 1776 

bool vertex_pair_reduction = true; 

// bool averaging_f_estimates = true;
bool averaging_f_estimates = true;  

extern vector<vector<int>> up_options, lo_options; 


extern long double RR_time, server_side_time, naive_server_side, 
	local_count_time, deg_esti_time;


long double gamma__;

bool efficient_RR = true, skip_neg_deg = false ;
// why cannot we skip negative vertices

// convert your BiGraph instance g2 to the biGraph struct 
biGraph convertBiGraphTobiGraph(BiGraph& oldGraph) {
    biGraph newGraph;
    
    // Initialize basic properties
    newGraph.n1 = oldGraph.num_v1;
    newGraph.n2 = oldGraph.num_v2;
    newGraph.m = oldGraph.num_edges;


    // Allocate space
    newGraph.edges.resize(newGraph.m);
    newGraph.e1.resize(newGraph.m);
    newGraph.e2.resize(newGraph.m);
    newGraph.pU.resize(newGraph.n1 + 5);
    newGraph.pV.resize(newGraph.n2 + 5);
    
    int edge_index = 0;
    for (vid_t u = 0; u < oldGraph.num_v1; ++u) {
        // for each upper vertex.
        assert(oldGraph.is_upper(u));

        auto& neighbors = oldGraph.neighbor[u];
        for (const auto& v__ : neighbors) {
            assert(oldGraph.is_lower(v__));
            vid_t v = v__ - oldGraph.num_v1;
            assert(v < oldGraph.num_v2); 

            newGraph.edges[edge_index].u = u;         
            newGraph.edges[edge_index].v = v;          
            edge_index++;
        }
        // Clear memory for this neighbor list
        std::vector<vid_t>().swap(neighbors);
    }

    // cout<<"edge_index = "<<edge_index<<endl;
    assert(edge_index== newGraph.m );

    // compute the degree ordering in new graph
    newGraph.changeToDegreeOrder();

    return newGraph;
}

void randomized_response_single_bit(int u, int v, BiGraph& g, BiGraph& g2) {
	double keep_probability = g.has(u, v) ? 1 - p : p;
	// if(sampling_noisy_graph){
	// 	keep_probability *= p____;
	// }
	assert(u!=v);
	g2.edge_vector[min(u, v)][max(u, v)] = (genrand_real2() < keep_probability);
	// it is either 1 or 0
}

void private_estimate_of_degrees(BiGraph& g) {
    // private estimate degrees.
    priv_dmax_1 = 0;
    priv_dmax_2 = 0;
    priv_deg.resize(g.num_nodes());
    for (int i = 0; i < g.num_nodes(); i++) {
        priv_deg[i] = g.degree[i] + stats::rlaplace(0.0, 1 / (Eps0), engine);
        if (edge_clipping) {
            priv_deg[i] += alpha;
        }
        if (g.is_upper(i)) {
            priv_dmax_1 = priv_dmax_1 > priv_deg[i] ? priv_dmax_1 : priv_deg[i];
        } else {
            priv_dmax_2 = priv_dmax_2 > priv_deg[i] ? priv_dmax_2 : priv_deg[i];
        }
    }
}

double my_genrand_real2() { return genrand_real2(); }

void my_init_genrand(unsigned long seed) { init_genrand(seed); }

void construct_noisy_graph(BiGraph& g, BiGraph& g2, unsigned long seed) {
    const int range_from = g2.num_v1;
    const int range_to = g2.num_nodes() - 1;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<int> distr(range_from, range_to);

    init_genrand(seed);
    int flip1 = 0;
    int visited_vertices = 0;
    int total_vertices = g2.num_v1;
    int ten_percent = total_vertices / 5;
    long double max_time_per_user = -1;

    int max_num_noisy_edges_per_vertex = -1;
    for (int i = 0; i < g2.num_v1; i++) {
        visited_vertices++;
        // if (visited_vertices % ten_percent == 0) {
        // 	int progress = visited_vertices * 100 / total_vertices;
        // 	cout << "Processed " << progress << "% of vertices" << endl;
        // }
        int num_noisy_edges_per_i = 0;
        double tx = omp_get_wtime();
        for (int j = g2.num_v1; j < g2.num_nodes(); j++) {
            if (std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) !=
                g.neighbor[i].end()) {
                if (genrand_real2() >= p) {  // 1  --> 1
                    g2.addEdge(i, j);
                    flip1++;
                    num_noisy_edges_per_i++;
                }
            } else {
                if (genrand_real2() < p) {  // 0 --> 1
                    g2.addEdge(i, j);
                    flip1++;
                    num_noisy_edges_per_i++;
                }
            }
        }
        if(num_noisy_edges_per_i > max_num_noisy_edges_per_vertex){
            max_num_noisy_edges_per_vertex =  num_noisy_edges_per_i; 
        }

        double ty = omp_get_wtime();
        max_time_per_user =max_time_per_user > (ty - tx) ? max_time_per_user : (ty - tx);
    }
    // RR_time += max_time_per_user;
    // the dominating cost is incurred in server side butterfly counting on the
    // dense noisy graph

    cout << "noisy edges = " << flip1 << endl;

    if (eva_comm) {
        cout<<"computing the of Randomized responses (1)"<<endl;
        double byte_per_edge = 8*(log2(g.num_v1) + log2(g.num_v2));
        communication_cost += max_num_noisy_edges_per_vertex * byte_per_edge;
        // communication_cost += flip1 * sizeof(int);
    }

    long double expected_E =g.num_edges * (1 - p) + (g.num_v1 * g.num_v2 - g.num_edges) * p;
    cout << "expected E = " << expected_E << endl;

    g2.computePriority();
}

// applying RR to lower vertices 
void construct_noisy_graph_2(BiGraph& g, BiGraph& g2, unsigned long seed) {
    const int range_from = g2.num_v1;
    const int range_to = g2.num_nodes() - 1;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<int> distr(range_from, range_to);

    init_genrand(seed);
    int flip1 = 0;
    int visited_vertices = 0;
    int total_vertices = g2.num_v1;
    int ten_percent = total_vertices / 5;
    long double max_time_per_user = -1;

    // applying RR on U(G)
    // for (int i = 0; i < g2.num_v1; i++) {

    // now applying RR to L(G)
    int max_num_noisy_edges_per_vertex = -1;
    for (int i = g2.num_v1; i < g2.num_nodes(); i++) {
        visited_vertices++;
        double tx = omp_get_wtime();

        // for (int j = g2.num_v1; j < g2.num_nodes(); j++) {
        int num_noisy_edges_per_i = 0;
        for (int j = 0; j < g2.num_v1; j++) {
            if (std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) !=
                g.neighbor[i].end()) {
                if (genrand_real2() >= p) {  // 1  --> 1
                    g2.addEdge(i, j);
                    flip1++;
                    num_noisy_edges_per_i++;
                }
            } else {
                if (genrand_real2() < p) {  // 0 --> 1
                    g2.addEdge(i, j);
                    flip1++;
                    num_noisy_edges_per_i++;
                }
            }
        }

        if(num_noisy_edges_per_i > max_num_noisy_edges_per_vertex){
            max_num_noisy_edges_per_vertex =  num_noisy_edges_per_i; 
        }
        double ty = omp_get_wtime();
        max_time_per_user =
            max_time_per_user > (ty - tx) ? max_time_per_user : (ty - tx);
    }
    // RR_time += max_time_per_user;
    // the dominating cost is incurred in server side butterfly counting on the
    // dense noisy graph

    cout << "noisy edges = " << flip1 << endl;

    // if (eva_comm) {
    //     cout<<"computing the of Randomized responses (1)"<<endl;
    //     communication_cost += flip1 * sizeof(int);
    // }
    if (eva_comm) {
        cout<<"computing the of Randomized responses (1)"<<endl;
        double byte_per_edge = 8*(log2(g.num_v1) + log2(g.num_v2));
        communication_cost += max_num_noisy_edges_per_vertex * byte_per_edge;
        // communication_cost += flip1 * sizeof(int);
    }
    long double expected_E =g.num_edges * (1 - p) + (g.num_v1 * g.num_v2 - g.num_edges) * p;
    cout << "expected E = " << expected_E << endl;

    g2.computePriority();
}

long double two_round_btf(BiGraph& g, unsigned long seed) {
    // Phase 0. deg_esti_time records the maximum degree perturbation time.
    double t0 = omp_get_wtime();
    cout << "private_estimate_of_degrees(g); " << endl;
    Eps0 = Eps * 0.1;
    private_estimate_of_degrees(g);

    // upload noisy degrees
    if (eva_comm) communication_cost += g.num_nodes() * sizeof(int);

    // Phase 1. RR
    double t1 = omp_get_wtime();
    cout << "construct_noisy_graph(g); " << endl;
    Eps1 = Eps * 0.5;
    p = 1.0 / (exp(Eps1) + 1.0);

    BiGraph g2(g);
    if (sampling_noisy_graph) {
        cout << "sampling ratio = " << p____ << endl;
        const int range_from = g2.num_v1;
        const int range_to = g2.num_nodes() - 1;
        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_int_distribution<int> distr(range_from, range_to);
        init_genrand(seed);
        int flip1 = 0;
        int visited_vertices = 0;
        int total_vertices = g2.num_v1;
        int ten_percent = total_vertices / 5;
        // long double max_time_per_user = -1;
        for (int i = 0; i < g2.num_v1; i++) {
            // double tx = omp_get_wtime();
            for (int j = g2.num_v1; j < g2.num_nodes(); j++) {
                if (std::find(g.neighbor[i].begin(), g.neighbor[i].end(), j) !=
                    g.neighbor[i].end()) {
                    if (genrand_real2() >= p) {  // 1  --> 1
                        if (genrand_real2() >= p____)
                            continue;  // keep with probability p____
                        g2.addEdge(i, j);
                        flip1++;
                    }
                } else {
                    if (genrand_real2() < p) {  // 0 --> 1
                        if (genrand_real2() >= p____)
                            continue;  // keep with probability p____
                        g2.addEdge(i, j);
                        flip1++;
                    }
                }
            }
        }
        cout << "noisy edges = " << flip1 << endl;
        long double expected_E =
            g.num_edges * (1 - p) + (g.num_v1 * g.num_v2 - g.num_edges) * p;
        cout << "expected E = " << expected_E << endl;
        g2.computePriority();  // is this necessary?
    } else {
        construct_noisy_graph(g, g2, seed);  // upload noisy edges
    }

    // Phase 2. local counting records the counting time
    double t2 = omp_get_wtime();
    cout << "local counting" << endl;
    if (eva_comm) {
        // for each vertex, it needs to download all vertex degrees.
        communication_cost += g.num_nodes() * g.num_nodes() * sizeof(int);

        // for each vertex, it needs to the whol noisy graph
        communication_cost += g.num_nodes() * g2.num_edges * sizeof(int);

        // for each vertex, it needs to upload local count wi
        communication_cost +=
            g.num_nodes() * sizeof(long double);  // upload local counts
        return 0;
    }
    Eps2 = Eps * 0.4;

    // use eps2
    long double global_sensitivity, sum = 0;

    for (int u = 0; u < g.num_nodes(); u++) {
        // if(g.degree[u]==0) continue;

        if (edge_clipping && priv_deg[u] <= 0){
            continue;
        }

        // when edge clipping is inplace, we can only visit at most priv_deg[u]
        // neighbors for u
        long double s1 = 0, s2 = 0, s3 = 0;

        unordered_map<vid_t, int> count_wedge(0);

        long double du = g.degree[u];

        if (edge_clipping && (du > priv_deg[u])) {
            du = priv_deg[u];
        }

        s3 += (du * (du - 1) / 2) * ((g.is_upper(u) ? g.num_v1 : g.num_v2) - 1);

        long double sum_deg_v = 0;
        int visited_nb = 0;
        for (auto v : g.neighbor[u]) {
            if (edge_clipping && visited_nb == priv_deg[u]) {
                break;
            }
            sum_deg_v += g2.degree[v];
            for (auto w : g2.neighbor[v]) {
                if (u != w) {
                    count_wedge[w]++;
                }
            }
            if (edge_clipping) visited_nb++;
        }

        for (auto ele : count_wedge) {
            // only execute this for butterfly counting
            if ((!count_cate) && (ele.second >= 2)) {
                s1 += (ele.second - 1) * ele.second / 2; // noisy butterfly counting 
            }
            s2 += ele.second * (du - 1); // noisy caterpillar counting 
        }

        int deg_up = edge_clipping ? priv_deg[u] : g.degree[u];
        if (count_cate) {
            // for caterpillar counting: we only need to count s2, and s3.
            if (g.is_upper(u)) {
                global_sensitivity = 3 * deg_up * g2.v2_max_degree + sum_deg_v +
                                     2 * p * deg_up * (g.num_v1 - 1);
            } else {
                global_sensitivity = 3 * deg_up * g2.v1_max_degree + sum_deg_v +
                                     2 * p * deg_up * (g.num_v2 - 1);
            }
            sum += s2 - 2 * p * s3;
        } else {
            // for buttrfly counting:
            if (g.is_upper(u)) {
                global_sensitivity = (1 - p) * deg_up * g2.v2_max_degree +
                                     p * sum_deg_v +
                                     p * p * deg_up * (g.num_v1 - 1);
            } else {
                global_sensitivity = (1 - p) * deg_up * g2.v1_max_degree +
                                     p * sum_deg_v +
                                     p * p * deg_up * (g.num_v2 - 1);
            }

            // question: why cannot we make GS smaller? 
            // i.e., g2.v1_max_degree seems pretty large 
            sum += s1 - p * s2 + p * p * s3;
        }

        // I  think it is effectively counting the number of butterflies containing u. 
        sum += stats::rlaplace(0.0, (global_sensitivity / Eps2), engine);  // add calibrated noise
        // communication_cost += sizeof(long double);
    }
    // return sum/(4*(1-2*p)*(1-2*p));
    double t3 = omp_get_wtime();

    RR_time += t2 - t1;
    deg_esti_time += t1 - t0;
    local_count_time += t3 - t2;
    
    if (count_cate) {
        if (sampling_noisy_graph) sum /= p____;
        return sum / (2 * (1 - 2 * p));
    } else {
        if (sampling_noisy_graph){
            sum /= p____ * p____;
        }
        return sum / (4 * (1 - 2 * p) * (1 - 2 * p));
    }
}

void compute_m3_m2(long double& m4, long double& m3, long double& m2,
                   long double& m1, long double& m0, BiGraph& g2) {
    long double caterpillar = 0;
    long double chopsticks = 0;
    long double wedges = 0;

    long double g2_num_edges = g2.num_edges;

    long double g2_num_v1 = g2.num_v1;
    long double g2_num_v2 = g2.num_v2;

#pragma omp parallel for reduction(+ : caterpillar, chopsticks)
    for (int i = 0; i < g2.num_v1; i++) {
        for (auto j : g2.neighbor[i]) {
            // for each edge (i,j)
            caterpillar += (g2.degree[i] - 1) * (g2.degree[j] - 1);
            chopsticks += g2_num_edges - g2.degree[i] - g2.degree[j] + 1;
        }
    }

    chopsticks /= 2;  // each is counted twice

#pragma omp parallel for reduction(+ : wedges)
    for (int i = 0; i < g2.num_nodes(); i++) {
        long double deg_i = g2.degree[i];
        if (deg_i == 0) continue;
        if (g2.is_upper(i)) {
            wedges += (g2_num_v1 - 1) * deg_i * (deg_i - 1) /
                      2;  // is this correct? I am curious.
        } else {
            wedges += (g2_num_v2 - 1) * deg_i * (deg_i - 1) / 2;
        }
    }

    m3 = caterpillar - 4 * m4;

    long double m21 = wedges - 4 * m4 - 2 * m3;

    long double m22 = chopsticks - 2 * m4 - m3;

    m2 = m21 + m22;

    m1 = g2_num_edges * (g2_num_v1 - 1) * (g2_num_v2 - 1) - 4 * m4 - 3 * m3 -
         2 * m2;

    m0 = (g2_num_v1 * (g2_num_v1 - 1) / 2) * (g2_num_v2 * (g2_num_v2 - 1) / 2) -
         m4 - m3 - m2 - m1;

    cout << "m4 = " << m4 << endl;
    cout << "m3 = " << m3 << endl;
    cout << "m2 = " << m2 << endl;
    cout << "m1 = " << m1 << endl;
    cout << "m0 = " << m0 << endl;
}

long double one_round_btf(BiGraph& g, unsigned long seed) {
    long double t0 = omp_get_wtime();

    std::mt19937 rng(std::random_device{}());

    BiGraph g2(g);

    // user side:
    construct_noisy_graph(g, g2, seed);

    long double t1 = omp_get_wtime();

    RR_time += t1 - t0;

    // server side:
    long double m4;

    if (sampling_noisy_graph) {
        cout << "sampling ratio = " << p____ << endl;
        long double sum__ = 0;
        int num_itr = 1;
        for (int xxx = 0; xxx < num_itr; xxx++) {
            // get a sampled subgraph:
            BiGraph g__(g);  // sampled subgraph
            init_genrand(rng());
            for (int i = 0; i < g2.num_v1; i++) {
                for (auto j : g2.neighbor[i]) {
                    if (genrand_real2() < p____) {
                        g__.addEdge(i, j);
                    }
                }
            }
            g__.computePriority();
            cout << "num edges in sampled noisy graph = " << g__.num_edges
                 << endl;
            m4 = BFC_EVP(g__);
            m4 /= (p____ * p____ * p____ * p____);
            sum__ += m4;
        }
        m4 = sum__ / num_itr;

        cout << "estimated m4' = " << m4 << endl;

    } else {
        cout<<"computing m4"<<endl;
        m4 = BFC_EVP(g2);
    }

    cout << "m4 is ready" << endl;

    long double t1x = omp_get_wtime();

    // should we estimate cate and other things independently?

    long double m3 = 0, m2 = 0, m1 = 0, m0 = 0, estimate, mu = exp(Eps);

    compute_m3_m2(m4, m3, m2, m1, m0, g2);

    cout << "computing m3, m2, m1, m0" << endl;

    // this actually means getting bipartite clustering coefficient
    /*
    if(count_cc){
        long double cate_estimate = -4 * mu * mu * mu * m4 + mu * mu * (mu * mu + 3) * m3 -
                   2 * mu * (mu * mu + 1) * m2 + (3 * mu * mu + 1) * m1 -
                   4 * mu * m0;
        cate_estimate /= ((mu - 1) * (mu - 1) * (mu - 1) * (mu - 1));
        // 
        long double btf_estimate = mu * mu * mu * mu * m4 - mu * mu * mu * m3 + mu * mu * m2 -
                   mu * m1 + m0;
        btf_estimate /= ((mu - 1) * (mu - 1) * (mu - 1) * (mu - 1));
        estimate = 4 * btf_estimate / cate_estimate; 

        naive_estis[iteration] = m4 * 4/(m4 * 4 + m3);


    }else 
    */
    if (count_cate) {
        cout << "\tcaterpillar estimation: " << endl;
        estimate = -4 * mu * mu * mu * m4 
                   + mu * mu * (mu * mu + 3) * m3 -
                   2 * mu * (mu * mu + 1) * m2 
                   + (3 * mu * mu + 1) * m1 
                   -4 * mu * m0;
        estimate /= ((mu - 1) * (mu - 1) * (mu - 1) * (mu - 1));
        naive_estis[iteration] = m4 * 4 + m3;
    } 
    else {
        cout << "\tbtf estimation: " << endl;
        estimate = mu * mu * mu * mu * m4 - mu * mu * mu * m3 + mu * mu * m2 -
                   mu * m1 + m0;
        estimate /= ((mu - 1) * (mu - 1) * (mu - 1) * (mu - 1));
        naive_estis[iteration] = m4;
        // the time needed to get m4 on G2 alone.
        naive_server_side += t1x - t1;
    }
    long double t2 = omp_get_wtime();
    // record time elapsed.
    // RR_time += t1-t0; // compute the total time from randomized response.
    server_side_time += t2 - t1;

    return estimate;
}

// the challenge lies in how to compute deg(u, w) = {v \in N(u), v < w } for
// each u < w combination.
long double BFC_EVP(BiGraph& g) {
    long double BTF = 0;
#pragma omp parallel for reduction(+ : BTF)
    for (int u = 0; u < g.num_nodes(); u++) {
        if (g.degree[u] <= 1) continue;

        // cout<<"u  = "<<u<<endl;
        unordered_map<vid_t, int> count_wedge(0);
        for (auto v : g.neighbor[u]) {
            // cout<<"\t v = "<<v<<endl;
            // u -> v
            for (auto w : g.neighbor[v]) {
                // u->v->w
                // cout<<"u, v, w= "<<u<<" "<<"v"<<" "<<w<<endl;
                // this step is not invoked for g3.
                // what if we just use id.
                if (g.com_p(w, v) & g.com_p(w, u)) {  // this is a lot faster.
                    count_wedge[w] = count_wedge[w] + 1;
                }
            }
        }
        // long double btf_u = 0;
        for (auto ele : count_wedge) {
            // cout<<"wedge count = " <<ele.second <<endl;
            if (ele.second >= 2) {
                BTF += (ele.second - 1) * ele.second / 2;
            }
        }
    }

    // cout<<"BTF = "<<BTF<<endl;
    return BTF;
}

long double get_wedges(BiGraph& g) {
    long double wedges = 0;
    for (int i = 0; i < g.num_nodes(); i++) {
        long double deg_i = g.degree[i];
        if (deg_i == 0) continue;
        if (g.is_upper(i)) {
            wedges += deg_i * (deg_i - 1) / 2;
        } else {
            wedges += deg_i * (deg_i - 1) / 2;
        }
    }
    return wedges;
}

long double get_laplace(long double parameter) {
    return stats::rlaplace(0.0, parameter, engine);
}

long double get_cate(BiGraph& g) {
    long double caterpillar = 0;
    for (int i = 0; i < g.num_v1; i++) {
        for (auto j : g.neighbor[i]) {
            caterpillar += (g.degree[i] - 1) * (g.degree[j] - 1);
        }
    }
    return caterpillar;
}

// new approach: 
// we might also need to adopt some budget optimization strategy here. 
// this is a two phase algorithm. 
// first construct the whole noisy graph, and then use it
// is it possible to combine this with the two-round algorithm? 


long double wedge_based_two_round_2_K_biclique(BiGraph& g, unsigned long seed) {
    // Phase 0. deg_esti_time records the maximum degree perturbation time.
    // double t0 = omp_get_wtime();
    // cout << "private_estimate_of_degrees(g); " << endl;
    Eps0 = Eps * 0.05;
    // private_estimate_of_degrees(g);
	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
	for(int i=0;i<g.num_nodes();i++){
		deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
	}
    // if (eva_comm) {
    //     cout<<"communication cots of uploading degree estimates"<<endl;
    //     communication_cost += g.num_nodes() * sizeof(int);
    // }

    // upload noisy degrees

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Phase 1. RR
    double t1 = omp_get_wtime();
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);
    cout<<"constructing g2\n";
	construct_noisy_graph(g, g2, seed);  // upload noisy edges
    // unfortunately, this step cannot be run in parallel


    BiGraph g3(g);
    if(two_noisy_graph_switch){
        cout<<"constructing g3\n";
        construct_noisy_graph_2(g, g3, seed);  // upload noisy edges
    }

    // Phase 2. local counting
    double t2 = omp_get_wtime();
    cout << "local counting" << endl;

    if (eva_comm) {
        // for each vertex, download all vertex degrees.
        // communication_cost += g.num_nodes() * g.num_nodes() * sizeof(int);



        // for each vertex, download the whole noisy graph
        double byte_per_edge = 8*(log2(g.num_v1) + log2(g.num_v2));
        communication_cost += g2.num_edges * byte_per_edge;
        if(two_noisy_graph_switch){
            communication_cost += g3.num_edges * byte_per_edge;
        }

        // upload common neighbor estimates:
        size_t pairwise_count = g.num_v1 * (g.num_v1 - 1) / 2;
        communication_cost += sizeof(long double) * pairwise_count;
        if (multi_estimator_switch) {
            communication_cost += sizeof(long double) * pairwise_count;
        }

        return 0;
    }

	Eps2 = Eps - Eps1 - Eps0;
    
	// cout<<"using Eps2 = "<<Eps2 <<endl;
	gamma__ = (1-p) / (1-2*p);
    // remember to use eps2
    // long double global_sensitivity, sum = 0;
	long double res___ = 0; 

	// what if we only consider upper vertices ?  --> better efficiency and effect  
	int K = K___;  // we are considering (2, K)-biclique right now

    cout<<"K___ = "<<K___ <<endl;

	int start__, end__;

	start__ = g.num_v1 < g.num_v2 ? 0 : g.num_v1; 
	end__ = g.num_v1 < g.num_v2 ? g.num_v1 : g.num_nodes(); 

	#pragma omp parallel
	{
	#pragma omp for schedule(static)
		for(int u =start__ ; u <end__ ; u++) {
			for(int w =start__ ; w <end__ ; w++) {
                if(u<=w) continue;

                long double f_u_w, f_w_u;    
                if(two_noisy_graph_switch){
                    // when this switch is on, by default we expect multiple estimators.
                    f_u_w = locally_compute_f_given_q_and_x_two_graphs(u, w, g, g2, g3);
                    f_w_u = locally_compute_f_given_q_and_x_two_graphs(w, u, g, g2, g3);

                    // basically getting the same thing using two noisy graphs.
                }else{
                    f_u_w = locally_compute_f_given_q_and_x(u, w, g, g2);
                    if(multi_estimator_switch){
                        f_w_u = locally_compute_f_given_q_and_x(w, u, g, g2);
                    }
                }

                long double diff1 =0, diff2 = 0;
                
                long double local_res = 0;

                /*
                bool new_method = false; 
                if(new_method){
                    // averaging the biclique estimators (average in the end)
                    long double variance_f_u = 2 * pow(gamma__, 2) / pow(Eps2, 2) + 
                                            p * (1 - p) * deg_estis[u] / pow(1 - 2 * p, 2);

                    long double variance_f_w = 2 * pow(gamma__, 2) / pow(Eps2, 2) + 
                                            p * (1 - p) * deg_estis[w] / pow(1 - 2 * p, 2);

                    long double moment_2 = pow(f_u_w, 2) - variance_f_u; // f^2

                    long double local_res_1 =  (moment_2 - f_u_w)/2 ;


                    long double moment_2_ = pow(f_w_u, 2) - variance_f_w; // f^2

                    long double local_res_2 =  (moment_2_ - f_w_u)/2 ;  

                    local_res =    (local_res_1 + local_res_2)/2  ;

                }else{
                */

                // define some variables
                long double esti_var_f, variance_f_u, variance_f_w, main_fu, main_fw;

                if(!multi_estimator_switch){
                    // single source estimator: 
                    esti_var_f = 2 * pow(gamma__,2)  / pow(Eps2,2) + p * (1 - p) * deg_estis[u] / pow(1-2*p,2); 
                }else{
                    // multi source estimator:  
                    f_u_w = (f_u_w + f_w_u)/2;
                    // esitmate the variance of f_u_w
                    esti_var_f = 0;
                    // these are variance from laplace, not affected by two_noisy_graph_switch
                    variance_f_u = 2 * pow(gamma__,2)  / pow(Eps2,2);
                    variance_f_w = 2 * pow(gamma__,2)  / pow(Eps2,2);
                    // main_fu and main_fw are the variance from local counts
                    main_fu =  p * (1 - p) * deg_estis[u] / pow(1-2*p,2);
                    main_fw =  p * (1 - p) * deg_estis[w] / pow(1-2*p,2);
                    if(two_noisy_graph_switch){
                        // maybe we should increase epsilon 2 to reduce the impact of laplace.
                        main_fu/=2;
                        main_fw/=2;
                    }
                    variance_f_u += main_fu;
                    variance_f_w += main_fw;
                    esti_var_f = (variance_f_u + variance_f_w) / 4;
                }
                // (2, K)-biclique need these moments of the unbiased estimate of f^2:
                local_res = compute_local_res(K, f_u_w, esti_var_f);

				#pragma omp critical
				res___ += local_res; 
			}

		}
	}
    return res___;
}


// the new version of the wedge_based_two_round_2_K_biclique algorithm 
// with layer-based optimization 
long double layer_based_wedge_based_two_round_2_K_biclique(BiGraph& g, unsigned long seed) {

    // estimated dmax1 and dmax2
    long double dmax1 = 0, dmax2 = 0;

    long double real_dmax1 = 0, real_dmax2 = 0;

    Eps0 = Eps * 0.05;
    // private_estimate_of_degrees(g);
	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());

    for(int i = 0; i < g.num_nodes(); i++) {
        deg_estis[i] = g.degree[i] + stats::rlaplace(0.0, 1/(Eps0), engine); 
        dmax1 = g.is_upper(i) ? (deg_estis[i] > dmax1 ? deg_estis[i] : dmax1) : dmax1;
        real_dmax1 = g.is_upper(i) ? (static_cast<long double>(g.degree[i]) > real_dmax1 ? static_cast<long double>(g.degree[i]) : real_dmax1) : real_dmax1;
        dmax2 = !g.is_upper(i) ? (deg_estis[i] > dmax2 ? deg_estis[i] : dmax2) : dmax2;
        real_dmax2 = !g.is_upper(i) ? (static_cast<long double>(g.degree[i]) > real_dmax2 ? static_cast<long double>(g.degree[i]) : real_dmax2) : real_dmax2;
    }

    cout<<"dmax1 = "<<real_dmax1<<", esti = "<<dmax1 <<endl;
    cout<<"dmax2 = "<<real_dmax2<<", esti = "<<dmax2 <<endl;
    // exit(1);


    cout<<"estimate the cost associatd with p-tuple enumeration in U(G) and q-tuple eumeration in L(G)"<<endl;
    // Calculate S1 and S2
    // Assuming n1, n2, p, q are defined (e.g., n1 and n2 are sizes of partitions)
    long double S1 = binomial(g.num_v1, 2) * dmax1;
    long double S2 = binomial(g.num_v2, K___) * dmax2;
    cout<<"S1 = "<<S1 <<endl;
    cout<<"S2 = "<<S2 <<endl;

    bool use_upper_layer_for_enumeration = false; 
    if(S1 < S2){
        use_upper_layer_for_enumeration = true; 
    }


    // right now, it simply look at which layer is smaller.
	int start__, end__;
	start__ = g.num_v1 < g.num_v2 ? 0 : g.num_v1; 
	end__ = g.num_v1 < g.num_v2 ? g.num_v1 : g.num_nodes(); 

    cout<<"start__ = "<<start__ <<endl;
    cout<<"end__ = "<<end__ <<endl;

    if(start__==0){
        cout<<"n1 n2 tells me to use U(G)"<<endl;
    }else{
        cout<<"n1 n2 tells me to use L(G)"<<endl;
    }


    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Phase 1. RR
    double t1 = omp_get_wtime();
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);
    cout<<"constructing g2\n";
	construct_noisy_graph(g, g2, seed);  // upload noisy edges
    // unfortunately, this step cannot be run in parallel





    BiGraph g3(g);
    if(two_noisy_graph_switch){
        cout<<"constructing g3\n";
        construct_noisy_graph_2(g, g3, seed);  // upload noisy edges
    }

    // Phase 2. local counting
    double t2 = omp_get_wtime();
    cout << "local counting" << endl;

	Eps2 = Eps - Eps1 - Eps0;
    
	// cout<<"using Eps2 = "<<Eps2 <<endl;
	gamma__ = (1-p) / (1-2*p);

	long double res___ = 0; 

	// what if we only consider upper vertices ?  --> better efficiency and effect  
	int K = K___;  // we are considering (2, K)-biclique right now

    cout<<"K___ (q) = "<<K___ <<endl;






    // if()



	#pragma omp parallel
	{
	#pragma omp for schedule(static)
		for(int u =start__ ; u <end__ ; u++) {
			for(int w =start__ ; w <end__ ; w++) {
                if(u<=w) continue;

                long double f_u_w, f_w_u;    
                if(two_noisy_graph_switch){
                    f_u_w = locally_compute_f_given_q_and_x_two_graphs(u, w, g, g2, g3);
                    f_w_u = locally_compute_f_given_q_and_x_two_graphs(w, u, g, g2, g3);

                    // basically getting the same thing using two noisy graphs.
                }else{
                    f_u_w = locally_compute_f_given_q_and_x(u, w, g, g2);
                    if(multi_estimator_switch){
                        f_w_u = locally_compute_f_given_q_and_x(w, u, g, g2);
                    }
                }

                long double diff1 =0, diff2 = 0;
                
                long double local_res = 0;

                // define some variables
                long double esti_var_f, variance_f_u, variance_f_w, main_fu, main_fw;

                if(!multi_estimator_switch){
                    // single source estimator: 
                    esti_var_f = 2 * pow(gamma__,2)  / pow(Eps2,2) + p * (1 - p) * deg_estis[u] / pow(1-2*p,2); 
                }else{
                    // multi source estimator:  
                    f_u_w = (f_u_w + f_w_u)/2;
                    // esitmate the variance of f_u_w
                    esti_var_f = 0;
                    // these are variance from laplace, not affected by two_noisy_graph_switch
                    variance_f_u = 2 * pow(gamma__,2)  / pow(Eps2,2);
                    variance_f_w = 2 * pow(gamma__,2)  / pow(Eps2,2);
                    // main_fu and main_fw are the variance from local counts
                    main_fu =  p * (1 - p) * deg_estis[u] / pow(1-2*p,2);
                    main_fw =  p * (1 - p) * deg_estis[w] / pow(1-2*p,2);
                    if(two_noisy_graph_switch){
                        // maybe we should increase epsilon 2 to reduce the impact of laplace.
                        main_fu/=2;
                        main_fw/=2;
                    }
                    variance_f_u += main_fu;
                    variance_f_w += main_fw;
                    esti_var_f = (variance_f_u + variance_f_w) / 4;
                }
                // (2, K)-biclique need these moments of the unbiased estimate of f^2:
                local_res = compute_local_res(K, f_u_w, esti_var_f);

				#pragma omp critical
				res___ += local_res; 
			}

		}
	}
    return res___;
}

// this function is here to handle when p = 3
long double wedge_based_two_round_3_K_biclique(BiGraph& g, unsigned long seed) {
    double t1 = omp_get_wtime();

    Eps0 = Eps * 0.05;

	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
	for(int i=0;i<g.num_nodes();i++){
		deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
	}

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Phase 1. RR
    
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);
	construct_noisy_graph(g, g2, seed);  


    // two noisy graph technique
    BiGraph g3(g);
    if(two_noisy_graph_switch){
        cout<<"constructing g3\n";
        construct_noisy_graph_2(g, g3, seed);  // upload noisy edges
    }

	Eps2 = Eps - Eps1 - Eps0;
    
	long double res___ = 0; 


	int K = K___;  // we are considering (2, K)-biclique right now

    cout<<"p = "<<3 <<endl;
    cout<<"q = "<<K___ <<endl;


    // Calculate the size of the smaller partition
    int smaller_partition_size = std::min(g.num_v1, g.num_v2);

    // Calculate the total number of possible triples in the smaller partition
    long long total_triples = static_cast<long long>(smaller_partition_size) * 
                            (smaller_partition_size - 1) * 
                            (smaller_partition_size - 2) / 6; 


    // Determine how many triples to sample based on the fraction
    // double sample_fraction = 1e-4;

    // we can increase this for better effectiveness
    long double T = pow(10,6);

    double sample_fraction = T / total_triples;
    
    long long num_triples_to_sample = static_cast<long long>(total_triples * sample_fraction);

    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Sample the triples
    cout<<"total triplets: "<<total_triples <<endl;
    cout << "Sample size = " << num_triples_to_sample <<endl;

    bool is_upper_smaller = (g.num_v1 < g.num_v2 );
    gamma__ = (1-p) / (1-2*p);

    long double res = 0; 

    // Determine the start and end bounds based on partition sizes
    // by default, right now we only consider upper layer.
    int start__ = 0; 
    int end__ = g.num_v1; 

    #pragma omp parallel
	{
        // #pragma omp for schedule(static)
        // for (size_t i = 0; i < tripletVec.size(); ++i) {
        //     const auto& triplet = tripletVec[i];
        //     auto [v1, v2, v3] = triplet;
    #pragma omp for schedule(static)
    for (int v1 = start__; v1 < end__ - 2; ++v1) {
        for (int v2 = v1 + 1; v2 < end__ - 1; ++v2) {
            for (int v3 = v2 + 1; v3 < end__; ++v3) {

                if (dis(gen) >= sample_fraction) continue;


                // from the neighbors of v1:
                long double f1 = 0, f2= 0, f3 = 0, f12 = 0, f13=0;
                long double local_res = 0, esti_var_f_uvw = 0, fuvw = 0 ;

                if(multi_estimator_switch){
                    // multi-source estimator
                    for(auto nb: g.neighbor[v1]){
                        long double A1 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                        long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);

                        if (two_noisy_graph_switch) {
                            A1 = (A1 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                            A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                        }
                        f1 += A1 * A2; 
                        f12 += A1; 
                        f13 += A2; 
                    }
                    // two_noisy_graph_switch does not change GS and Lap noise
                    f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
                    f12 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    f13 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    
                    long double f21 = 0, f23=0;
                    for(auto nb: g.neighbor[v2]){
                        long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                        long double A2 = (static_cast<long double>(g2.has(nb, v3)) - p) / (1 - 2 * p);
                        if (two_noisy_graph_switch) {
                            A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                            A2 = (A2 + (static_cast<long double>(g3.has(nb, v3)) - p) / (1 - 2 * p)) / 2;
                        }
                        f2 += A1 * A2; 
                        f21 += A1; 
                        f23 += A2;          
                    }
                    // two_noisy_graph_switch does not change GS and Lap noise
                    f2 += stats::rlaplace(0.0,  (gamma__*gamma__/Eps2), engine); 
                    f21 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    f23 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    long double f31 = 0, f32=0;
                    for(auto nb: g.neighbor[v3]){
                        long double A1 = (static_cast<long double>(g2.has(nb, v1)) - p) / (1 - 2 * p);
                        long double A2 = (static_cast<long double>(g2.has(nb, v2)) - p) / (1 - 2 * p);
                        if (two_noisy_graph_switch) {
                            A1 = (A1 + (static_cast<long double>(g3.has(nb, v1)) - p) / (1 - 2 * p)) / 2;
                            A2 = (A2 + (static_cast<long double>(g3.has(nb, v2)) - p) / (1 - 2 * p)) / 2;
                        }
                        f3 += A1 * A2; 
                        f31 += A1; 
                        f32 += A2;    

                    }
                    f3 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 
                    f31 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 
                    f32 += stats::rlaplace(0.0, (gamma__/Eps2), engine); 

                    // averaging
                    fuvw = (f1 + f2 + f3 )/3 ; 

                    long double var_phi = p * (1-p) / pow(1-2*p, 2); 
                    if (two_noisy_graph_switch) {
                        var_phi/=2;
                    }

                    long double esti_var_f1, esti_var_f2, esti_var_f3; 

                    // this should always be true
                    bool improvement = true;
                    if(improvement){
                        // averaging f12 and f21 is useful too!
                        esti_var_f1 = var_phi * (f12 + f21 + f13 + f31)/2 ; 
                        esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                        esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2); // lap noise

                        esti_var_f2 = var_phi * (f12 + f21 + f23 + f32)/2 ; 
                        esti_var_f2 += deg_estis[v2] * pow(var_phi,2);
                        esti_var_f2 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                        esti_var_f3 = var_phi * (f13 + f31 + f23 + f32)/2 ; 
                        esti_var_f3 += deg_estis[v3] * pow(var_phi,2);
                        esti_var_f3 += 2 * pow(gamma__,4) / pow(Eps2, 2);
                    }else{ 
                        // this will be true when switch = 3
                        esti_var_f1 = var_phi * (f12 + f13) ; 
                        esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                        esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2); // lap noise

                        esti_var_f2 = var_phi * (f12 + f23) ; 
                        esti_var_f2 += deg_estis[v2] * pow(var_phi,2);
                        esti_var_f2 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                        esti_var_f3 = var_phi * (f13 + f23) ; 
                        esti_var_f3 += deg_estis[v3] * pow(var_phi,2);
                        esti_var_f3 += 2 * pow(gamma__,4) / pow(Eps2, 2);
                    }

                    
                    // include C2(vi, vj)
                    f12 = (f12 + f21)/2;
                    f13 = (f13 + f31)/2;
                    f23 = (f23 + f32)/2;

                    // compute terms from var(f1), var(f2), and var(f3)
                    esti_var_f_uvw = (esti_var_f1 + esti_var_f2 + esti_var_f3);

                    // consider the co-variance of f1, f2, and f3.
                    esti_var_f_uvw += 2 * var_phi * (f12 + f13 + f23); 
                    esti_var_f_uvw /= 9;
                }else{
                    // single-source estimator: f1
                    for(auto nb: g.neighbor[v1]){
                        long double A1 = g2.has(nb, v2) ? 1 : 0 ; 
                        A1 = (A1-p) / (1-2*p); 

                        long double A2 = g2.has(nb, v3) ? 1 : 0 ; 
                        A2 = (A2-p) / (1-2*p); 

                        f1 += A1 * A2;
                    }
                    f1 += stats::rlaplace(0.0, (gamma__*gamma__/Eps2), engine); 

                    // no averaging
                    fuvw = f1;
                    
                    // estimate the variance of f(u, v, w)

                    long double var_phi = p * (1-p) / pow(1-2*p, 2); 

                    long double esti_var_f1 = var_phi * (f12 + f13) ; 
                    esti_var_f1 += deg_estis[v1] * pow(var_phi,2);
                    esti_var_f1 += 2 * pow(gamma__,4) / pow(Eps2, 2);

                    esti_var_f_uvw = esti_var_f1;
                }

                local_res = compute_local_res(K, fuvw, esti_var_f_uvw);

                #pragma omp critical
                res += local_res; 
            }
        }
    }
        
    }

    // assert(count == num_triples_to_sample);
    double t2 = omp_get_wtime();
    cout<<"time  = "<<t2 - t1 <<endl;
    return res / sample_fraction;

}

// need to implement the version for P in general. 
// let's first test whether this works for P = 2 and 3.
long double wedge_based_two_round_general_biclique(BiGraph& g, 
    unsigned long seed, int P___, int K___ ) {
    double t1 = omp_get_wtime();

    Eps0 = Eps * 0.05;

	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
	for(int i=0;i<g.num_nodes();i++){
		deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
	}

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    // Phase 1. RR
    
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);
    gamma__ = (1-p) / (1-2*p);
    BiGraph g2(g);
	construct_noisy_graph(g, g2, seed);  

	Eps2 = Eps - Eps1 - Eps0;
    
	long double res___ = 0; 

    int N = g.num_v1 ; 

    // prepare the subsets in advance
    vector<vector<int>> subsets;
    vector<int> subset(P___);
    for (int i = 0; i < P___; ++i) subset[i] = i;
    // generating all combinations of size P___ from N elements. 
    while (true) {
        subsets.push_back(subset); // Store the current subset
        int i;
        for (i = P___ - 1; i >= 0; --i) {
            if (subset[i] != i + N - P___) {
                ++subset[i];
                for (int j = i + 1; j < P___; ++j) {
                    subset[j] = subset[j - 1] + 1;
                }
                break;
            }
        }
        if (i < 0) break; // Finished all combinations
    }



    #pragma omp parallel
	{
    #pragma omp for schedule(static)
    for (int i__ = 0;  i__ < subsets.size(); ++i__) {
        const auto& subset = subsets[i__];


        // picking the min degree vertex as source
        int v1 = subset[0];
        long double min_deg = deg_estis[v1];
        for (int i = 1; i < P___; ++i) {
            if (deg_estis[subset[i]] < min_deg) {
                min_deg = deg_estis[subset[i]];
                v1 = subset[i];
            }
        }
        
        // std::cout << "Min deg vertex (v1): " << v1 << ", deg_est: " << min_deg << "\n";

        // given a central vertex v1 and construct its sup set: X

        // instead loop through vertices in subset
        // long double aggregated_local_res = 0;
        // for(auto v1:subset){
            // if(v1!=subset[0]) continue;

            // construct the X set based on v1 (everything minus v1)
            vector<int> X;
            for (int i = 0; i < P___; ++i) {
                if (subset[i] != v1){
                    X.emplace_back(subset[i]);
                }
            }
            // construct estimators using v1
            // step 1: estimate the number of common neighbors among all vertices in subset
            long double f1 = 0;
            // long double real_f1 = 0; 
            for(auto nb: g.neighbor[v1]){
                long double fv1 = 1.0; 

                // computing a product
                for(auto xx : X){
                    long double A_ = g2.has(nb, xx) ? 1 : 0 ; 
                    A_ = (A_-p) / (1-2*p); 
                    fv1 = fv1 * A_;
                    // real_fv1 = real_fv1 * (g.has(nb, xx) ? 1 : 0 );
                }
                f1 += fv1; 
                // real_f1 += real_fv1;
            }
            f1 += stats::rlaplace(0.0, (pow(gamma__,X.size())/Eps2), engine); 
            
            
            // step 2: estimate the variance of f1.
            long double esti_var_f = 0; 
            long double theta = p * (1-p) / pow(1-2*p, 2); 
            // X is the set excluding \ v1.

            long double esti_var_f_noisy = 0; 
            int X_size = X.size();
            // for (int mask = 0; mask < (1 << X_size); ++mask) { 
            for (int mask = 0; mask < (1 << X_size) - 1; ++mask) {
                // Include all subsets, from 0 to (1 << X_size) - 1
                vector<int> Y;
                for (int i = 0; i < X_size; ++i) {
                    if (mask & (1 << i)) {
                        Y.push_back(X[i]);
                    }
                }
                // processing subset Y:
                // estimate the number of common neighbors among Y \cup v1.
                long double f1Y = 0; 
                if(Y.size()==0){
                    f1Y = 1; 
                }else{
                    // what if I make f1Y deterministic. 
                    for(auto nb: g.neighbor[v1]){
                        long double fv1 = 1.0; 
                        // long double real_fv1 = 1.0; 
                        for(auto xx : Y){
                            long double A_ = g2.has(nb, xx) ? 1 : 0 ; 
                            A_ = (A_-p) / (1-2*p); 
                            fv1 = fv1 * A_; 

                            // real_fv1 *= (g.has(nb, xx) ? 1 : 0);
                        }
                        f1Y += fv1; 
                        // f1Y += real_fv1; 
                    }
                    // f1Y += stats::rlaplace(0.0, (pow(gamma__,Y.size())/Eps2), engine); 
                }
                esti_var_f_noisy += pow(theta, P___-1-Y.size()) * f1Y;

            }
            esti_var_f_noisy += 2 *pow(gamma__,2*P___-2) / pow(Eps2, 2);
            esti_var_f = esti_var_f_noisy;






            // to do: need to debug why variance estimation is off.

            // naive_estis[iteration] = esti_var_f ;
            naive_estis[iteration] = esti_var_f_noisy ;

            // special handling

            long double local_res = compute_local_res(K___, f1, esti_var_f);

        #pragma omp critical
        res___ += local_res;
        // res___ += aggregated_local_res;
    }
    }

    return res___ ;
}

long double compute_local_res(int K, long double f_u_w, long double esti_var_f) {
    std::vector<long double> moment(K + 1, 0);
    moment[1] = f_u_w;
    moment[2] = pow(f_u_w, 2) - esti_var_f;

    if (K == 2) {
        return (moment[2] - f_u_w) / 2;
    }

    if (K >= 3) {
        moment[3] = pow(f_u_w, 3) - 3 * f_u_w * esti_var_f;
        if (K == 3) {
            return (moment[3] - 3 * moment[2] + 2 * f_u_w) / 6;
        }
    }

    if (K >= 4) {
        moment[4] = pow(f_u_w, 4)
                    - 6 * moment[2] * esti_var_f
                    - 3 * pow(esti_var_f, 2);
        if (K == 4) {
            return (moment[4] - 6 * moment[3] + 11 * moment[2] - 6 * f_u_w) / 24;
        }
    }

    if (K >= 5) {
        moment[5] = pow(f_u_w, 5)
                    - 10 * moment[3] * esti_var_f
                    - 15 * f_u_w * pow(esti_var_f, 2);
        if (K == 5) {
            return (moment[5] - 10 * moment[4] + 35 * moment[3] - 50 * moment[2] + 24 * f_u_w) / 120;
        }
    }

    if (K >= 6) {
        moment[6] = pow(f_u_w, 6)
                    - 15 * moment[4] * esti_var_f
                    - 45 * moment[2] * pow(esti_var_f, 2)
                    - 15 * pow(esti_var_f, 3);
        if (K == 6) {
            return (moment[6] - 15 * moment[5] + 85 * moment[4] - 225 * moment[3]
                    + 274 * moment[2] - 120 * f_u_w) / 720;
        }
    }

    if (K >= 7) {
        moment[7] = pow(f_u_w, 7)
                    + 21 * moment[5] * esti_var_f
                    + 105 * moment[3] * pow(esti_var_f, 2)
                    + 105 * f_u_w * pow(esti_var_f, 3);
        if (K == 7) {
            return (moment[7] - 21 * moment[6] + 105 * moment[5] - 210 * moment[4]
                    + 252 * moment[3] - 140 * moment[2] + 24 * f_u_w) / 5040;
        }
    }

    if (K >= 8) {
        moment[8] = pow(f_u_w, 8)
                    - 28 * moment[6] * esti_var_f
                    - 140 * moment[4] * pow(esti_var_f, 2)
                    - 210 * moment[2] * pow(esti_var_f, 3)
                    - 105 * pow(esti_var_f, 4);
        if (K == 8) {
            return (moment[8] - 28 * moment[7] + 140 * moment[6] - 364 * moment[5]
                    + 560 * moment[4] - 560 * moment[3] + 336 * moment[2] - 70 * f_u_w) / 40320;
        }
    }

    if (K >= 9) {
        moment[9] = pow(f_u_w, 9)
                    - 36 * moment[7] * esti_var_f
                    - 210 * moment[5] * pow(esti_var_f, 2)
                    - 420 * moment[3] * pow(esti_var_f, 3)
                    - 315 * moment[1] * pow(esti_var_f, 4);
        if (K == 9) {
            return (moment[9] - 36 * moment[8] + 168 * moment[7] - 504 * moment[6]
                    + 1260 * moment[5] - 2520 * moment[4] + 3024 * moment[3]
                    - 2016 * moment[2] + 504 * f_u_w) / 362880;
        }
    }

    if (K == 10) {
        moment[10] = pow(f_u_w, 10);
        return (moment[10] - 45 * moment[9] + 210 * moment[8] - 630 * moment[7]
                + 1260 * moment[6] - 2520 * moment[5] + 3024 * moment[4]
                - 2520 * moment[3] + 1260 * moment[2] - 210 * f_u_w) / 3628800;
    }

    // Unsupported K
    return 0;
}

double locally_compute_f_given_q_and_x(int q, int x, BiGraph& g, BiGraph& g2) {

	// cout<<"using Eps2 = "<<Eps2 <<endl;
	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;


	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;
	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;

    // cout<<"q = "<<q<<endl;
	for(auto nb: g.neighbor[q]){

        // only consider when priority nb < priority q

        // cout<<"nb of q: "<<nb<<endl;
		if(nb == x){
			x_is_a_nb_of_q = true;
		}
		if(g2.has(nb, x)){
			Nx_cap_Nq_minus_x++;
		}else{
			Nq_minus_Nx_minus_x++;
		}
	}
    // cout<<endl;

	// if x is not a neighbor of q, then 
	// if(g2.is_bipartite){
		// in the case of bipartite graphs
    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }
	// }else{
	// 	// in the case of general graphs 
	// 	if(x_is_a_nb_of_q){
	// 		assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x-1);
	// 	}else{
	// 		assert(Nq_minus_Nx_minus_x == g.degree[q] - Nx_cap_Nq_minus_x); 
	// 	}
	// }
	// locally the degree of q is known. 
	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

	long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	res += noise; 

	return res;
}

double locally_compute_f_given_q_and_x_ad_hoc(int q, int x, BiGraph& g, BiGraph& g2) {

	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;


	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;
	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;
    // cout<<"x = "<<x <<endl;q
    // cout<<"visiting nb of q:"<< q<<endl;
	for(auto nb: g.neighbor[q]){

		if(nb == x){
			x_is_a_nb_of_q = true;
		}

        // need to check (nb, x) \in g2.
        if (!g2.has_computed(nb, x)){
            randomized_response_single_bit(nb, x, g, g2);
        }

        // if(g2.edge_vector[min(nb,x)][max(nb,x)]){
		// 	Nx_cap_Nq_minus_x++;
		// }else{
		// 	Nq_minus_Nx_minus_x++;
		// }
        unsigned int smaller = (nb < x) ? nb : x;
        unsigned int larger = (nb < x) ? x : nb;

        if (g2.edge_vector[smaller][larger]) {
            Nx_cap_Nq_minus_x++;
        } else {
            Nq_minus_Nx_minus_x++;
        }
	}


    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }

	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

	long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	res += noise; // why without noise it is even bigger? 

	return res;
}


double locally_compute_f_given_q_and_x_two_graphs(int q, int x, BiGraph& g, BiGraph& g2, BiGraph& g3) {

	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;

	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;

    double Nx_cap_Nq_minus_x_dup =0, Nq_minus_Nx_minus_x_dup =0;

	bool x_is_a_nb_of_q = false;

    // cout<<"q = "<<q<<endl;
	for(auto nb: g.neighbor[q]){

		if(nb == x){
			x_is_a_nb_of_q = true;
		}
		if(g2.has(nb, x)){
			Nx_cap_Nq_minus_x++;
		}else{
			Nq_minus_Nx_minus_x++;
		}
		if(g3.has(nb, x)){
			Nx_cap_Nq_minus_x_dup++;
		}else{
			Nq_minus_Nx_minus_x_dup++;
		}
	}

    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }

	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

    double res2 = Nx_cap_Nq_minus_x_dup * gamma__ + Nq_minus_Nx_minus_x_dup * (-p)/(1-2*p); 

    res = (res + res2)/2;

    // the GS is the same.
	res += stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	return res;
}


// count the common neighbor of q and x 
// where prority is less than q 
double locally_compute_f_given_q_and_x_vp(int q, int x, BiGraph& g, BiGraph& g2, int& res__) {

	// cout<<"using Eps2 = "<<Eps2 <<endl;
	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;


	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;
	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;

    // cout<<"q = "<<q<<endl;
	for(auto nb: g.neighbor[q]){

        // only consider when priority nb < priority q
        if(g.prio[q] <= g.prio[nb] ) 
            continue;

        assert(g.prio[q] > g.prio[nb]);

        // cout<<"nb of q: "<<nb<<endl;
		if(nb == x){
			x_is_a_nb_of_q = true;
		}
		if(g2.has(nb, x)){
			Nx_cap_Nq_minus_x++;
		}else{
			Nq_minus_Nx_minus_x++;
		}
	}
    // cout<<endl;

	// if x is not a neighbor of q, then 
	// if(g2.is_bipartite){
		// in the case of bipartite graphs
    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }

	// locally the degree of q is known. 
	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

    res__ = Nx_cap_Nq_minus_x + Nq_minus_Nx_minus_x;

	long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	res += noise; 

	return res;
}

double locally_compute_f_given_q_and_x_vp_2(int q, int x, BiGraph& g, BiGraph& g2, int& res__) {

	// cout<<"using Eps2 = "<<Eps2 <<endl;
	double res =-1;
	int start, end;

    start = ( g2.is_upper(q) ) ?   g2.num_v1 : 0 ;
    end =   ( g2.is_upper(q) ) ?   g2.num_nodes() : g2.num_v1;


	double Nx_cap_Nq_minus_x =0, Nq_minus_Nx_minus_x =0;
	// it looks like edge clipping makes it worse
	bool x_is_a_nb_of_q = false;

    // cout<<"q = "<<q<<endl;
	for(auto nb: g.neighbor[q]){

        // only consider when priority nb < priority q
        if(g.prio[x] <= g.prio[nb] ) 
            continue;

        assert(g.prio[x] > g.prio[nb]);

        // cout<<"nb of q: "<<nb<<endl;
		if(nb == x){
			x_is_a_nb_of_q = true;
		}
		if(g2.has(nb, x)){
			Nx_cap_Nq_minus_x++;
		}else{
			Nq_minus_Nx_minus_x++;
		}
	}
    // cout<<endl;

	// if x is not a neighbor of q, then 
	// if(g2.is_bipartite){
		// in the case of bipartite graphs
    if(x_is_a_nb_of_q){
        cout<<"x cannot be a neighbor of q. They are from the same layer."<<endl;
        cout<<"x = "<<x <<endl;
        cout<<"q = "<<q <<endl;
        exit(1);
    }

	// locally the degree of q is known. 
	res = Nx_cap_Nq_minus_x * gamma__ + Nq_minus_Nx_minus_x * (-p)/(1-2*p); 

    res__ = Nx_cap_Nq_minus_x + Nq_minus_Nx_minus_x;

	long double noise = stats::rlaplace(0.0, (gamma__/Eps2), engine); 

	res += noise; 

	return res;
}

// this works
long double weighted_pair_sampling_non_DP(BiGraph& g, unsigned long seed) {
    // non-DP version.
    
    int T = 10000; // we want to sample 100 vertex pairs from L(G)
    double m = 0;
    for (int i = 0; i < g.num_v1; ++i){
        m += g.degree[i];
    }
    cout<<"m = "<<m<<endl;
    assert(m==g.num_edges);

    double res = 0; 

    init_genrand(seed);
    for (int k = 0; k < T; k++) {
        // cout << "pair # " << k << endl;
        // Sample vertex u
        double rand_val_u = genrand_real2();  // Generate a single random number
        double cumulative_prob_u = 0.0;
        int u = -1;
        for (int i = 0; i < g.num_v1; ++i) {
            cumulative_prob_u += (double)g.degree[i] / m;  // Accumulate probabilities

            if (rand_val_u < cumulative_prob_u) {
                u = i;  // Select the current vertex if rand_val_u falls within its probability range
                break;
            }
        }
        double rand_val_w = genrand_real2();  // Generate a single random number
        double cumulative_prob_w = 0.0;
        int w = -1;
        for (int i = 0; i < g.num_v1; ++i) {
            cumulative_prob_w += (double)g.degree[i] / m;  // Accumulate probabilities

            if (rand_val_w < cumulative_prob_w) {
                w = i;  // Select the current vertex if rand_val_w falls within its probability range
                break;
            }
        }
        // Ensure u and w are not -1
        assert((u != -1) && (w != -1));

        if (u != w) {
            // compute the common neighbors 
            double common_u_w = 0;
            for (int xxx : g.neighbor[u]) {
                if (find(g.neighbor[w].begin(), g.neighbor[w].end(), xxx) != g.neighbor[w].end()) {
                    common_u_w++;
                }
            }
            // cout<<"common neighbor = "<< common_u_w <<endl;
            if(common_u_w>0){
                double common_u_w_choose_two = common_u_w * (common_u_w - 1) /2 ; 
                res += common_u_w_choose_two * m * m *1.0 / (2 * g.degree[u] * g.degree[w]);
            }
        }
    }
    return res / T; 
}

// long double binomial(int n, int k) {
//     if (k < 0 || n < k) return 0.0;
//     if (k == 0) return 1.0;
//     long double res = 1.0;
//     for (int i = 0; i < k; ++i) {
//         res *= (n - i);
//         res /= (i + 1);
//     }
//     return res;
// }

long double binomial(int n, int k) {
    if (n < 0 || k < 0 || k > n) return 0.0;
    if (k == 0 || k == n) return 1.0;
    k = std::min(k, n - k); // Optimize by using smaller k
    long double res = 1.0;
    for (int i = 0; i < k; ++i) {
        res *= static_cast<long double>(n - i) / (i + 1);
    }
    return res;
}



// one-round biclique counting: 
// _switch = btf:0, cate:1, biclique:2, quasi-biclique: 3.
long double one_round_biclique(BiGraph& g, unsigned long seed, 
    int p__, int q__){

	BiGraph g2(g); 
    
	construct_noisy_graph(g, g2, seed);

    cout<<"p__ = "<<p__ <<endl;
    cout<<"q__ = "<<q__ <<endl;


    std::vector<int> U, L; 
    for (int i = 0; i < g.num_v1; ++i) U.push_back(i); 
    
    for (int i = g.num_v1; i < g.num_nodes(); ++i) L.push_back(i);
    

    // Generate all combinations of p vertices from U
    generate_combinations(U, p__, up_options);
    generate_combinations(L, q__, lo_options);
	
	cout<<"counting biclique on noisy graph"<<endl;

    
    std::vector<long double> m__(p__ * q__ + 1, 0);
    // Convert adjacency lists to unordered_set for O(1) edge lookup
    std::vector<std::unordered_set<int>> adj(g2.neighbor.size());
    for (size_t i = 0; i < g2.neighbor.size(); i++) {
        
        adj[i] = std::unordered_set<int>( g2.neighbor[i].begin(), g2.neighbor[i].end());

    }


    std::cout << "Old way: Counting mi numbers\n";
    #pragma omp parallel
    {
        std::vector<long double> local_m(p__ * q__ + 1, 0);  // Private array for each thread

        #pragma omp for collapse(2) nowait
        for (size_t up_idx = 0; up_idx < up_options.size(); up_idx++) {
            for (size_t lo_idx = 0; lo_idx < lo_options.size(); lo_idx++) {
                const auto& xxx = up_options[up_idx];
                const auto& yyy = lo_options[lo_idx];

                int num_edges = 0;
                for (int u : xxx) {
                    for (int v : yyy) {
                        if (adj[u].count(v)) num_edges++; // O(1) edge lookup
                    }
                }
                local_m[num_edges]++;
            }
        }
        // Reduce results
        #pragma omp critical
        for (size_t i = 0; i <= p__ * q__; i++) { // Start from 0
            #pragma omp atomic
            m__[i] += local_m[i];
        }
    }
    

    
    /*
    vector<long double> m__(p__ * q__ + 1, 0);
	cout<<"counting mi numbers"<<endl;
	#pragma omp parallel for collapse(2)
	for(auto xxx:up_options){ 
		for(auto yyy:lo_options){

            // for each motif, check how many incident edges are there
			int num_edges = 0;
			for(auto i:xxx){
				for(auto j:yyy){
                    // they all need to be checked using g2. 
					if(std::find(g2.neighbor[i].begin(), g2.neighbor[i].end(), j) != g2.neighbor[i].end() ){
						num_edges++;
					}
				}
			}
			#pragma omp atomic
			m__[num_edges]++;
		}
	}
    */
    

	for(int i=0;i<m__.size();i++){
		cout<<"edge = "<<i<<" num = "<<m__[i]<<endl;
	}



	long double res = 0, mu = exp(Eps); 
	
	// if(_switch==2){
    // hard to obtain the other motif counts
    for(int i=0;i<m__.size();i++){
        res +=  power(-mu, i) * m__[i]; 
    }
    res /= power(1-mu, p__*q__);
    naive_estis[iteration] = m__[m__.size()-1];
	// }

    /*
	if(_switch==3){

		res += -6*mu* m__[0];

		res += (5*power(mu,2) + 1) * m__[1];

		res += -2*mu*(2*power(mu,2) + 1) * m__[2];

		res += 3*power(mu,2)*(power(mu,2) + 1) * m__[3];

		res += -2*power(mu,3)*(power(mu,2) + 2) * m__[4];

		res += power(mu,4)*(power(mu,2) + 5) * m__[5];

		res += -6*power(mu,5) * m__[6];

		res /= power(mu-1, 6);

		naive_estis[iteration] = m__[m__.size()-2];
	}
    */
	// double t2 = omp_get_wtime();

	// record time elapsed.
	// RR_time += t1-t0;
	// server_side_time += t2-t1;
	//
    cout<<"naive esti = "<<naive_estis[iteration] <<endl;
	return res;
}


long double one_round_biclique_2_3(BiGraph& g, unsigned long seed) {
    BiGraph g2(g);
    construct_noisy_graph(g, g2, seed);

    int p__ = 2;
    int q__ = 3;
    std::cout << "p__ = " << p__ << "\n";
    std::cout << "q__ = " << q__ << "\n";

    // Get upper and lower vertex indices
    std::vector<int> U, L;
    for (int i = 0; i < g.num_v1; ++i) U.push_back(i);
    for (int i = g.num_v1; i < g.num_nodes(); ++i) L.push_back(i);
    int n1 = U.size();
    int n2 = L.size();

    std::cout << "Counting biclique on noisy graph\n";

    // Convert adjacency lists to unordered_set for O(1) lookup
    std::vector<std::unordered_set<int>> adj(g2.neighbor.size());
    for (size_t i = 0; i < g2.neighbor.size(); ++i) {
        adj[i] = std::unordered_set<int>(g2.neighbor[i].begin(), g2.neighbor[i].end());
    }

    // Array to store motif counts (B_i for i = 0 to 6)
    std::vector<long double> m__(p__ * q__ + 1, 0);


    long double sum___ =  0 ; 
    /**
    std::cout << "Old way: Counting mi numbers\n";

    // Generate all combinations of p vertices from U
    cout<<"Generate all combinations of p vertices from U "<<endl;
    generate_combinations(U, p__, up_options);
    generate_combinations(L, q__, lo_options);

    #pragma omp parallel
    {
        std::vector<long double> local_m(p__ * q__ + 1, 0);  // Private array for each thread

        #pragma omp for collapse(2) nowait
        for (size_t up_idx = 0; up_idx < up_options.size(); up_idx++) {
            for (size_t lo_idx = 0; lo_idx < lo_options.size(); lo_idx++) {
                const auto& xxx = up_options[up_idx];
                const auto& yyy = lo_options[lo_idx];

                int num_edges = 0;
                for (int u : xxx) {
                    for (int v : yyy) {
                        if (adj[u].count(v)) num_edges++; // O(1) edge lookup
                    }
                }
                local_m[num_edges]++;
            }
        }
        // Reduce results
        #pragma omp critical
        for (size_t i = 0; i <= p__ * q__; i++) { // Start from 0
            #pragma omp atomic
            m__[i] += local_m[i];
        }
    }
    // Output motif counts

    sum___ =  0 ; 
    for (size_t i = 0; i < m__.size(); ++i) {
        std::cout << "# edge = " << i << " num = " << m__[i] << "\n";
        sum___ += m__[i] ; 
    }

    cout<<"sum = "<< sum___ <<endl;
    // exit(0);
    */


    // reset these numbers: 
    for (size_t i = 0; i < m__.size(); ++i) {
        m__[i] = 0;
    }

    std::cout << "New way: Counting Bi numbers (optimized)\n";
    #pragma omp parallel
    {
        std::vector<long double> local_m(p__ * q__ + 1, 0); // Private array for each thread

        #pragma omp for collapse(2) nowait
        for (size_t u1_idx = 0; u1_idx < n1; ++u1_idx) {
            for (size_t u2_idx = u1_idx + 1; u2_idx < n1; ++u2_idx) {
                int u1 = U[u1_idx];
                int u2 = U[u2_idx];

                // Compute s2 = |N(u1) ∩ N(u2)|
                int s2 = 0;
                for (int v : g2.neighbor[u1]) {
                    if (adj[u2].count(v)) ++s2;
                }

                // Compute degrees
                int deg_u1 = g2.neighbor[u1].size();
                int deg_u2 = g2.neighbor[u2].size();

                // Compute s1 = |N(u1) ∪ N(u2)| - s2 = deg(u1) + deg(u2) - 2 * s2
                int s1 = deg_u1 + deg_u2 - 2 * s2;

                // Compute s0 = n2 - |N(u1) ∪ N(u2)| = n2 - (deg(u1) + deg(u2) - s2)
                int s0 = n2 - (deg_u1 + deg_u2 - s2);

                // Compute B_i for i = 0 to 6
                local_m[0] += binomial(s0, 3);
                local_m[1] += binomial(s0, 2) * binomial(s1, 1);
                local_m[2] += binomial(s0, 1) * binomial(s1, 2) + binomial(s0, 2) * binomial(s2, 1);
                local_m[3] += binomial(s1, 3) + binomial(s0, 1) * binomial(s1, 1) * binomial(s2, 1);
                local_m[4] += binomial(s1, 2) * binomial(s2, 1) + binomial(s0, 1) * binomial(s2, 2);
                local_m[5] += binomial(s1, 1) * binomial(s2, 2);
                local_m[6] += binomial(s2, 3);
            }
        }

        // Reduce results
        #pragma omp critical
        for (size_t i = 0; i <= p__ * q__; ++i) {
            #pragma omp atomic
            m__[i] += local_m[i];
        }
    }

    // Output motif counts
    sum___ =  0 ; 
    for (size_t i = 0; i < m__.size(); ++i) {
        std::cout << "# edge = " << i << " num = " << m__[i] << "\n";
        sum___ += m__[i];
    }

    cout<<"sum = "<< sum___ <<endl;

    long double target = binomial(n1, 2) * binomial(n2, 3) ; 
    cout<<"real totoal = "<<target  <<endl;

    assert(sum___ == target);

    // Compute unbiased estimate using Theorem 2
    long double res = 0, mu = std::exp(Eps);
    for (size_t i = 0; i < m__.size(); ++i) {
        res += std::pow(-mu, i) * m__[i];
    }
    res /= std::pow(1 - mu, p__ * q__);
    naive_estis[iteration] = m__[m__.size() - 1];

    std::cout << "naive esti = " << naive_estis[iteration] << "\n";
    return res;
}

long double one_round_biclique_2_K(BiGraph& g, int K, unsigned long seed) {
    BiGraph g2(g);
    construct_noisy_graph(g, g2, seed);

    int p__ = 2;
    int q__ = K;
    std::cout << "p__ = " << p__ << "\n";
    std::cout << "q__ = " << q__ << "\n";

    // Get upper and lower vertex indices
    std::vector<int> U, L;
    for (int i = 0; i < g.num_v1; ++i) U.push_back(i);
    for (int i = g.num_v1; i < g.num_nodes(); ++i) L.push_back(i);
    int n1 = U.size();
    int n2 = L.size();

    std::cout << "Counting biclique on noisy graph\n";

    // Convert adjacency lists to unordered_set for O(1) lookup
    std::vector<std::unordered_set<int>> adj(g2.neighbor.size());
    for (size_t i = 0; i < g2.neighbor.size(); ++i) {
        adj[i] = std::unordered_set<int>(g2.neighbor[i].begin(), g2.neighbor[i].end());
    }

    // Array to store motif counts (B_i for i = 0 to 2*K)
    std::vector<long double> m__(2 * K + 1, 0.0);

    std::cout << "Counting Bi numbers (optimized)\n";
    #pragma omp parallel
    {
        std::vector<long double> local_m(2 * K + 1, 0.0); // Private array for each thread

        #pragma omp for collapse(2) nowait
        for (size_t u1_idx = 0; u1_idx < n1; ++u1_idx) {
            for (size_t u2_idx = u1_idx + 1; u2_idx < n1; ++u2_idx) {
                int u1 = U[u1_idx];
                int u2 = U[u2_idx];

                // Compute s2 = |N(u1) ∩ N(u2)|
                int s2 = 0;
                for (int v : g2.neighbor[u1]) {
                    if (adj[u2].count(v)) ++s2;
                }

                // Compute degrees
                int deg_u1 = g2.neighbor[u1].size();
                int deg_u2 = g2.neighbor[u2].size();

                // Compute s1 = deg(u1) + deg(u2) - 2 * s2
                int s1 = deg_u1 + deg_u2 - 2 * s2;

                // Compute s0 = n2 - (deg(u1) + deg(u2) - s2)
                int s0 = n2 - (deg_u1 + deg_u2 - s2);

                // Compute contributions to m__[i] for i = 0 to 2*K
                for (int c = 0; c <= K; ++c) {
                    for (int b = 0; b <= K - c; ++b) {
                        int a = K - b - c;
                        if (a >= 0) {
                            int i = 2 * c + b;
                            long double contrib = binomial(s0, a) * binomial(s1, b) * binomial(s2, c);
                            local_m[i] += contrib;
                        }
                    }
                }
            }
        }

        // Reduce results
        #pragma omp critical
        for (size_t i = 0; i <= 2 * K; ++i) {
            #pragma omp atomic
            m__[i] += local_m[i];
        }
    }

    cout<<"motif count distribution"<<endl;
    for (size_t i = 0; i <= 2 * K; ++i) {
        cout<<"i = "<<i<<" mi = "<< m__[i] <<endl;
    }

    // Output motif counts and verify sum
    long double sum___ = std::accumulate(m__.begin(), m__.end(), 0.0);
    long double target = binomial(n1, 2) * binomial(n2, K);
    std::cout << "Sum of motif counts: " << sum___ << "\n";
    std::cout << "Expected total: " << target << "\n";
    if (std::abs(sum___ - target) > 1e-6 * target) {
        // when this happens, there is the issue of overflow
        std::cout << "Sum of motif counts: " << sum___ << "\n";
        std::cout << "Expected total: " << target << "\n";
        std::cerr << "Warning: Sum of motif counts does not match expected total.\n";
    }

    // Compute unbiased estimate using Theorem 2
    long double res = 0, mu = std::exp(Eps);
    for (size_t i = 0; i <= 2 * K; ++i) {
        res += std::pow(-mu, i) * m__[i];
    }
    res /= std::pow(1 - mu, 2 * K);
    naive_estis[iteration] = m__[2 * K];

    std::cout << "naive esti = " << naive_estis[iteration] << "\n";
    return res;
}



// Naive biclique count
long double naive_biclique(BiGraph& g, unsigned long seed, 
    int p__, int q__){

	BiGraph g2(g); 
    
    long double res = 0; 

	construct_noisy_graph(g, g2, seed);

    // if(p__ == 2 && q__ ==2){
    //     res = BFC_EVP(g2);
    //     cout<<"btf res = "<<res<<endl;
    //     return res;
    // }

    // only do this when larger biclique
    biGraph convertedGraph = convertBiGraphTobiGraph(g2);

    cout << "Converted graph: n1=" << convertedGraph.n1 << ", n2=" << convertedGraph.n2 << ", m=" << convertedGraph.m << std::endl;
    
    // Create and use BCListPlusPlus
    BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, p__, q__);

    res = counter->exactCount();

    cout<<"res = "<<res<<endl;

	return res;
}
// todo: implement the vertex-priority-based wedge butterfly counting 

// TODO: (1) use estimated priority instead
// (2) use estimated value of number of priority-obeying neighbors
// (3) if we do not care about degree, just use vertex id, what will happen?
long double VP_wedge_based_two_round_btf(BiGraph& g, unsigned long seed) {

    Eps0 = Eps * 0.05;
    // private_estimate_of_degrees(g);
	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
	for(int i=0;i<g.num_nodes();i++){
		deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
	}

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;
    // Eps1 = (Eps -Eps0) * 0.6 ;
    // Eps2 = (Eps -Eps0) * 0.4 ;

    // Phase 1. RR
    double t1 = omp_get_wtime();
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);
	construct_noisy_graph(g, g2, seed);  // upload noisy edges
    // unfortunately, this step cannot be run in parallel

    // Phase 2. local counting, this step can benefit from parallism 
	// each vertex u download the noisy graph. 
    double t2 = omp_get_wtime();
    cout << "local counting" << endl;
	Eps2 = Eps - Eps1 - Eps0;
    
	// cout<<"using Eps2 = "<<Eps2 <<endl;
	gamma__ = (1-p) / (1-2*p);
    // use eps2
    // long double global_sensitivity, sum = 0;
	long double res___ = 0; 

	#pragma omp parallel
	{
	#pragma omp for schedule(static)

        // for each u in V
        // looks like this becomes much slower.
        // for imbalanced layers, we can just go ahead and choose the smaller layer 
        // for normal layers, we apply vertex-priority
		for(int u =0 ; u < g.num_nodes() ; u++) {

            // for each w with priority smaller than u
            // cout<<"u = "<<u<<endl;
            int start__, end__; 
            if(g.is_upper(u)){
                start__ = 0;
                end__ = g.num_v1;
            } else {  // Ensure this is a proper else block
                start__ = g.num_v1;
                end__ = g.num_nodes();
            }
			for(int w =start__ ; w <end__ ; w++) {

				if(u==w) 
					continue;

                // cout<<"w = "<<w<<endl;
                assert(g.same_layer(u,w));

				if(vertex_pair_reduction &&  g.prio[u] <= g.prio[w]) 
					continue; 
                
                assert(g.prio[u] > g.prio[w]);
                // we only consider each pair once
                // how do we get the common neighbors of u and we in g? 


                // computing the ground truth
                long double real_f_u_w = 0;
                // real_f_u_w = count_if(g.neighbor[u].begin(), g.neighbor[u].end(), [&](auto xx){ return g.has(xx,w); });

                // the number of priority-obeying wedges estimated from N(u)
                int x1, x2;
				long double f_u_w = locally_compute_f_given_q_and_x_vp(u, w, g, g2, x1);

                // the number of priority-obeying wedges estimated from N(w)
                long double f_w_u = locally_compute_f_given_q_and_x_vp_2(w, u, g, g2, x2);

                // averaging always gives pretty good result 
                f_u_w = (f_u_w + f_w_u)/2;

				long double local_res = f_u_w * f_u_w - f_u_w; 


                // deg_estis[u] should be the number of priority-obeying neighbor of u 
                // deg_estis[w] should be the number of priority-obeying neighbor of w
                long double esti_var_f  = 0;
                long double esti_var_f_1 = 2 * gamma__ * gamma__ / (Eps2 * Eps2); 

                // need to change x1 and x2 to be esitimates
                esti_var_f_1 += p * (1 - p) * x1 / ((1 - 2 * p) * (1 - 2 * p));

                long double esti_var_f_2 = 2 * gamma__ * gamma__ / (Eps2 * Eps2); 

                esti_var_f_2 += p * (1 - p) * x2 / ((1 - 2 * p) * (1 - 2 * p));

                esti_var_f = (esti_var_f_1 + esti_var_f_2)/4;
				
				local_res -= esti_var_f; 

				// if the degree estimate of u is more accurate then it's better
				#pragma omp critical
				res___ += local_res/2; // incrementing butterfly(u,w)

			}

		}
	}
    return res___;

}

void fetch_or_compute_biclique_count(int P___, int K___, 
    string dataset, BiGraph& g){

    sqlite3* db;
    if (sqlite3_open("../biclq_counts.db", &db) != SQLITE_OK) {
        std::cerr << "Error opening database: " << sqlite3_errmsg(db) << std::endl;
        exit(1);
    }
    // Dataset, p, and q values to filter
    size_t found = dataset.find_last_of("/\\");  // Find the last slash
    std::string dataset_to_find = dataset.substr(found + 1);  // Extract part after the last slash

    // std::string dataset_to_find = "unicode";  // Example dataset
    int p_to_find = P___;                   // Example p value
    int q_to_find = K___;                   // Example q value

    // Query to retrieve one row based on dataset, p, and q
    const char* sql = "SELECT dataset, count FROM pqbiclique_counts WHERE dataset = ? AND p = ? AND q = ?;";
    sqlite3_stmt* stmt;

    // Prepare statement
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
        std::cerr << "Error preparing statement: " << sqlite3_errmsg(db) << std::endl;
        sqlite3_close(db);
        exit(1);
    }

    // Bind parameters
    sqlite3_bind_text(stmt, 1, dataset_to_find.c_str(), -1, SQLITE_STATIC);
    sqlite3_bind_int(stmt, 2, p_to_find);
    sqlite3_bind_int(stmt, 3, q_to_find);

    // Execute the query
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        std::string dataset = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        unsigned long long count = static_cast<unsigned long long>(sqlite3_column_int64(stmt, 1));
        std::cout << "Dataset: " << dataset << ", p = " << P___ << ", q = " << K___ << ", biclique count = " << count << std::endl;
        real = count;
    } else {
        std::cout << "No matching data found." << std::endl;

        // when this happens, we need to compute count ad hoc.
        biGraph convertedGraph = convertBiGraphTobiGraph(g);
        std::cout << "Converted graph: n1=" << convertedGraph.n1 
                  << ", n2=" << convertedGraph.n2 
                  << ", m=" << convertedGraph.m << std::endl;
        BCListPlusPlus* counter = new BCListPlusPlus(&convertedGraph, P___, K___);
        real = counter->exactCount();
        cout<<"cliq count = "<<real<<endl;

        // Insert the new value into the database
        const char* insert_sql = "INSERT INTO pqbiclique_counts (dataset, p, q, count) VALUES (?, ?, ?, ?);";
        sqlite3_stmt* insert_stmt;

        // Prepare the INSERT statement
        if (sqlite3_prepare_v2(db, insert_sql, -1, &insert_stmt, nullptr) != SQLITE_OK) {
            std::cerr << "Error preparing insert statement: " << sqlite3_errmsg(db) << std::endl;
            sqlite3_close(db);
            exit(1);
        }

        // Bind parameters for the insert statement
        sqlite3_bind_text(insert_stmt, 1, dataset_to_find.c_str(), -1, SQLITE_STATIC);
        sqlite3_bind_int(insert_stmt, 2, P___);
        sqlite3_bind_int(insert_stmt, 3, K___);
        sqlite3_bind_int64(insert_stmt, 4, real);

        // Execute the insert statement
        if (sqlite3_step(insert_stmt) != SQLITE_DONE) {
            std::cerr << "Error inserting data: " << sqlite3_errmsg(db) << std::endl;
        } else {
            std::cout << "Inserted new biclique count into database." << std::endl;
        }
        // Cleanup the insert statement
        sqlite3_finalize(insert_stmt);

    }
    // Cleanup
    sqlite3_finalize(stmt);
    sqlite3_close(db);
}

long double weighted_pair_sampling(BiGraph& g, unsigned long seed) {

    // question: should we select left side or right side? 
    init_genrand(seed);

    // Phase 0. deg_esti_time records the maximum degree perturbation time.
    Eps0 = Eps*0.1 ;

    double m = 0;

    bool sample_from_upper = (g.num_v1 > (g.num_nodes() - g.num_v1));

    // sample_from_upper = false ; 

    // we choose the layer with more vertices, seems better.
    // what if we take avg from both options 

	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());

	vector<long double> deg_weights; 
	deg_weights.resize(g.num_nodes());

    // issue: when degree is negative, how to sample vertices? 
    // can we just skip the vertices with negative degrees? 
    // (1) skipp the vertices with negative degrees
    // okay, we can set a threshold. 
    // if the degree is less than 0, set it to be 1.
    int start = sample_from_upper ? 0 : g.num_v1;
    int end = sample_from_upper ? g.num_v1 : g.num_nodes();

    for (int i = start; i < end; ++i) {
        deg_estis[i] = g.degree[i] + stats::rlaplace(0.0, 1 / Eps0, engine);

        deg_weights[i] =deg_estis[i];
        // if(deg_weights[i]<=0) deg_weights[i] = 0; 

        // why we cannot even use degrees anymore? 
        deg_weights[i] = 1.0 * pow(g.degree[i],2);
        // deg_weights[i] = 1.0*g.degree[i] ;

        // we want to use estimates from vertices with higher vertices 
        m += deg_weights[i];
    }

    // obtain T vertex pairs with replacemnt, sampling 10K pairs 
    int T = 10000;

    // Maps to store occurrences
    map<int, map<int, int>> pair_count_u;  // Map from vertex u to map of vertex w and count

    // what if I always use both upper and lower vertex pairs? 
    // sample T pairs of vertices 
    for (int k = 0; k < T; k++) {
        int u = -1, w = -1;
        // Sample u from the appropriate set
        double rand_val_u = genrand_real2();  
        // cout<<"rand_val_u = "<<rand_val_u <<endl;
        double cumulative_prob_u = 0.0;  
        int start = sample_from_upper ? 0 : g.num_v1;
        int end = sample_from_upper ? g.num_v1 : g.num_nodes();
        for (int i = start; i < end; ++i) {
            cumulative_prob_u += (double) deg_weights[i] / m;
            if (rand_val_u < cumulative_prob_u) {
                u = i;
                break;
            }
        }
        // cout<<"u  = "<<u<<endl;
        double rand_val_w = genrand_real2();
        // cout<<"rand_val_w = "<<rand_val_w <<endl;
        double cumulative_prob_w = 0.0;
        for (int i = start; i < end; ++i) {
            cumulative_prob_w += (double) deg_weights[i] / m;
            if (rand_val_w < cumulative_prob_w) {
                w = i;
                break;
            }
        }
        // cout<<"w  = "<<w<<endl;
        assert(u != -1 && w != -1);
        if (u != w) {
            pair_count_u[u][w]++;
        }
        // looks like we do not need to store them and re-process, 
        // we could just do computations here 
    }

    // // Phase 1. RR
    // this can be made faster by computing g2(u,v) adhoc.
    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;

    p = 1.0 / (exp(Eps1) + 1.0);
    BiGraph g2(g);

    g2.edge_vector.clear();
    assert(g2.edge_vector.size()==0);
    g2.edge_vector.resize(g.num_v1 + g.num_v2);
    
    if(!efficient_RR){
        cout << "constructing noisy_graph(g); " << endl;
        construct_noisy_graph(g, g2, seed);  // upload noisy edges
    }
    // we do not construct the whole graph here 


    // Output results
    // cout << "Pairs for each u:" << endl;
    double res = 0;
    gamma__ = (1-p) / (1-2*p);

    // counting phase (apply RR in an ad-hoc manner.)
    cout<<"// counting phase (apply RR in an ad-hoc manner.) "<<endl;
    for (const auto& entry_u : pair_count_u) {
        int u = entry_u.first;
        // cout << "Vertex u = " << u << ":" << endl;

        for (const auto& pair : entry_u.second) {
            int w = pair.first;

            // cout << "\tVertex w = " << w << ":" << endl;

            int count = pair.second;

            assert(u!=w);

            // why do we need to assert the estimated degrees are positive?

            long double f_u_w = -1; 
            if(efficient_RR){
                // cout<<"computing from NB of u"<<endl;
                f_u_w = locally_compute_f_given_q_and_x_ad_hoc(u, w, g, g2); 
                // cout<<"computing from NB of w"<<endl;
                f_u_w += locally_compute_f_given_q_and_x_ad_hoc(w, u, g, g2); 
            }else{
                f_u_w = locally_compute_f_given_q_and_x(u, w, g, g2); 
                f_u_w = locally_compute_f_given_q_and_x(w, u, g, g2); 
            }

            f_u_w /= 2; 

            long double local_res = f_u_w * f_u_w - f_u_w; 

            // compute the variance of tilde(f)
            long double esti_var_f = 2 * gamma__ * gamma__ / (Eps2 * Eps2); 
            esti_var_f += p * (1 - p) * deg_estis[u] / ((1 - 2 * p) * (1 - 2 * p));
            if(averaging_f_estimates){
                long double esti_var_f_2 = 2 * gamma__ * gamma__ / (Eps2 * Eps2); 
                esti_var_f_2 += p * (1 - p) * deg_estis[w] / ((1 - 2 * p) * (1 - 2 * p));	
                esti_var_f = (esti_var_f + esti_var_f_2)/4;
            }
            local_res -= esti_var_f; 
            // local_res is an estimator for (common choose 2)
            // cout<<"common esti = "<<local_res<<endl;
            // res += count * local_res * m * m / (2 * deg_estis[u] * deg_estis[w]) ;
            assert(deg_weights[u]!=0);
            assert(deg_weights[w]!=0);
            res += count * local_res * m * m / (2 * deg_weights[u] * deg_weights[w]) ;
        }
    }

    res /=2;

    return res/T;
}
// this function builds two noisy grahs to improve accuracy
long double wedge_based_btf_avg(BiGraph& g, unsigned long seed) {

    Eps0 = Eps * 0.05;

	vector<long double> deg_estis; 
	deg_estis.resize(g.num_nodes());
	for(int i=0;i<g.num_nodes();i++){
		deg_estis[i] = g.degree[i]+stats::rlaplace(0.0, 1/(Eps0), engine); 
        // there is a chance that deg_estis is negative 
	}

    Eps1 = Eps * 0.6;
    Eps2 = Eps - Eps1 - Eps0;


    // Phase 1. RR
    double t1 = omp_get_wtime();
    cout << "construct_noisy_graph(g); " << endl;
    
    p = 1.0 / (exp(Eps1) + 1.0);

    // build one noisy graph by applying RR on U(G)
    BiGraph g2(g);
	construct_noisy_graph(g, g2, seed);  // upload noisy edges

    BiGraph g3(g);
	construct_noisy_graph_2(g, g3, seed);  // upload noisy edges


    // unfortunately, this step cannot be run in parallel

    // Phase 2. local counting, this step can benefit from parallism 
	// each vertex u download the noisy graph. 
    double t2 = omp_get_wtime();
    cout << "local counting" << endl;
	Eps2 = Eps - Eps1 - Eps0;
    
	// cout<<"using Eps2 = "<<Eps2 <<endl;
	gamma__ = (1-p) / (1-2*p);
    // use eps2
    // long double global_sensitivity, sum = 0;
	long double res___ = 0; 

	// what if we only consider upper vertices ?  --> better efficiency and effect  
	// #pragma omp parallel for reduction(+:res___)


	
	int start__, end__; 
	start__ = g.num_v1 < g.num_v2 ? 0 : g.num_v1; 
	end__ = g.num_v1 < g.num_v2 ? g.num_v1 : g.num_nodes(); 
	#pragma omp parallel
	{
	#pragma omp for schedule(static)
		for(int u =start__ ; u <end__ ; u++) {

			for(int w =start__ ; w <end__ ; w++) {
				if(u==w) 
					continue;

				if(vertex_pair_reduction && u<w) 
					continue; 
                
                // we only consider each pair once
                // how do we get the common neighbors of u and we in g? 


                // computing the ground truth
                long double real_f_u_w = 0;
                // real_f_u_w = count_if(g.neighbor[u].begin(), g.neighbor[u].end(), [&](auto xx){ return g.has(xx,w); });

                // phi is used in here. 

                // computing the local wedge estimators using two noisy graphs: 
                long double f_u_w = locally_compute_f_given_q_and_x_two_graphs(u, w, g, g2, g3);
                long double f_w_u = locally_compute_f_given_q_and_x_two_graphs(w, u, g, g2, g3);

				// long double f_u_w = locally_compute_f_given_q_and_x(u, w, g, g2);
                // long double f_w_u = locally_compute_f_given_q_and_x(w, u, g, g2);

                long double diff1 =0, diff2 = 0;

                // averaging always gives pretty good result 
				if(averaging_f_estimates){// using the average of fu and fw.
					f_u_w = (f_u_w + f_w_u)/2;
				}
				long double local_res = f_u_w * f_u_w - f_u_w; 

                
                long double esti_var_f  = 0;
                long double esti_var_f_1 = 2 * gamma__ * gamma__ / (Eps2 * Eps2); 
                esti_var_f_1 += p * (1 - p) * deg_estis[u] / (2 * (1 - 2 * p) * (1 - 2 * p));

                long double esti_var_f_2 = 2 * gamma__ * gamma__ / (Eps2 * Eps2); 
                esti_var_f_2 += p * (1 - p) * deg_estis[w] / (2 * (1 - 2 * p) * (1 - 2 * p));

				if(averaging_f_estimates){
					esti_var_f = (esti_var_f_1 + esti_var_f_2)/4;
				}

				local_res -= esti_var_f; 

				// if the degree estimate of u is more accurate then it's better
				#pragma omp critical
				res___ += local_res; // incrementing butterfly(u,w)

			}

		}
	}
	if (vertex_pair_reduction){
        // we are running this version
		return res___/2; 
	}else{
		return res___/4; 
        // this is because we have computed the butterfly count for u, w 
        // and w, u
	}

}


// biclique related code
bool BCListPlusPlus::costEstimateRaw() {
    std::default_random_engine e(10007);
    std::uniform_int_distribution<unsigned> chooseU(0, g->n1 - 1);
    std::uniform_int_distribution<unsigned> chooseV(0, g->n2 - 1);
    std::vector<uint32_t> sum(std::max(g->n1, g->n2));
    std::vector<uint32_t> tmp(std::max(g->n1, g->n2));
    int l = 0;
    uint32_t rd = 10;
    uint32_t Du = 0;
    uint32_t u = chooseU(e);
    for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
        uint32_t v = g->e1[i];
        for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
            uint32_t w = g->e2[j];

            if(w != u) {
                if(sum[w] == 0) tmp[l++] = w;
                sum[w]++;
            }
        }
    }
    for(int i = 0; i < l; i++) {
        if(sum[tmp[i]] >= (uint32_t)q) {
            Du++;
        }
        sum[tmp[i]] = 0;
    }
    l = 0;
    
    // t = rd;
    uint32_t Dv = 0;
    // while(t--) {
    uint32_t v = chooseV(e);
    // for(uint32_t v = 0; v < g->n2; v += rd) {
    for(uint32_t i = g->pV[v]; i < g->pV[v + 1]; i++) {
        uint32_t u = g->e2[i];

        for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
            uint32_t w = g->e1[j];

            if(w != v) {
                if(sum[w] == 0) tmp[l++] = w;
                sum[w]++;
            }
        }
    }
    

    for(int i = 0; i < l; i++) {
        if(sum[tmp[i]] >= (uint32_t)p) {
            Dv++;
        }
        sum[tmp[i]] = 0;
    }
    l = 0;
    
    

    double costU = pow(1.0 * Du / (g->n1 / 10), p - 2) * Du * g->maxDu;
    double costV = pow(1.0 * Dv / (g->n2 / 10), q - 2) * Dv * g->maxDv;
    // printf("cost:%.2f %.2f, du %u, dv %u\n", costU, costV, Du, Dv);

    return costU > costV;
}

bool BCListPlusPlus::costEstimate() {

    std::vector<uint32_t> sum(std::max(g->n1, g->n2));
    std::vector<uint32_t> tmp(std::max(g->n1, g->n2));
    int l = 0;

    uint32_t rd = 100;
    // uint32_t t = rd;
    uint32_t Du = 0, maxDu, rdu = 0;
    double sumDu = 0.0;
    // while(t--) {
        // uint32_t u = chooseU(e);
    for(uint32_t u = 0; u < g->n1; u += rd) {
        rdu++;
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];

                if(w > u) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }

        for(int i = 0; i < l; i++) {
            if(sum[tmp[i]] >= (uint32_t)q) {
                Du++;
            }
            sum[tmp[i]] = 0;
        }
        l = 0;

        maxDu = std::max(maxDu, Du);
        sumDu += Du;
    }
    
    uint32_t Dv = 0, maxDv, rdv = 0;
    double sumDv = 0.0;
    // while(t--) {
        // uint32_t v = chooseV(e);
    for(uint32_t v = 0; v < g->n2; v += rd) {
        rdv++;
        for(uint32_t i = g->pV[v]; i < g->pV[v + 1]; i++) {
            uint32_t u = g->e2[i];
    // if(u >= g->n1) {
    //     printf("n2 %u, %u, %u, i %u\n", g->n2, v, u, i);fflush(stdout);
    // }
            for(uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
                uint32_t w = g->e1[j];

                if(w > v) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }

        for(int i = 0; i < l; i++) {
            if(sum[tmp[i]] >= (uint32_t)p) {
                Dv++;
            }
            sum[tmp[i]] = 0;
        }
        l = 0;

        maxDv = std::max(maxDv, Dv);
        sumDv += Dv;
    }
    
    sumDu = sumDu / rdu * g->n1;
    sumDv = sumDv / rdv * g->n2;

    double avgDu = std::max(2.0, sumDu / g->n1);
    double avgDv = std::max(2.0, sumDv / g->n2);

    double costU = pow(avgDu, p - 2) * sumDu;
    double costV = pow(avgDv, q - 2) * sumDv;
    // printf("cost:%.2f %.2f, du %u, dv %u\n", costU, costV, Du, Dv);

    return costU > costV;
}

long double BCListPlusPlus::exactCount() {

    collect2HopNeighbors();
    // printf("collect 2\n"); fflush(stdout);

    S.resize(g->n2);

    uint32_t maxDegree = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        maxDegree = std::max(maxDegree, (uint32_t)H.lists[u].size());
    }


    tmpNodes.resize(p + 1);
    tmpNodes[0].resize(g->n1);
    for(int i = 1; i <= p; i++) {
        tmpNodes[i].resize(maxDegree);
    }
    H.nodes.resize(g->n1);

    H.d.resize(p + 1);
    for(int i = 0; i <= p; i++) {
        H.d[i].resize(g->n1);
    }
    for(uint32_t u = 0; u < g->n1; u++) {
        H.d[0][u] = H.lists[u].size();
    }

    ans = 0;

    layerBasedListing(0, g->n1, g->n2);

    printf("ans %.2Lf\n", ans);

    // unsigned long long int res = ans;

    return ans ; 
    // fflush(stdout);
}

void BCListPlusPlus::layerBasedListing(int l, int pH, int pS) {


    if(l == p) {
        ans += C[pS][q];
        return;
    }
    H.nodes.copy(tmpNodes[l].data(), pH);
// for(int i = 0; i < pH; i++) {
//     printf("%u ", tmpNodes[l][i]);
// }printf("\n");fflush(stdout);

auto t1 = std::chrono::steady_clock::now();
    for(int j = 0; j < pH; j++) {

        uint32_t u = tmpNodes[l][j];
        // uint32_t u = H.nodes[l][j];
        if(H.lists[u].size() < uint32_t(p - l - 1)) {
            continue;
        }

        int pSS = 0;
        if(l == 0) {
            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                S.changeTo(g->e1[i], pSS++);
            }
        }
        else {
            for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                if(S.idx(g->e1[i]) < (uint32_t)pS) {
                    S.changeTo(g->e1[i], pSS++);
                }
            }
        }

        if(pSS < q) continue;
        
        int pHH = 0;

        for(int i = 0; i < H.d[l][u]; i++) {
            auto v = H.lists[u][i];
            if(H.nodes.idx(v) < (uint32_t)pH) {
                H.nodes.changeTo(v, pHH++);
            }
            // if(H.nodes[l].idx(v) < (uint32_t)pH) {
            //     H.nodes[l + 1].changeTo(v, pHH++);
            // }
        }

        if(l + 1 < p)
        for(int i = 0; i < pHH; i++) {
            // uint32_t u = H.nodes[l + 1][i];
            uint32_t u = H.nodes[i];
            int d = H.d[l][u];
            for(int k = 0; k < d; k++) {
                auto v = H.lists[u][k];
                if(H.nodes.idx(v) >= pHH) {
                    std::swap(H.lists[u][k], H.lists[u][--d]);
                    --k;
                }
                // if(H.nodes[l + 1].idx(v) >= pHH) {
                //     std::swap(H.lists[u][k], H.lists[u][--d]);
                //     --k;
                // }
            }
            H.d[l + 1][u] = d;
        }

        layerBasedListing(l + 1, pHH, pSS);
    }
}

void BCListPlusPlus::collect2HopNeighbors() {
    H.lists.resize(g->n1);

    std::vector<uint32_t> sum(g->n1);
    std::vector<uint32_t> tmp(g->n1);
    uint32_t l = 0;

double twotwo = 0;
    for(uint32_t u = 0; u < g->n1; u++) {
        for(uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for(uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                if(w > u) {
                    if(sum[w] == 0) tmp[l++] = w;
                    sum[w]++;
                }
            }
        }

        for(uint32_t i = 0; i < l; i++) {
            uint32_t w = tmp[i];
            if(sum[w] >= q) {
twotwo += C[sum[w]][q];
                H.lists[u].push_back(w);
            }
            sum[w] = 0;
        }
        l = 0;
    }
// printf("2-2 clique %.0f\n", twotwo);
}

