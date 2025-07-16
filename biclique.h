#pragma once
#ifndef BICLIQUE
#define BICLIQUE
#include <sstream>

#include "bigraph.h"

void printMemoryUsage() ;

long double binomial(int n, int k) ;

long double one_round_biclique_2_3(BiGraph& g, unsigned long seed);

long double one_round_biclique_2_K(BiGraph& g, int K, unsigned long seed);

// convert my graph to biclique-counting expected graph:
biGraph convertBiGraphTobiGraph(BiGraph& oldGraph); 

long double weighted_pair_sampling_non_DP(BiGraph& g, unsigned long seed); 
    
double locally_compute_f_given_q_and_x_ad_hoc(int q, int x, BiGraph& g, BiGraph& g2);

void randomized_response_single_bit(int u, int v, BiGraph& g, BiGraph& g2); 

// weighted sampling based approach:
long double weighted_pair_sampling(BiGraph& g, unsigned long seed); 

// (2, K)
long double wedge_based_two_round_2_K_biclique(BiGraph& g, unsigned long seed) ; 

// layer based:
long double layer_based_wedge_based_two_round_2_K_biclique(BiGraph& g, unsigned long seed) ;

// (3, K)
long double wedge_based_two_round_3_K_biclique(BiGraph& g, unsigned long seed) ;

// general shape: 
long double wedge_based_two_round_general_biclique(BiGraph& g, 
    unsigned long seed, int P___, int K___ );

long double compute_local_res(int K, long double f_u_w, long double esti_var_f) ;

void fetch_or_compute_biclique_count(int P___, int K___, string dataset, BiGraph& g);

// improved new approach
long double wedge_based_btf_avg(BiGraph& g, unsigned long seed);

// VP based 
long double VP_wedge_based_two_round_btf(BiGraph& g, unsigned long seed);

double locally_compute_f_given_q_and_x_vp(int q, int x, BiGraph& g, BiGraph& g2, int& res__);

double locally_compute_f_given_q_and_x_vp_2(int q, int x, BiGraph& g, BiGraph& g2, int& res__);

double locally_compute_f_given_q_and_x(int q, int x, BiGraph& g, BiGraph& g2); 

double locally_compute_f_given_q_and_x_two_graphs(int q, int x, BiGraph& g, BiGraph& g2, BiGraph& g3);

void private_estimate_of_degrees(BiGraph& g);

long double BFC_EVP(BiGraph& g);

bool satisfy_bound(long double upper, long double lower, int u, int x, int v,
                   int w, BiGraph& g);

double my_genrand_real2();

void my_init_genrand(unsigned long seed);

// long double BFC_EVP_sample(BiGraph& g, long double sampling_prob);

void compute_m3_m2(long double& m4, long double& m3, long double& m2,
                   long double& m1, long double& m0, BiGraph& g2);


long double naive_biclique(BiGraph& g, unsigned long seed, int p__, int q__);

long double one_round_btf(BiGraph& g, unsigned long seed);

long double two_round_btf(BiGraph& g, unsigned long seed);

// naive noisy motif counts
void get_noisy_naive(BiGraph& g, BiGraph& g2, long double& local_btfs,
                     long double& local_cate, long double& res);

// process butterflies in batch
void BFC_EVP_noisy(BiGraph& g, BiGraph& g2, long double& BTF, long double& cate,
                   long double& res);

void construct_noisy_graph(BiGraph& g, BiGraph& g2, unsigned long seed);

void construct_noisy_graph_2(BiGraph& g, BiGraph& g2, unsigned long seed); 

long double get_wedges(BiGraph& g);

long double get_cate(BiGraph& g);

int findPercentile(std::vector<int>& data);

void test_random_number(unsigned long seed);

long double get_laplace(long double parameter);

void compute_m3_m2_2(long double& m4, long double& m3, long double& m2,
                     long double& m1, long double& m0, BiGraph& g2);


// one-round biclique counting
long double one_round_biclique(BiGraph& g, unsigned long seed, 
    int p__, int q__); 


// Custom hash function for std::vector<int>
struct VectorHasher {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t seed = v.size();
        for (const auto& i : v) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};



// biclique counting related code:
class BCListPlusPlus {
public:
    const int maxD = 100000;
    int p, q;
    // great, this specifies input P and Q.

    biGraph * g;

    double ** C, *bf3;
    void computeC() {
        int maxPQ = std::max(p, q) + 1;
        C = new double*[maxD];
        bf3 = new double[maxPQ * maxD];
        for(int i = 0; i < maxD; i++) {
            C[i] = bf3 + i * maxPQ;
        }
        C[0][0] = 1;
        C[1][0] = 1;
        C[1][1] = 1;
        for(int i = 2; i < maxD; i++) {
            C[i][0] = 1;
            if(i < maxPQ) C[i][i] = 1;
            for(int j = 1; j < i && j < maxPQ; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
    }
    long double ans;

private:
    struct twoHopGraph {
        std::vector<std::vector<uint32_t>> lists;
        std::vector<std::vector<int>> d;
        // std::vector<LinearSet> nodes;
        LinearSet nodes;
        std::vector<bool> contained;
        std::vector<std::vector<uint32_t>> containNodes;

        void print() {
            for(uint32_t i = 0; i < lists.size(); i++) {
                printf("%u:", i);
                for(auto v:lists[i]) printf("%u ", v);
                printf("\n");
            }
        }

        bool contain(uint32_t u, uint32_t v) {
            return std::binary_search(containNodes[u].begin(), containNodes[u].end(), v);
        }
    } H;
    void collect2HopNeighbors();
    void collect2HopNeighborsContained();
    void collect2HopNeighborsContainedCDAG();

    LinearSet S;
    std::vector< std::vector<uint32_t> > tmpNodes;


public:
    ~BCListPlusPlus() {
        delete [] bf3;
        delete [] C;
    }


    // this constructor function will be called 
    BCListPlusPlus(biGraph* g_, int p_, int q_) {

        auto t1 = std::chrono::steady_clock::now();
        p = p_;
        q = q_;
        computeC();

        // g = new biGraph(filePath);
        g = g_;
        g->coreReduction(p, q);


        auto t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        // std::cout << "core reduction time:" << duration.count() << "ms" << std::endl;
        // std::cout << "n1 n2 m:" << g->n1 << ' '<< g->n2 << ' ' << g->m<< std::endl;

        t2 = t1;

        // hey, maybe we dont bother doing this. 
        if(costEstimate()) {
            g->swapUV(); 
            std::swap(p, q);
            printf("swap\n");
        }

        t2 = std::chrono::steady_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        printf("costEstimate %lld ms\n", (long long)duration.count());
        printf("load graph\n");
        fflush(stdout);
    }
    
    long double exactCount();
    void layerBasedListing(int l, int pH, int pS);

    bool costEstimate();
    bool costEstimateRaw();
};

#endif