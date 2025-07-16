#pragma once
#ifndef __BIGRAPH_H
#define __BIGRAPH_H
#include "utility.h"

using namespace std;
using namespace std::chrono;

class BiGraph {
   public:
    BiGraph(std::string dir);
    BiGraph(const BiGraph& other);
    BiGraph();
    ~BiGraph() {
    }
    // note that addEdgeRaw is processing raw ids.
    void addEdgeRaw(vid_t u, vid_t v);

    // addEdge is processing real ids.
    void addEdge(vid_t u, vid_t v);
    // deleteEdge is processing real ids.
    void deleteEdge(vid_t u, vid_t v);

   public:
    void computePriority();
    void init(unsigned int num_v1, unsigned int num_v2);
    void loadGraph(std::string dir);
    void print_graph();
    void show_nodes();
    // vector<LinearHeap> make_heaps();
    bool same_layer(vid_t u, vid_t v) {
        if (is_upper(u) && is_upper(v)) return true;
        if (is_lower(u) && is_lower(v)) return true;
        return false;
    }

    bool diff_layer(vid_t u, vid_t v) { return !same_layer(u, v); }

    bool is_upper(vid_t u) const { return (u < this->num_v1); }
    bool is_lower(vid_t u) const { return (u >= this->num_v1); }
    // void generalGraph(string dir);

    bool com_p(vid_t u, vid_t v) {
        // return prio[u] > prio[v];
        return prio[u] < prio[v];
    }

    bool is_active(int vertex) { return degree[vertex] > 0; }
    int num_nodes() { return this->num_v1 + this->num_v2; }

    // checks the neighbors 
    // check whether i and j are connected 
    bool has(int i, int j) {
        return std::find(neighbor[i].begin(), neighbor[i].end(), j) !=
               neighbor[i].end();
    }

    // for each pair of vertices u, v, it checks whether this pair has been computed. 
	bool has_computed(int u, int v){
		assert(u!=v);
		int smaller = (u < v) ? u : v;
		int larger = (u < v) ? v : u;
		return edge_vector[smaller].find(larger) != edge_vector[smaller].end();
	}



    bool hasEdge(int u, int v) {
        if (u < 0 || u >= num_nodes() || v < 0 || v >= num_nodes()) {
            cout << "Invalid vertex index\n";
            // return false;
            exit(1);
        }
        // Check if u and v belong to different layers
        if (same_layer(u, v)) {
            // Invalid edge between upper and lower vertices
            cout << "Invalid edge between upper and lower vertices" << endl;
            // return false;
            exit(1);
        }
        // Get the degrees of vertices u and v
        int degreeU = degree[u];
        int degreeV = degree[v];

        // Determine which vertex to scan first based on their degrees
        int vertexToScanFirst = (degreeU < degreeV) ? u : v;
        int vertexToScanSecond = (vertexToScanFirst == u) ? v : u;

        // Check if vertexToScanFirst is connected to vertexToScanSecond
        for (const auto& neighbor : neighbor[vertexToScanFirst]) {
            if (neighbor == vertexToScanSecond) {
                // Edge found
                return true;
            }
        }
        // Edge not found
        return false;
    }

    std::string dir;
    num_t num_v1;
    num_t num_v2;
    long num_edges;
    vector<vector<vid_t>> neighbor;
    vector<int> degree;
    vector<int> prio;
    vector<int> U, V;
    vector<int> nodes;
   public:
    // max and min degrees
    int v1_max_degree;
    int v2_max_degree;

    vector<map<int, bool>> edge_vector;
};



struct biGraph {
    uint32_t n1, n2, m, maxDu, maxDv;
    uint32_t n[2];
    bool connect(uint32_t u, uint32_t v, int i) {
        if(i == 0) return connectUV(u, v);
        else return connectUV(v, u);
    }

    struct Edge{
        uint32_t u, v;
    };
    std::vector<Edge> edges;
    
    std::vector<uint32_t> pU, e1, pV, e2;

    biGraph() {}

    // reading from my datafiles, this is modified by me.
    // I am not using this function in this project tho. 
    /*
    biGraph(const std::string &filePath) {
        std::string metaFile = filePath + ".meta";
        std::string edgeFile = filePath + ".e";
    
        fastIO meta(metaFile, "r");
        fastIO in(edgeFile, "r");

        // Read n1, n2, m from meta file
        n1 = meta.getUInt();
        n2 = meta.getUInt();
        m = meta.getUInt();

        printf("graph size: n1 = %u, n2 = %u, m = %u\n", n1, n2, m);
        fflush(stdout);

        // Allocate space
        edges.resize(m);
        e1.resize(m);
        e2.resize(m);
        pU.resize(n1 + 5);
        pV.resize(n2 + 5);

        // printf("done allocating space \n");
        
        for (uint32_t i = 0; i < m; i++) {
            uint32_t u = in.getUInt();
            uint32_t v = in.getUInt();

            // Print input values before processing
            // printf("Reading edge %u: (u=%u, v=%u)\n", i, u, v);
            fflush(stdout); // Ensure output is printed immediately

            edges[i].u = u;          // Upper node (U) remains the same
            edges[i].v = v;          // Shift lower node (V) by n1
            // maybe this is unnecessary
        }

        // printf("done reading edges space \n");

        // Reorder graph for efficiency
        changeToDegreeOrder();

        // printf("// Reorder graph for efficiency \n");
    }
    */

    void coreReductionFast22() {
        uint32_t n = n1 + n2;
        ListLinearHeap heap(n, std::max(maxDu, maxDv) + 1);
        
        std::vector<uint32_t> ids(n);
        std::vector<uint32_t> keys(n);
        std::vector<uint32_t> labelsL(n1);
        std::vector<uint32_t> labelsR(n2);
        std::vector<bool> visL(n1 + 1), visR(n2 + 1);

        for(uint32_t i = 0; i < n1; i++) {
            ids[i] = i;
            keys[i] = deg1(i) + 1;
        }
        for(uint32_t i = n1; i < n; i++) {
            ids[i] = i;
            keys[i] = deg2(i - n1) + 1;
        }

        heap.init(n, std::max(maxDu, maxDv) + 1, ids.data(), keys.data());
        
        for(uint32_t i = 0; i < n; i++) {
            uint32_t u, degU;

            if(!heap.pop_min(u, degU)) printf("errorLheap\n");
            if(degU >= 2 + 1) {
                if(u < n1) visL[u] = true;
                else visR[u - n1] = true;
                
                continue; 
            }

            if(u < n1) {
                // labelsL[u] = lId++;
                for(uint32_t j = pU[u]; j < pU[u + 1]; j++) {
                    heap.decrement(e1[j]+n1);
                }
            }
            else {
                // labelsR[u - n1] = rId++;
                for(uint32_t j = pV[u - n1]; j < pV[u - n1 + 1]; j++) {
                    heap.decrement(e2[j]);
                }
            }
        }

        uint32_t pL = 1, pR = 1;
        for(uint32_t u = 0; u < n1; u++) {
            if(visL[u]) labelsL[u] = pL++;
        }
        for(uint32_t v = 0; v < n2; v++) {
            if(visR[v]) labelsR[v] = pR++;
        }
        
        uint32_t pm = 0;
        for(uint32_t u = 0; u < n1; u++) {
            if(!visL[u]) continue;
            for(uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                uint32_t v = e1[i];
                if(visR[v]) {
                    edges[pm].u = labelsL[u] - 1;
                    edges[pm].v = labelsR[v] - 1;
                    ++pm;
                } 
            }
        }
        m = pm;

        n1 = pL - 1;
        n2 = pR - 1;
        // printf("n1 %u, n2 %u, m %u\n", n1, n2, m);
        
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        // std::fill(d1.begin(), d1.begin() + n1 + 1, 0);
        // std::fill(d2.begin(), d2.begin() + n2 + 1, 0);
        changeToDegreeOrder();

        

        printf("core reduction\n");
    }

    void coreReduction(int p, int q) {
        std::queue<uint32_t> qL, qR;
        std::vector<uint32_t> d1(n1 + 1), d2(n2 + 1);
        std::vector<uint32_t> labelsL(n1 + 1), labelsR(n2 + 1);
        std::vector<bool> visL(n1 + 1), visR(n2 + 1);
                
        for(uint32_t i = 0; i < n1; i++) {
            d1[i] = deg1(i);
            if(deg1(i) < q) {
                qL.push(i);
                visL[i] = true;
            }
        }
        for(uint32_t i = 0; i < n2; i++) {
            d2[i] = deg2(i);
            if(deg2(i) < p) {
                qR.push(i);
                visR[i] = true;
            }
        }

        while(!qL.empty() || !qR.empty()) {
            while(!qL.empty()) {
                uint32_t u = qL.front(); qL.pop();

                for(uint32_t i = 0; i < d1[u]; i++) {
                    uint32_t v = e1[pU[u] + i];
                    // if(d2[v] < q) continue;

                    for(uint32_t j = pV[v]; j < pV[v] + d2[v]; j++) {
                        if(e2[j] == u) {
                            --d2[v];
                            std::swap(e2[j], e2[pV[v] + d2[v]]);

                            if(d2[v] == p - 1 && !visR[v]) {
                                qR.push(v);
                                visR[v] = true;
                            }
                            break;
                        }
                    }
                }
            }

            while(!qR.empty()) {
                uint32_t v = qR.front(); qR.pop();

                for(uint32_t i = 0; i < d2[v]; i++) {
                    uint32_t u = e2[pV[v] + i];
                    // if(d1[u] < p) continue;

                    for(uint32_t j = pU[u]; j < pU[u] + d1[u]; j++) {
                        if(e1[j] == v) {
                            --d1[u];
                            std::swap(e1[j], e1[pU[u] + d1[u]]);

                            if(d1[u] == q - 1 && !visL[u]) {
                                qL.push(u);
                                visL[u] = true;
                            }
                            break;
                        }
                    }
                }
            }
        }

        uint32_t pL = 1, pR = 1;
        for(uint32_t u = 0; u < n1; u++) {
            if(!visL[u]) labelsL[u] = pL++;
        }
        for(uint32_t v = 0; v < n2; v++) {
            if(!visR[v]) labelsR[v] = pR++;
        }
        
        uint32_t pm = 0;
        for(uint32_t u = 0; u < n1; u++) {
            if(visL[u]) continue;
            for(uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                uint32_t v = e1[i];
                if(!visR[v]) {
                    edges[pm].u = labelsL[u] - 1;
                    edges[pm].v = labelsR[v] - 1;
                    ++pm;
                } 
            }
        }
        m = pm;

        n1 = pL - 1;
        n2 = pR - 1;
        // printf("n1 %u, n2 %u, m %u\n", n1, n2, m);
        
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        std::fill(d1.begin(), d1.begin() + n1 + 1, 0);
        std::fill(d2.begin(), d2.begin() + n2 + 1, 0);
        changeToDegreeOrder();

        // std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        // std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        // std::fill(d1.begin(), d1.begin() + n1 + 1, 0);
        // std::fill(d2.begin(), d2.begin() + n2 + 1, 0);

        // for(uint32_t i = 0; i < m; i++) {
        //     ++d1[edges[i].u];
        //     ++d2[edges[i].v];
        // }

        // maxDu = 0;
        // for(uint32_t i = 0; i < n1; i++) {
        //     maxDu = std::max(maxDu, d1[i]);
        // }

        // maxDv = 0;
        // for(uint32_t i = 0; i < n2; i++) {
        //     maxDv = std::max(maxDv, d2[i]);
        // }

        // for(uint32_t i = 0; i < n1; i++) {
        //     pU[i + 1] = pU[i] + d1[i];
        // }
        // for(uint32_t i = 0; i < n2; i++) {
        //     pV[i + 1] = pV[i] + d2[i];
        // }

        // for(uint32_t i = 0; i < m; i++) {
        //     e1[ pU[edges[i].u]++ ] = edges[i].v;
        // }
        // for(uint32_t i = 0; i < m; i++) {
        //     e2[ pV[edges[i].v]++ ] = edges[i].u;
        // }

        // pU[0] = pV[0] = 0;
        // for(uint32_t i = 0; i < n1; i++) {
        //     pU[i + 1] = pU[i] + d1[i];
        // }
        // for(uint32_t i = 0; i < n2; i++) {
        //     pV[i + 1] = pV[i] + d2[i];
        // }

        // for(uint32_t i = 0; i < n1; i++) {
        //     std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        // }
        // for(uint32_t i = 0; i < n2; i++) {
        //     std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        // }

        // printf("core reduction\n");
    }

    void changeToCoreOrder() {
        std::vector<uint32_t> d1, d2;
        
        d1.resize(n1);
        d2.resize(n2);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        maxDu = 0;
        for(uint32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }
        maxDv = 0;
        for(uint32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }

        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }
        for(uint32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v; 
        }
        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        } 
        
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u; 
        }
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }

        ListLinearHeap lheap(n1, maxDu + 1), rheap(n2, maxDv + 1);
        uint32_t n = std::max(n1, n2);
        std::vector<uint32_t> ids(n);
        std::vector<uint32_t> keys(n);
        std::vector<uint32_t> labelsL(n1);
        std::vector<uint32_t> labelsR(n2);

        for(uint32_t i = 0; i < n1; i++) {
            ids[i] = i;
            keys[i] = d1[i] + 1;
        }
        lheap.init(n1, maxDu + 1, ids.data(), keys.data());
        for(uint32_t i = 0; i < n2; i++) {
            ids[i] = i;
            keys[i] = d2[i] + 1;
        }
        rheap.init(n2, maxDv + 1, ids.data(), keys.data());
        
        uint32_t minN = std::min(n1, n2);
        for(uint32_t i = 0; i < minN; i++) {
            uint32_t u, degU;
            uint32_t v, degV;

            if(!lheap.pop_min(u, degU)) printf("errorLheap\n");
            if(!rheap.pop_min(v, degV)) printf("errorRheap\n");

            labelsL[u] = i;
            labelsR[v] = i;
            for(uint32_t j = pU[u]; j < pU[u + 1]; j++) {
                rheap.decrement(e1[j]);
            }
            for(uint32_t j = pV[v]; j < pV[v + 1]; j++) {
                lheap.decrement(e2[j]);
            }
        }
        if(n1 < n2) {
            for(uint32_t j = n1; j < n2; j++) {
                uint32_t v, degV;
                if(!rheap.pop_min(v, degV)) printf("errorRheap\n");
                labelsR[v] = j;
            }
        }
        else if(n1 > n2) {
            for(uint32_t j = n2; j < n1; j++) {
                uint32_t u, degU;
                if(!lheap.pop_min(u, degU)) printf("errorRheap\n");
                labelsL[u] = j;
            }
        }

        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = labelsL[edges[i].u];
            edges[i].v = labelsR[edges[i].v];
        }

        std::fill(d1.begin(), d1.begin() + n1, 0);
        std::fill(d2.begin(), d2.begin() + n2, 0);
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < m; i++) {
            e1[ pU[edges[i].u]++ ] = edges[i].v;
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[ pV[edges[i].v]++ ] = edges[i].u;
        }

        pU[0] = pV[0] = 0;
        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for(uint32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }

        // print();
        // fflush(stdout);
    }

    void rawOrder() {
        std::vector<uint32_t> d1, d2;

        d1.resize(n1);
        d2.resize(n2);

        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }

        maxDu = 0;
        for(uint32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }
        maxDv = 0;
        for(uint32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }

        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        }
        for(uint32_t i = 0; i < m; i++) {
            e1[pU[edges[i].u]++] = edges[i].v; 
        }
        pU[0] = 0;
        for(uint32_t u = 0; u < n1; u++) {
            pU[u + 1] = d1[u] + pU[u];
        } 
        
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }   
        for(uint32_t i = 0; i < m; i++) {
            e2[pV[edges[i].v]++] = edges[i].u; 
        }
        pV[0] = 0;
        for(uint32_t v = 0; v < n2; v++) {
            pV[v + 1] = d2[v] + pV[v];
        }

        for(uint32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for(uint32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }
    }

    void changeToDegreeOrder() {
        std::vector<uint32_t> d1, d2;
        std::vector<uint32_t> label1, label2;

        d1.resize(std::max(n1, n2) + 1);
        d2.resize(std::max(n1, n2) + 1);
        label1.resize(n1);
        label2.resize(n2);
        
        // printf("line 547 \n");
        for(uint32_t i = 0; i < m; i++) {
            ++d1[edges[i].u];
            ++d2[edges[i].v];
        }
        
        // printf("line 553 \n");
        maxDu = 0;
        for(uint32_t i = 0; i < n1; i++) {
            maxDu = std::max(maxDu, d1[i]);
        }

        maxDv = 0;
        for(uint32_t i = 0; i < n2; i++) {
            maxDv = std::max(maxDv, d2[i]);
        }
        // printf("line 563 \n");
        for(uint32_t i = 0; i < n1; i++) {

            pV[d1[i] + 1]++;
        // printf("%u-%u-%u\n", i, d1[i] + 1, pU[d1[i] + 1]);
        }

        for(uint32_t i = 0; i < maxDu; i++) {
            pV[i + 1] += pV[i];
        }
        // printf("line 578 \n");
        // printf("labels:");
        for(uint32_t i = 0; i < n1; i++) {
            label1[i] = pV[d1[i]]++;
        // printf("%u-%u ", i,label1[i]);
        }
        // printf("\n");
        // printf("there\n");fflush(stdout);
        for(uint32_t i = 0; i < n2; i++) {
            pU[d2[i] + 1]++;
        }
        for(uint32_t i = 0; i < maxDv; i++) {
            pU[i + 1] += pU[i];
        }
        // printf("line 592 \n");
        for(uint32_t i = 0; i < n2; i++) {
            label2[i] = pU[d2[i]]++;
        }
        // printf("line 596 \n");
        for(uint32_t i = 0; i < m; i++) {
            edges[i].u = label1[edges[i].u];
            edges[i].v = label2[edges[i].v];
        }
        // printf("line 601 \n");
        std::fill(d1.begin(), d1.begin() + std::max(n1, n2) + 1, 0);
        std::fill(d2.begin(), d2.begin() + std::max(n1, n2) + 1, 0);
        std::fill(pU.begin(), pU.begin() + n1 + 1, 0);
        std::fill(pV.begin(), pV.begin() + n2 + 1, 0);
        // memset(buffer, 0, sizeof(uint32_t) * bufferSize);
        // printf("line 607 \n");
        for(uint32_t i = 0; i < m; i++) {
            // printf("i  %u:%u-%u \n", i, edges[i].u, edges[i].v);
            // fflush(stdout);
            ++d1[edges[i].u];
            ++d2[edges[i].v];

            // i626:190-10545 // why edgesi.v is 10545?
        }
        // printf("line 611 \n");
        // printf("\n");
        // printf("there2\n");fflush(stdout);
        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        // printf("line 619 \n");
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }
        // printf("line 623 \n");
        for(uint32_t i = 0; i < m; i++) {
            e1[ pU[edges[i].u]++ ] = edges[i].v;
        }
        for(uint32_t i = 0; i < m; i++) {
            e2[ pV[edges[i].v]++ ] = edges[i].u;
        }
        // printf("there3\n");fflush(stdout);
        pU[0] = pV[0] = 0;
        for(uint32_t i = 0; i < n1; i++) {
            pU[i + 1] = pU[i] + d1[i];
        }
        for(uint32_t i = 0; i < n2; i++) {
            pV[i + 1] = pV[i] + d2[i];
        }

        for(uint32_t i = 0; i < n1; i++) {
            std::sort(e1.begin() + pU[i], e1.begin() + pU[i + 1]);
        }
        for(uint32_t i = 0; i < n2; i++) {
            std::sort(e2.begin() + pV[i], e2.begin() + pV[i + 1]);
        }
        // printf("thereed\n");
        // print();
        // delete [] buffer;
        // printf("line 645 \n");
    }

    void swapUV() {
        std::swap(n1, n2);
        std::swap(maxDu, maxDv);
        
        pU.swap(pV);
        e1.swap(e2);
    }

    uint32_t deg1(uint32_t u) {
        return pU[u + 1] - pU[u];
    }
    uint32_t deg2(uint32_t u) {
        return pV[u + 1] - pV[u];
    }

    bool connectUV(uint32_t u, uint32_t v) {
        return std::binary_search(
            e1.begin() + pU[u], e1.begin() + pU[u + 1],
            v
        );
    }
    bool connectVU(uint32_t v, uint32_t u) {
        return std::binary_search(
            e2.begin() + pV[v], e2.begin() + pV[v + 1],
            u
        );
    }

    void print() {
        printf("U:\n");
        for(uint32_t u = 0; u < n1; u++) {
            printf("%u:", u);
            for(uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                printf("%u ", e1[i]);
            }
            printf("\n");
        }
        printf("V:\n");
        for(uint32_t v = 0; v < n2; v++) {
            printf("%u:", v);
            for(uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                printf("%u ", e2[i]);
            }
            printf("\n");
        }
    }


    //SIMD
    std::vector<hopstotchHash> hashTablesL, hashTablesR;

    void createHashTables() {
        createHashTableL();
        createHashTableR();
    }

    void createHashTableL() {
        hashTablesL.resize(n1);
        
        for(uint32_t u = 0; u < n1; u++) {
            if(deg1(u) > H) {//H is the bucket size of hopstotch hash
                hashTablesL[u].build(e1.data() + pU[u], deg1(u));
            }
        }
    }

    void createHashTableR() {
        hashTablesR.resize(n2);
        
        for(uint32_t v = 0; v < n2; v++) {
            if(deg2(v) > H) {
                hashTablesR[v].build(e2.data() + pV[v], deg2(v));
            }
        }
    }

    bool connectUVFast(uint32_t u, uint32_t v) {
        // return connectUV(u, v);
        if(hashTablesL[u].n > 0) return hashTablesL[u].contain(v);
        else {
            return connectUV(u, v);

            // if(deg1(u) < 8) return connectUV(u, v);
            // return ccSIMD(e1.data(), pU[u], v) || std::binary_search(
            //         e1.begin() + pU[u] + 8, e1.begin() + pU[v + 1], v);
        }
    }

    bool connectVUFast(uint32_t v, uint32_t u) {
        // return connectVU(v, u);
        if(hashTablesR[v].n > 0) return hashTablesR[v].contain(u);
        else {
            return connectVU(v, u);
        }
    }

};

#endif /* __BIGRAPH_H */