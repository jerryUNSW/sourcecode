#include "bigraph.h"

using namespace std;

extern bool one_round;

// a copy constructor:
BiGraph::BiGraph(const BiGraph& other) {
    this->dir = other.dir;
    this->init(other.num_v1, other.num_v2);
    // this->num_edges = 0;
    this->v1_max_degree = 0;
    this->v2_max_degree = 0;
}

BiGraph::BiGraph(string dir) {
    num_v1 = 0;
    num_v2 = 0;
    num_edges = 0;
    v1_max_degree = 0;
    v2_max_degree = 0;
    neighbor.clear();
    degree.clear();
    this->dir = dir;
    // check dir name to enter big graph mode automatically.
    loadGraph(dir);
}

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        return std::hash<T1>()(p.first) ^ std::hash<T2>()(p.second);
    }
};

// default constructor
BiGraph::BiGraph() {
    dir = "";
    num_v1 = 0;
    num_v2 = 0;
    num_edges = 0;
    v1_max_degree = 0;
    v2_max_degree = 0;
    neighbor.clear();
    degree.clear();
}
void BiGraph::print_graph() {
    print_dash(50);
    cout << "print_graph() called\n";
    for (int i = 0; i < degree.size(); i++) {
        vector<vid_t> NB = neighbor[i];
        cout << i << ": ";
        vector_show<vid_t>(NB);
    }
}
void BiGraph::show_nodes() {
    cout << "show_nodes() called\n";
    // cout<<"upper nodes: ";
    for (int i = 0; i < num_v1 + num_v2; i++) {
        if (degree[i] > 0) {
            cout << i << " ";
        }
    }
    cout << endl;
}
// initialize the BiGraph with size requirements
void BiGraph::init(unsigned int num1, unsigned int num2) {
    num_v1 = num1;
    num_v2 = num2;
    num_edges = 0;
    neighbor.resize(num_v1 + num_v2);
    degree.resize(num_v1 + num_v2);

    prio.resize(num_v1 + num_v2);

    fill_n(degree.begin(), num_v1 + num_v2, 0);
    // neighborHash.resize(num_v1+num_v2);

    // edge_vector.resize(num_v1 + num_v2);
}

void BiGraph::computePriority() {
    std::vector<std::pair<int, int>> vertexDegrees(num_nodes());
    for (int i = 0; i < num_nodes(); i++)
        vertexDegrees[i] = std::make_pair(i, degree[i]);
    // Sort the vertexDegrees based on degrees and IDs
    std::sort(vertexDegrees.begin(), vertexDegrees.end(),
              [](const auto& a, const auto& b) {
                  // the vertex with higher priority has lower degree.
                  if (a.second != b.second) {
                      return a.second > b.second;
                  } else {
                      return a.first > b.first;
                  }
              });
    for (int i = 0; i < num_nodes(); i++) {
        neighbor[i].shrink_to_fit();
        sort(neighbor[i].begin(), neighbor[i].end());
        // cout<<"rank = "<<i<<",  id = "<<vertexDegrees[i].first<<",  deg =
        // "<<vertexDegrees[i].second<<endl;
        prio[vertexDegrees[i].first] = i;
    }
}

void BiGraph::loadGraph(string dir) {

    unsigned int n1, n2;
    unsigned int edges = 0;
    int u, v, r;
    string metaFile = dir + ".meta";
    string edgeFile = dir + ".e";
    FILE* metaGraph = fopen(metaFile.c_str(), "r");
    FILE* edgeGraph = fopen(edgeFile.c_str(), "r");

    // bool include_amat = true ;
    // scan the meta file and read number of nodes
    if (fscanf(metaGraph, "%d\n%d", &n1, &n2) != 2) {
        fprintf(stderr, "Bad file format: n1 n2 incorrect\n");
        exit(1);
    }
    fprintf(stdout, "n1: %d, n2: %d\n", n1, n2);
    init(n1, n2);

    // if(include_amat){a_mat = new map<int, int>[n1+n2+1];}

    while ((r = fscanf(edgeGraph, "%d %d", &u, &v)) != EOF) {
        // fprintf(stderr, "%d, %d\n", u, v);
        if (r != 2) {
            fprintf(stderr, "Bad file format: u v incorrect\n");
            exit(1);
        }
        // addEdgeRaw(u, v);

        // cout<<" u = "<<u<<endl;
        // cout<<" v = "<<u<<endl;
        // assert(!this->same_layer(u,v));

		// neighbor[u].push_back(v + num_v1);
		// neighbor[v].push_back(u);	
        neighbor[u].push_back(v + num_v1);
        neighbor[v + num_v1].push_back(u);

        num_edges++;

        // edge vector based:
        // looks like we do not need to make use of edge_vector for g. just g2. 
		// edge_vector[u][v] = true;
		// edge_vector[v][u] = true;

    }
    cout << "|E| = " << num_edges << endl;

	// computing the maximum degree.
    cout<<"computing degrees"<<endl;
	for(int u=0;u<num_nodes();u++){
		degree[u] = neighbor[u].size();
        if(is_upper(u)){
            v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u]; 
        }else{
            v2_max_degree = v2_max_degree > degree[u] ? v2_max_degree : degree[u]; 
        }
		
	}

    fclose(metaGraph);
    fclose(edgeGraph);

    cout<<"computing priority"<<endl;
    computePriority();
}

bool satisfy_bound(long double upper, long double lower, int u, int x, int v,
                   int w, BiGraph& g) {
    // return (g.degree[u] <= upper) & (g.degree[x] <= upper) & (g.degree[v] <=
    // lower) & (g.degree[w] <= lower);

    return (g.degree[u] >= upper) & (g.degree[x] >= upper) &
           (g.degree[v] >= lower) & (g.degree[w] >= lower);
    // all heavy vertices.
}

// u,v are raw ids, convert them into real ids, and then update neighbor and
// degree
void BiGraph::addEdgeRaw(vid_t u, vid_t v) {
    neighbor[u].push_back(v + num_v1);

    degree[u]++;

    neighbor[v + num_v1].push_back(u);

    degree[v + num_v1]++;

    num_edges++;
    // v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u];
    // v2_max_degree = v2_max_degree > degree[v + num_v1] ? v2_max_degree : degree[v + num_v1];
}

// maybe this function takes too long?
void BiGraph::addEdge(vid_t u, vid_t v) {
    // it looks like only these two lines are necessary
    neighbor[u].push_back(v);
    neighbor[v].push_back(u);

    degree[u]++;
    degree[v]++;
    num_edges++;
    v1_max_degree = v1_max_degree > degree[u] ? v1_max_degree : degree[u];
    v2_max_degree = v2_max_degree > degree[v] ? v2_max_degree : degree[v];
}

// u1, u2 are real ids
void BiGraph::deleteEdge(vid_t u1, vid_t u2) {
    vid_t upper, lower;
    if (same_layer(u1, u2)) {
        fprintf(stderr, "Bad input for delete edge: u v on the same layer\n");
        exit(1);
    }
    for (int i = 0; i < degree[u1]; ++i) {
        int vv = neighbor[u1][i];
        if (vv == u2) {
            swap(neighbor[u1][i], neighbor[u1][degree[u1] - 1]);
            --degree[u1];
            neighbor[u1].pop_back();
            num_edges--;
            break;
        }
    }
    for (int i = 0; i < degree[u2]; ++i) {
        int uu = neighbor[u2][i];
        if (uu == u1) {
            swap(neighbor[u2][i], neighbor[u2][degree[u2] - 1]);
            --degree[u2];
            neighbor[u2].pop_back();
            break;
        }
    }
}