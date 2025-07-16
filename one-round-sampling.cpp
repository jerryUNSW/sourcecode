// this is not useful
long double one_round_biclique_sampling(BiGraph& g, unsigned long seed, int p__, int q__){


    assert(p__ == 2 && q__ == 3); // We support only (2,3) motifs here

    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist_u(0, g.num_v1 - 1);

    std::vector<std::unordered_set<int>> adj(g.neighbor.size());
    for (int i = 0; i < g.num_nodes(); ++i)
        adj[i] = std::unordered_set<int>(g.neighbor[i].begin(), g.neighbor[i].end());

    std::vector<long double> m__(7, 0); // 0 to 6 edges possible

    int num_samples  = 10000;

    #pragma omp parallel
    {
        std::mt19937 local_rng(seed ^ omp_get_thread_num());
        std::vector<long double> local_m(7, 0);

        #pragma omp for
        for (int sample = 0; sample < num_samples; ++sample) {
            // Sample 2 distinct upper nodes
            int u1 = dist_u(local_rng), u2 = dist_u(local_rng);
            while (u2 == u1) u2 = dist_u(local_rng);

            // Compute intersection of neighbors
            std::vector<int> common_neighbors;
            for (int v : g.neighbor[u1]) {
                if (adj[u2].count(v)) common_neighbors.push_back(v);
            }

            if (common_neighbors.size() < 3) continue; // skip if too few common neighbors

            // Sample 3 distinct lower nodes from common neighbors
            std::shuffle(common_neighbors.begin(), common_neighbors.end(), local_rng);
            int v1 = common_neighbors[0];
            int v2 = common_neighbors[1];
            int v3 = common_neighbors[2];

            // Count number of edges among the 6 possible u-v pairs
            int count = 0;
            for (int u : {u1, u2}) {
                for (int v : {v1, v2, v3}) {
                    if (adj[u].count(v)) count++;
                }
            }

            local_m[count]++;
        }

        // Reduce
        #pragma omp critical
        {
            for (int i = 0; i <= 6; ++i)
                m__[i] += local_m[i];
        }
    }

    // Print estimated histogram
    for (int i = 0; i < m__.size(); ++i)
        std::cout << "edge = " << i << " est count = " << m__[i] << std::endl;

    // Normalize if needed
    long double total = std::accumulate(m__.begin(), m__.end(), 0.0L);
    if (total > 0) {
        for (auto& val : m__) val /= total;
    }

    // Optional: return anything derived from histogram
    naive_estis[iteration] = m__[6]; // or any estimator you want
    // return m__[6];



	long double res = 0, mu = exp(Eps); 
	
    for(int i=0;i<m__.size();i++){
        res +=  power(-mu, i) * m__[i]; 
    }
    res /= power(1-mu, p__*q__);
    naive_estis[iteration] = m__[m__.size()-1];

	// double t2 = omp_get_wtime();

	// record time elapsed.
	// RR_time += t1-t0;
	// server_side_time += t2-t1;
	//
    cout<<"naive esti = "<<naive_estis[iteration] <<endl;
	return res;

    /*
    long double sampling_ratio  = 0.1 ; 

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
	


    // Sampling 0.1% of the combinations
    std::shuffle(up_options.begin(), up_options.end(), std::default_random_engine(seed));
    std::shuffle(lo_options.begin(), lo_options.end(), std::default_random_engine(seed + 1)); // Different seed for diversity
    size_t up_sample_size = std::max<size_t>(1, up_options.size() * sampling_ratio);
    size_t lo_sample_size = std::max<size_t>(1, lo_options.size() * sampling_ratio);
    up_options.resize(up_sample_size);
    lo_options.resize(lo_sample_size);



	cout<<"counting biclique on noisy graph"<<endl;


    std::vector<long double> m__(p__ * q__ + 1, 0);
    // Convert adjacency lists to unordered_set for O(1) edge lookup
    std::vector<std::unordered_set<int>> adj(g2.neighbor.size());
    for (size_t i = 0; i < g2.neighbor.size(); i++) {
        adj[i] = std::unordered_set<int>( g2.neighbor[i].begin(), g2.neighbor[i].end());
    }


    std::cout << "Counting mi numbers (sampling-based)\n";
    #pragma omp parallel
    {
        std::vector<long double> local_m(p__ * q__ + 1, 0);  // Private array for each thread

        #pragma omp for collapse(2) nowait
        for (size_t up_idx = 0; up_idx < up_options.size(); up_idx++) {
            for (size_t lo_idx = 0; lo_idx < lo_options.size(); lo_idx++) {
                const auto& xxx = up_options[up_idx];
                const auto& yyy = lo_options[lo_idx];
                int num_edges = 0;

                // for each commbination.
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
    
    // list the results
	for(int i=0;i<m__.size();i++){
        m__[i] = m__[i] / pow(0.001, 2);
		cout<<"edge = "<<i<<" num = "<<m__[i]<<endl;
	}

	long double res = 0, mu = exp(Eps); 
	
    for(int i=0;i<m__.size();i++){
        res +=  power(-mu, i) * m__[i]; 
    }
    res /= power(1-mu, p__*q__);
    naive_estis[iteration] = m__[m__.size()-1];

    cout<<"naive esti = "<<naive_estis[iteration] <<endl;
	return res;
    */
}

