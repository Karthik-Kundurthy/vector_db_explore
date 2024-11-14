#include <iostream>
#include <fstream>
#include <chrono>
#include <map>
#include <vector>
#include <list>
#include <unordered_map>
#include <string>
#include <memory>
#include <sstream>
#include <thread>
#include <queue>
#include <optional>
#include <mutex>
#include <shared_mutex>
#include <cassert>
#include <cstring> 
#include <exception>
#include <atomic>
#include <set>
#include <unordered_set>
#include <cmath>
#include <random>
#include <limits>
#include <algorithm>
#include <cassert>
#include <cmath>
#include "H5Cpp.h"



#define UNUSED(p)  ((void)(p))



/*
 * Scratch Pad: Implementation level details
 * Just implemnting the paper pseudocode directly
 * Pull dataset loading modules from brute_force_search.cpp
 * HNSWIndex contains the entry node and max level as well as insert, search and neighbor search functionality
 * DataPoint object representing a vector in euclidean space
 * HNSWNode: encapsulates the level that it is currently, pointer to the DataPoint, lateral connections and downward connection
 * Multiple level HNSW have pointer to same DataPoint object
 * Should I have a single unique pointer to the DataPoint at the top level and keep moving ownership as I move downwards and reset ownership when getting bottom?
 * Future design: serialization and deserialization can be managed by allocating page for every vector and overflow pages if needed and a metadata pg
 * Reusing BufferManger from BuzzDB
*/




// Represent vector field(s) and unstructured metadata
struct DataPoint {
    public:
        int id;
        std::vector<float> vector;

        DataPoint(int id, std::vector<float> vector) : id(id), vector(vector) {}
};

struct HNSWNode {
    public:

        DataPoint* data_point; // caution, this can mean multiple references to the same data point object. 
        uint64_t max_connections;
        uint64_t curr_level;

        std::vector<HNSWNode*> neighbors;
        HNSWNode* lower_node;

        HNSWNode(DataPoint* point, uint64_t max_conn, uint64_t level) : data_point(point), max_connections(max_conn), curr_level(level){
            lower_node = nullptr; //revist

        }

        bool is_full() {
            return neighbors.size() >= max_connections;
        }


};


class HNSW {
    public:

        int max_level = 0; // initially at 0 bc no nodes in the graph yet
        HNSWNode* ep = nullptr;
        float m_l; // normalization factor
    
        HNSW(float m_l) : m_l(m_l) {}

        float euclidean_distance(const std::vector<float>& x,const std::vector<float>& y) {
            UNUSED(x);
            UNUSED(y);

            float dist = 0.0f;

            assert((x.size() == y.size()) && "Error: x and y must be equal size");

            for (size_t i = 0; i < x.size(); i++) {
                dist += (x[i] - y[i]) * (x[i] - y[i]);
            }

            return std::sqrt(dist);
        }


        // floor(-1 * ln(unif(0...1)) * m_l)
        uint64_t getLevel() {
            return std::floor(-1 * std::log(std::uniform_real_distribution<>(0.0, 1.0)) * m_l);
        }

        /// ALGORITHM 1
        /// @param q: the query vector
        /// @param m: the number of established connections
        /// @param m_max: the maximum number of connectiosn for each element per layer
        /// @param ef_construction: size of dynamic candidate list
        /// update the HNSW graph
        /// 
        void insert(DataPoint q, uint64_t m, uint64_t m_max, size_t ef_construction) {
            UNUSED(q);
            UNUSED(m);
            UNUSED(m_max);
            UNUSED(ef_construction);


            


        }

        /// ALGORITHM 2
        /// @param q: the query vector
        /// @param ep: enter points ep
        /// @param ef: number of nearest to q elements to return
        /// @param l_c: layer number
        /// @return ef closest neighbors to q
        std::vector<HNSWNode*> search_layer(DataPoint q, std::vector<HNSW*> ep, uint64_t ef, size_t l_c) {
            UNUSED(q);
            UNUSED(ep);
            UNUSED(ef);
            UNUSED(l_c);

            
            return std::vector<HNSWNode*>();
        }

        /// ALGORITHM 3
        /// @param q: base element
        /// @param C:  candidate elements
        /// @param M: number of neighbors to return
        std::vector<HNSWNode*> select_neighbors_simple(DataPoint q, std::vector<HNSWNode*> C, uint64_t M) {
            UNUSED(q);
            UNUSED(C);
            UNUSED(M);

            return std::vector<HNSWNode*>();

            
        }

        /// ALGORITHM 4
        /// @param q query vector
        /// @param C candidate elements
        /// @param M: number of neighbors to return
        /// @param extendCandidates whether or not to extend candidate list
        /// @param keepPrunedConnections whether or not to add discarded elements
        /// @param 
        std::vector<HNSWNode*> select_neighbors_heuristic(DataPoint q, std::vector<HNSWNode> C, uint64_t M, size_t l_c, bool extendCandidates=false, bool keepPrunedConnections=false) {

            UNUSED(q);
            UNUSED(C);
            UNUSED(M);
            UNUSED(l_c);
            UNUSED(extendCandidates);
            UNUSED(keepPrunedConnections);


            return std::vector<HNSWNode*>();

        }

        /// ALGORITHM 5
        /// @param q: the query vector
        /// @param k: the number of nearest neighbors to return
        /// @param ef: size of dynamic candidate list
        /// @return K nearest elements to query vector in the dynamic candidate list
        ///
        std::vector<HNSWNode*> knn_search(DataPoint q, uint64_t k, uint64_t ef) {
            UNUSED(q);
            UNUSED(k);
            UNUSED(ef);

            return std::vector<HNSWNode*>();
        }   


        int getRandomLevel() {
            return std::floor(-std::log((static_cast<float>(rand())) / RAND_MAX) / std::log(m_l));
        }

        HNSWNode* getNearest(std::vector<HNSWNode*> W, DataPoint q) {
            UNUSED(W);
            UNUSED(q);
            return W[0];

        }

        
};


int main(int argc, char* argv[]) {
    UNUSED(argc);
    UNUSED(argv);

    return 0;

}