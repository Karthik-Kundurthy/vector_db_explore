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


#define UNUSED(p)  ((void)(p))



// Represent vector field(s) and unstructured metadata
struct DataPoint {
    std::string id;
    std::vector<float> vector_field;

    DataPoint(std::string &id, std::vector<float> vec) : id(id), vector_field(vec) {}

    DataPoint(const std::string &id, size_t dim) : id(id), vector_field(dim, 0.0f) {}

};

struct HNSWNode {
    DataPoint data_point;
    int max_level;
    std::vector<std::vector<HNSWNode*>> connections; 
    

    HNSWNode(const DataPoint& data_point, int max_level) :  data_point(data_point), max_level(max_level), connections(max_level + 1) {}

    bool layer_in_bounds(int layer) {
        return layer >= 0 && layer <= max_level;
    }

    void addConnection(HNSWNode* neighbor, uint64_t layer) {
        if (layer_in_bounds(layer)) {
            connections[layer].push_back(neighbor);
        }

    }

};

/// Manage the entire dataset: SIFT and GLOVE
class DataSet {
    public:
        std::string name;
        std::vector<DataPoint> dataset;

        DataSet(std::string name) : name(name) {}

        void addDataPoint(const DataPoint& data_point) {
            dataset.push_back(data_point);
        }

        // void removeDataPoint(const DataPoint&) {
            
        // }
};


class HNSW {
    public:

        int max_level = 0; // initially at 0 bc no nodes in the graph yet
        HNSWNode* entry_point = nullptr;
        float m_l; // normalization factor
    
        HNSW(float m_l) : m_l(m_l) {}

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


            int l = getRandomLevel();
            int L = max_level;
            HNSWNode* ep = new HNSWNode(q,l);

            // entry point has not yet been initialized
            if (entry_point == nullptr) {
                entry_point = ep;
                max_level = l; 
            } else {

                for (int l_c = L; l_c >= l + 1; l_c--) {
                    std::vector<HNSWNode*> W = search_layer(q,std::vector<HNSW*>(),1,l_c);
                    ep = getNearest(W,q);
                }
                
                for (int l_c = std::min(L,l); l_c >= 0; l_c --) {
                    std::vector<HNSWNode*> W = search_layer(q,std::vector<HNSW*>(),ef_construction,l_c);
                    std::vector<HNSWNode*> neighbors = select_neighbors_simple(q,W,m); // using the simple algorithm for now

                    // add bidirectionall connectionts from neighbors to q at layer lc
                    for (auto* neighbor: neighbors) {
                        ep -> addConnection(neighbor,l);
                        neighbor -> addConnection(ep,l);
                    }

                }

                if (l > max_level ) {
                    max_level = l;
                    entry_point = ep;
                }

            }



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

            // default value
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