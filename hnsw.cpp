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
#include <climits>
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

        int L = 0; // initially at 0 bc no nodes in the graph yet
        HNSWNode* entry_point = nullptr;
        float m_l; // normalization factor
        uint64_t MAX_CONN = 5;
    
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
            std::mt19937 generator(std::random_device{}());
            return std::floor(-1 * std::log((std::uniform_real_distribution<>(0.0, 1.0))(generator) * m_l));
        }

        /// ALGORITHM 1
        /// @param q: the query vector
        /// @param m: the number of established connections
        /// @param m_max: the maximum number of connectiosn for each element per layer
        /// @param ef_construction: size of dynamic candidate list
        /// update the HNSW graph
        /// 
        void insert(DataPoint* q, uint64_t m, uint64_t m_max, size_t ef_construction) {
            UNUSED(q);
            UNUSED(m);
            UNUSED(m_max);
            UNUSED(ef_construction);

            // List for the currently found nearest elements
            std::vector<HNSWNode*> W;

            int l = getLevel();
            std::cout << "INSERTION LEVEL: " << l << std::endl;

            // initialize case
            if (entry_point == nullptr && L == 0) {
                std::cout << "INITIALIZATION CASE" << std::endl;
                entry_point = new HNSWNode(q, MAX_CONN, l);
                L = l;
            }


            HNSWNode* ep = entry_point;
            

            // when initialized, this should only execute once because L == l
            for (int l_c = L; l_c >= l; l_c--) {
                W = search_layer(q, ep, 1, l_c); // W has closest neighbor now
                ep = W[0]; // ? right
            }

            // when initialized, this should only execute once because L == l
            HNSWNode* parent = nullptr;
            
            UNUSED(parent);
            bool extend_level_ep = false;
            for (int l_c = std::min(L,l); l_c >= 0; l_c--) {
                HNSWNode* new_node = new HNSWNode(q, m_max, l_c);

                // ONE TIME THING, meant 
                if (l > L && !extend_level_ep) {
                    entry_point = new_node;
                    extend_level_ep = true;
                }
                W = search_layer(q, ep, ef_construction, l_c);
                std::vector<HNSWNode*> neighbors = select_neighbors_simple(new_node, W, m, l_c);
                
                UNUSED(new_node);
                for (HNSWNode* node : neighbors) {
                    UNUSED(node);
                    ;
                    node -> neighbors.push_back(new_node);
                    new_node -> neighbors.push_back(node);

                    if (node -> neighbors.size() > node -> max_connections) {
                        std::vector<HNSWNode*> eNewConn = select_neighbors_simple(node, node->neighbors,node -> max_connections, l_c);
                        for (HNSWNode* selected_neighbor : eNewConn) {
                            selected_neighbor -> neighbors.push_back(node);
                        }
                        node -> neighbors = eNewConn;
                    }
                }
            }

            
            


        }

        /// ALGORITHM 2
        /// @param q: the query vector
        /// @param ep: enter point ep
        /// @param ef: number of nearest to q elements to return
        /// @param l_c: layer number
        /// @return ef closest neighbors to q
        std::vector<HNSWNode*> search_layer(DataPoint* q, HNSWNode* ep, int ef, uint64_t l_c) {
            UNUSED(q);
            UNUSED(ep);
            UNUSED(ef);
            UNUSED(l_c);


            HNSWNode* curr_ep = ep;

            while (curr_ep -> curr_level != l_c) {
                curr_ep = curr_ep -> lower_node;
            }


            // set of visited elements - should I track
            std::set<int> v;
            v.insert(curr_ep -> data_point -> id);

            // set of candidates that we have expanded
            std::set<HNSWNode*> C;
            C.insert(curr_ep);
            
            // list of the nearest neighbors so far in our search
            std::vector<HNSWNode*> W;
            W.push_back(curr_ep);

            while (C.size() > 0) {
                int min_dist = INT_MAX;
                int max_dist = INT_MIN;
                HNSWNode *c = nullptr;
                HNSWNode *f = nullptr;
                for (HNSWNode* candidate: C) {
                    if (euclidean_distance(candidate -> data_point -> vector,q->vector) < min_dist) {
                        c = candidate;
                    }

                    if (euclidean_distance(candidate -> data_point -> vector,q->vector) > max_dist) {
                        f = candidate;
                    }
                }

                if (euclidean_distance(c -> data_point -> vector, q->vector) > euclidean_distance(f -> data_point -> vector, q->vector)) {
                    break;
                }

                
                for (HNSWNode* e : c -> neighbors) {

                    // e element of v
                    if (v.find(e->data_point->id) == v.end()) {
                        v.insert(e -> data_point -> id);
                    }


                    // should create helper for getting the farthest value
                    f = nullptr;
                    max_dist = INT_MAX;
                    for (HNSWNode* candidate: C) {
                        if (euclidean_distance(candidate -> data_point -> vector,q->vector) > max_dist) {
                            f = candidate;
                        }
                    }
                    if ((euclidean_distance(e->data_point->vector,q->vector) < euclidean_distance(f->data_point->vector,q->vector)) || W.size() < static_cast<size_t>(ef)) {
                        C.insert(e);
                        W.push_back(e);
                        if (W.size() > static_cast<size_t>(ef)) {
                            // optmize remove by making W a map or set or smthg and pushing values to a vector at the end
                            W.erase(std::remove(W.begin(), W.end(), f), W.end());
                        }
                    }

                }


            }

            return W;
        }

        /// ALGORITHM 3
        /// @param q: base element
        /// @param C:  candidate elements
        /// @param M: number of neighbors to return
        std::vector<HNSWNode*> select_neighbors_simple(HNSWNode* q, std::vector<HNSWNode*> C, uint64_t M, int l_c) {
            UNUSED(q);
            UNUSED(C);
            UNUSED(M);
            UNUSED(l_c);


            auto comp = [](const std::pair<HNSWNode*, float>& a, const std::pair<HNSWNode*, float>& b) {
                return a.second > b.second;
            };


            std::priority_queue<std::pair<HNSWNode*, float>, std::vector<std::pair<HNSWNode*, float>>, decltype(comp)> priority_queue(comp);

            for (HNSWNode* node : C) {
                priority_queue.push({node, euclidean_distance(node->data_point->vector,q->data_point->vector)});
            }

            std::vector<HNSWNode*> ret;

            while (!priority_queue.empty()) {
                ret.push_back(priority_queue.top().first);
                priority_queue.pop();
            }

            return ret;
            
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
        std::vector<HNSWNode*> knn_search(DataPoint* q, uint64_t k, uint64_t ef) {
            UNUSED(q);
            UNUSED(k);
            UNUSED(ef);

            std::vector<HNSWNode*> W;
            

            HNSWNode* ep = entry_point;

            for (int l_c = L; l_c >= 1; l_c--) {
                W = search_layer(q,ep,1,l_c);
                ep=W[0];
            }

            // At entry point for the bottom most node, so just search that
            W = search_layer(q,ep,ef,0);

            auto comp = [](const std::pair<HNSWNode*, float>& a, const std::pair<HNSWNode*, float>& b) {
                return a.second > b.second;
            };
            std::priority_queue<std::pair<HNSWNode*, float>, std::vector<std::pair<HNSWNode*, float>>, decltype(comp)> priority_queue(comp);

            for (HNSWNode* node : W) {
                priority_queue.push({node, euclidean_distance(node->data_point->vector,q->vector)});
            }

            std::vector<HNSWNode*> ret;

            for (size_t i = 0; i < k; i++) {
                ret.push_back(priority_queue.top().first);
                priority_queue.pop();
            }

            return ret;
        }   


        
};

void listDatasets(H5::H5File &file) {
    std::string groupName = "/";
    H5::Group group = file.openGroup(groupName);
    hsize_t num_objects = group.getNumObjs();

    std::cout << "Listing datasets and groups from the root" << std::endl;

    for (hsize_t i = 0; i < num_objects; i++) {
        std::string objName = group.getObjnameByIdx(i);
        H5G_obj_t objType = group.getObjTypeByIdx(i);

        if (objType == H5G_DATASET) {
            std::cout <<"DATASET: " << objName << std::endl;
        }
    }

}


// Reference: Loading Hdf5 dataset into C++ vector - https://stackoverflow.com/questions/25568446/loading-data-from-hdf5-to-vector-in-c
H5::DataSet loadDataset(std::string dataset_split = "train") {
    // will set this parameter dynamically in future iterations
    
    std::string glove_file = "glove-25-angular.hdf5";

    

    // get the datset reference
    H5::H5File file(glove_file, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet(dataset_split);
   

    // Get the dataset size
    H5::DataSpace dataspace = dataset.getSpace(); // dataspace objects just helps to do eqv of np.shape, explroing data
    hsize_t dims[2];
    int ndims = dataspace.getSimpleExtentDims(dims, NULL);
    UNUSED(ndims);

    // Get dataset size
    std::cout << dataset_split << std::endl;
    std::cout << "DATSET DIMENSIONS: " << dims[0] << " x " << dims[1] << std::endl;
    std::cout << "___________________________________" << std::endl;


    // rather than loading the entire dataset into the buffer in one go, we read a certain number of rows per chunk
    hsize_t chunk_size = 1000;
    hsize_t num_chunks = (dims[0] + chunk_size - 1) / chunk_size;

    std::vector<float> single_chunk_buffer(chunk_size * dims[1]);


    /// Dump chunk to CSV file
    std::ofstream csvf(dataset_split + ".csv");

    for (hsize_t i = 0; i < num_chunks; i++) {
        hsize_t offset[2] = {i * chunk_size, 0};
        hsize_t count[2] = {std::min(chunk_size, dims[0] - i * chunk_size), dims[1]};
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

        // Defining memory dataspace matches the size of the chunk
        H5::DataSpace memspace(2, count);
        dataset.read(single_chunk_buffer.data(), H5::PredType::NATIVE_FLOAT, memspace, dataspace);


        


        for (hsize_t r=0; r < count[0]; r++) {
            for (hsize_t c=0; c < count[1]; c++) {
                csvf << single_chunk_buffer[r * dims[1] + c];
                if (c < count[1] - 1) {
                    csvf <<",";
                }
            }
            csvf << std::endl;
        }

    }


    csvf.close();
    file.close();

    return dataset;
}

// defaults
uint64_t default_k = 3;
const uint64_t VALIDATION_SIZE = 10000;


int main(int argc, char* argv[]) {
    UNUSED(argc);
    UNUSED(argv);

    return 0;

    std::string glove_file = "glove-25-angular.hdf5";
    
    // Listing the different datasets that could be used
    H5::H5File root_explore(glove_file, H5F_ACC_RDONLY);
    listDatasets(root_explore);
    root_explore.close();
    std::cout << "________________________________________" << std::endl;

}