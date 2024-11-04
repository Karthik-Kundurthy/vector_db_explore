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

// Manage the entire dataset: SIFT and GLOVE
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



};


int main(int argc, char* argv[]) {
    UNUSED(argc);
    UNUSED(argv);

    return 0;

}