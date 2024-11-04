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


int main(int argc, char* argv[]) {
    UNUSED(argc);
    UNUSED(argv);

    return 0;

}