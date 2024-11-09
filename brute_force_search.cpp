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
#include "H5Cpp.h"


#define UNUSED(p)  ((void)(p))




float angular_distance(const std::vector<float> &v1, const std::vector<float> &v2) {
    UNUSED(v1);
    UNUSED(v2);

    float dot_product = 0.0;
    float v1_normalize = 0.0;
    float v2_normalize = 0.0;

    for (uint64_t i = 0; i < v1.size(); i++) {
        dot_product += v1[i] * v2[i];
        v1_normalize= v1[i] * v1[i];
        v2_normalize= v2[i] * v2[i];
    }

    return 1.0 - (dot_product /(sqrt(v1_normalize) * sqrt(v2_normalize)));

}


int main(int argc, char* argv[]) {
    UNUSED(argc);
    UNUSED(argv);


    // will set this parameter dynamically in future iterations
    const std::string glove_file = "glove-25-angular.hdf5";
    const std::string dataset_split = "train";

    H5::H5File file(glove_file, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet(dataset_split);
    H5::DataSpace dataspace = dataset.getSpace(); // dataspace objects just helps to do eqv of np.shape, explroing data

    hsize_t dims[2];
    int ndims = dataspace.getSimpleExtentDims(dims, NULL);
    UNUSED(ndims);
    std::cout << "DATSET DIMENSIONS: " << dims[0] << " x " << dims[1] << std::endl;







    return 0;

} 