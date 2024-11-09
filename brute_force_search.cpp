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
int main(int argc, char* argv[]) {
    UNUSED(argc);
    UNUSED(argv);



    // will set this parameter dynamically in future iterations
    const std::string glove_file = "glove-25-angular.hdf5";
    
    // Listing the different datasets that could be used
    H5::H5File root_explore(glove_file, H5F_ACC_RDONLY);
    listDatasets(root_explore);
    root_explore.close();


    const std::string dataset_split = "train";

    

    // get the datset reference
    H5::H5File file(glove_file, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet(dataset_split);
   

    // Get the dataset size
    H5::DataSpace dataspace = dataset.getSpace(); // dataspace objects just helps to do eqv of np.shape, explroing data
    hsize_t dims[2];
    int ndims = dataspace.getSimpleExtentDims(dims, NULL);
    UNUSED(ndims);
    std::cout << "DATSET DIMENSIONS: " << dims[0] << " x " << dims[1] << std::endl;


    // rather than loading the entire dataset into the buffer in one go, we read a certain number of rows per chunk
    hsize_t chunk_size = 1000;
    hsize_t num_chunks = (dims[0] + chunk_size - 1) / chunk_size;

    std::vector<float> single_chunk_buffer(chunk_size * dims[1]);

    for (hsize_t i = 0; i < num_chunks; i++) {
        hsize_t offset[2] = {i * chunk_size, 0};
        hsize_t count[2] = {std::min(chunk_size, dims[0] - i * chunk_size), dims[1]};
        dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);

        // Defining memory dataspace matches the size of the chunk
        H5::DataSpace memspace(2, count);
        dataset.read(single_chunk_buffer.data(), H5::PredType::NATIVE_FLOAT, memspace, dataspace);
    }


    




    return 0;

} 