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
#include <sstream>
#include <chrono>
#include <filesystem>
#include "H5Cpp.h"


#define UNUSED(p)  ((void)(p))


struct DataPoint {
    public:
        int id;
        std::vector<float> vector;

        DataPoint(int id, std::vector<float> vector) : id(id), vector(vector) {}
    

};

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


int main(int argc, char* argv[]) {
    UNUSED(argc);
    UNUSED(argv);

    // Will be set dynamically in the future to select SIFT benchmarks
    std::string glove_file = "glove-25-angular.hdf5";
    
    // Listing the different datasets that could be used
    H5::H5File root_explore(glove_file, H5F_ACC_RDONLY);
    listDatasets(root_explore);
    root_explore.close();
    std::cout << "________________________________________" << std::endl;


    /** Data to insert */
    H5::DataSet train_data = loadDataset("train");


    /** Validation points and ground truth */
    H5::DataSet test_data = loadDataset("test"); 
    H5::DataSet neighbors = loadDataset("neighbors"); 
    H5::DataSet distances = loadDataset("distances"); 




    /*
    * Load Train Points as a flat array of Data Points 
    * Time with std::chrono, lecture 3 slides
    */


    auto start = std::chrono::high_resolution_clock::now();

    std::vector<DataPoint*> test_points;
    UNUSED(test_points);
    std::ifstream file("tr.csv");

    
    size_t ctr = 0;
    std::string line;
    
    while(std::getline(file,line)) {
        std::vector<float> row;
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss,value,',')) {
            row.push_back(std::stof(value));
        }


        DataPoint* point = new DataPoint(ctr,row);
        test_points.push_back(point);

        ++ctr;
    }

    file.close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    /// Getting current timestamp from stack overflow answer: https://stackoverflow.com/questions/6012663/get-unix-timestamp-with-c
    int64_t timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    UNUSED(timestamp);

    
    //** Writing the results */
    std::string csv_name= "baseline_benchmark.csv";
    bool csv_created = std::filesystem::exists(csv_name); // https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exists-using-standard-c-c11-14-17-c
    std::ofstream output_csv(csv_name, std::ios::app);

    if (!csv_created) {
        output_csv << "timestamp,id,time\n";
    }

    output_csv << timestamp  << ",all," << elapsed.count();







    /** Uncomment to verify all elements are inserted */

    // for (DataPoint* point : test_points) {
    //     std::vector<float> row = point -> vector;
    //     std::cout << "ID: " << point -> id << std::endl;
    //     for (float ele : row) {
    //         std::cout << ele << ", ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << std::endl;
    // }




    return 0;

} 