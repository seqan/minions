cmake_minimum_required (VERSION 3.8)

add_api_test (comparison_test.cpp)
target_use_datasources (comparison_test FILES example1.fasta example.ibf minimiser_hash_19_19_example1.out)

add_api_test (minimiser_distance_test.cpp)

add_api_test (modmer_test.cpp)
add_api_test (modmer_hash_test.cpp)
add_api_test (modmer_hash_distance_test.cpp)
