cmake_minimum_required (VERSION 3.8)

add_api_test (comparison_test.cpp)
target_use_datasources (comparison_test FILES example1.fasta)

add_api_test (minimiser_distance_test.cpp)

add_api_test (modmer_test.cpp)
add_api_test (modmer_hash_test.cpp)
add_api_test (modmer_hash_distance_test.cpp)

add_api_test (syncmer_test.cpp)
