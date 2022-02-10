cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/*
declare_datasource (FILE example.ibf
                    URL ${CMAKE_SOURCE_DIR}/test/data/example.ibf
                    URL_HASH SHA256=8d18ce55fdbb78acbd4f44d5414de1b55ac9964e00e391a5fd12bcc1622b1c6c)
declare_datasource (FILE minimiser_hash_19_19_example1.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/minimiser_hash_19_19_example1.out
                    URL_HASH SHA256=8086779dc7fb37a81f20a3d202d2f9c5a39c4611f67f2c3881e4ba394deef9e6)
declare_datasource (FILE example1.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/example1.fasta
                    URL_HASH SHA256=e7236e7b86303d84a86a4454044125005433857416183cdaac0f4cdf7ac34e06)
declare_datasource (FILE search.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/search.fasta
                    URL_HASH SHA256=abed0af7e29a07f5964239be77b46c427369f88cbd2e3c677f763cd7d94f9e4a)
declare_datasource (FILE expected_search_result.out
                    URL ${CMAKE_SOURCE_DIR}/test/data/expected_search_result.out
                    URL_HASH SHA256=7d51b8ac01dd8020bcb88353a5e93b9583bd352d10d38e229d669bbc196af898)
