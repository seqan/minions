cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/*
declare_datasource (FILE example1.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/example1.fasta
                    URL_HASH SHA256=e7236e7b86303d84a86a4454044125005433857416183cdaac0f4cdf7ac34e06)
