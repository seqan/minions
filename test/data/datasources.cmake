cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# copies file to <build>/data/*
declare_datasource (FILE example1.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/example1.fasta
                    URL_HASH SHA256=6e30fc35f908a36fe0c68a7a35c47f51f9570da16622fb0c072a20e6a9ba5b3e)
