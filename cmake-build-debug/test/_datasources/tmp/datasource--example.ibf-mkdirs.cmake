# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/_datasources/src/datasource--example.ibf"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/_datasources/src/datasource--example.ibf-build"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/_datasources"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/_datasources/tmp"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/_datasources/src/datasource--example.ibf-stamp"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/_datasources/src"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/_datasources/src/datasource--example.ibf-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/test/_datasources/src/datasource--example.ibf-stamp/${subDir}")
endforeach()
