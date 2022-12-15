# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/_deps/gtest_fetch_content-src"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/_deps/gtest_fetch_content-build"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/_deps/gtest_fetch_content-subbuild/gtest_fetch_content-populate-prefix"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/_deps/gtest_fetch_content-subbuild/gtest_fetch_content-populate-prefix/tmp"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/_deps/gtest_fetch_content-subbuild/gtest_fetch_content-populate-prefix/src/gtest_fetch_content-populate-stamp"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/_deps/gtest_fetch_content-subbuild/gtest_fetch_content-populate-prefix/src"
  "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/_deps/gtest_fetch_content-subbuild/gtest_fetch_content-populate-prefix/src/gtest_fetch_content-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/cmake-build-debug/_deps/gtest_fetch_content-subbuild/gtest_fetch_content-populate-prefix/src/gtest_fetch_content-populate-stamp/${subDir}")
endforeach()
