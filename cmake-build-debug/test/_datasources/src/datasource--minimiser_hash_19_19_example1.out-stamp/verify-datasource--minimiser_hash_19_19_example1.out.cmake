# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if("/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/minimiser_hash_19_19_example1.out" STREQUAL "")
  message(FATAL_ERROR "LOCAL can't be empty")
endif()

if(NOT EXISTS "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/minimiser_hash_19_19_example1.out")
  message(FATAL_ERROR "File not found: /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/minimiser_hash_19_19_example1.out")
endif()

if("SHA256" STREQUAL "")
  message(WARNING "File will not be verified since no URL_HASH specified")
  return()
endif()

if("8086779dc7fb37a81f20a3d202d2f9c5a39c4611f67f2c3881e4ba394deef9e6" STREQUAL "")
  message(FATAL_ERROR "EXPECT_VALUE can't be empty")
endif()

message(STATUS "verifying file...
     file='/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/minimiser_hash_19_19_example1.out'")

file("SHA256" "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/minimiser_hash_19_19_example1.out" actual_value)

if(NOT "${actual_value}" STREQUAL "8086779dc7fb37a81f20a3d202d2f9c5a39c4611f67f2c3881e4ba394deef9e6")
  message(FATAL_ERROR "error: SHA256 hash of
  /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/minimiser_hash_19_19_example1.out
does not match expected value
  expected: '8086779dc7fb37a81f20a3d202d2f9c5a39c4611f67f2c3881e4ba394deef9e6'
    actual: '${actual_value}'
")
endif()

message(STATUS "verifying file... done")
