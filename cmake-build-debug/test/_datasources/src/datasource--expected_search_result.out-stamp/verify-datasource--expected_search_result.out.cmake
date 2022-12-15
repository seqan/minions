# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if("/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/expected_search_result.out" STREQUAL "")
  message(FATAL_ERROR "LOCAL can't be empty")
endif()

if(NOT EXISTS "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/expected_search_result.out")
  message(FATAL_ERROR "File not found: /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/expected_search_result.out")
endif()

if("SHA256" STREQUAL "")
  message(WARNING "File will not be verified since no URL_HASH specified")
  return()
endif()

if("7d51b8ac01dd8020bcb88353a5e93b9583bd352d10d38e229d669bbc196af898" STREQUAL "")
  message(FATAL_ERROR "EXPECT_VALUE can't be empty")
endif()

message(STATUS "verifying file...
     file='/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/expected_search_result.out'")

file("SHA256" "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/expected_search_result.out" actual_value)

if(NOT "${actual_value}" STREQUAL "7d51b8ac01dd8020bcb88353a5e93b9583bd352d10d38e229d669bbc196af898")
  message(FATAL_ERROR "error: SHA256 hash of
  /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/expected_search_result.out
does not match expected value
  expected: '7d51b8ac01dd8020bcb88353a5e93b9583bd352d10d38e229d669bbc196af898'
    actual: '${actual_value}'
")
endif()

message(STATUS "verifying file... done")
