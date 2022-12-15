# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if("/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/search.fasta" STREQUAL "")
  message(FATAL_ERROR "LOCAL can't be empty")
endif()

if(NOT EXISTS "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/search.fasta")
  message(FATAL_ERROR "File not found: /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/search.fasta")
endif()

if("SHA256" STREQUAL "")
  message(WARNING "File will not be verified since no URL_HASH specified")
  return()
endif()

if("abed0af7e29a07f5964239be77b46c427369f88cbd2e3c677f763cd7d94f9e4a" STREQUAL "")
  message(FATAL_ERROR "EXPECT_VALUE can't be empty")
endif()

message(STATUS "verifying file...
     file='/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/search.fasta'")

file("SHA256" "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/search.fasta" actual_value)

if(NOT "${actual_value}" STREQUAL "abed0af7e29a07f5964239be77b46c427369f88cbd2e3c677f763cd7d94f9e4a")
  message(FATAL_ERROR "error: SHA256 hash of
  /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/search.fasta
does not match expected value
  expected: 'abed0af7e29a07f5964239be77b46c427369f88cbd2e3c677f763cd7d94f9e4a'
    actual: '${actual_value}'
")
endif()

message(STATUS "verifying file... done")
