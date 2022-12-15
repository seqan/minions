# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if("/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example1.fasta" STREQUAL "")
  message(FATAL_ERROR "LOCAL can't be empty")
endif()

if(NOT EXISTS "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example1.fasta")
  message(FATAL_ERROR "File not found: /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example1.fasta")
endif()

if("SHA256" STREQUAL "")
  message(WARNING "File will not be verified since no URL_HASH specified")
  return()
endif()

if("e7236e7b86303d84a86a4454044125005433857416183cdaac0f4cdf7ac34e06" STREQUAL "")
  message(FATAL_ERROR "EXPECT_VALUE can't be empty")
endif()

message(STATUS "verifying file...
     file='/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example1.fasta'")

file("SHA256" "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example1.fasta" actual_value)

if(NOT "${actual_value}" STREQUAL "e7236e7b86303d84a86a4454044125005433857416183cdaac0f4cdf7ac34e06")
  message(FATAL_ERROR "error: SHA256 hash of
  /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example1.fasta
does not match expected value
  expected: 'e7236e7b86303d84a86a4454044125005433857416183cdaac0f4cdf7ac34e06'
    actual: '${actual_value}'
")
endif()

message(STATUS "verifying file... done")
