# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

if("/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example.ibf" STREQUAL "")
  message(FATAL_ERROR "LOCAL can't be empty")
endif()

if(NOT EXISTS "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example.ibf")
  message(FATAL_ERROR "File not found: /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example.ibf")
endif()

if("SHA256" STREQUAL "")
  message(WARNING "File will not be verified since no URL_HASH specified")
  return()
endif()

if("8d18ce55fdbb78acbd4f44d5414de1b55ac9964e00e391a5fd12bcc1622b1c6c" STREQUAL "")
  message(FATAL_ERROR "EXPECT_VALUE can't be empty")
endif()

message(STATUS "verifying file...
     file='/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example.ibf'")

file("SHA256" "/Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example.ibf" actual_value)

if(NOT "${actual_value}" STREQUAL "8d18ce55fdbb78acbd4f44d5414de1b55ac9964e00e391a5fd12bcc1622b1c6c")
  message(FATAL_ERROR "error: SHA256 hash of
  /Users/matanatmammadli/Desktop/Bachelorarbeit/minions/test/data/example.ibf
does not match expected value
  expected: '8d18ce55fdbb78acbd4f44d5414de1b55ac9964e00e391a5fd12bcc1622b1c6c'
    actual: '${actual_value}'
")
endif()

message(STATUS "verifying file... done")
