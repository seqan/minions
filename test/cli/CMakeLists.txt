cmake_minimum_required (VERSION 3.8)

add_cli_test (minions_options_test.cpp)
add_cli_test (minions_accuracy_test.cpp FILES example.ibf)
add_cli_test (minions_coverage_test.cpp FILES example1.fasta)
add_cli_test (minions_speed_test.cpp FILES example1.fasta)
