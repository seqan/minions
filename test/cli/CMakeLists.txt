cmake_minimum_required (VERSION 3.8)

add_cli_test (minions_options_test.cpp)
add_cli_test (minions_accuracy_test.cpp FILES example.ibf expected_search_result.out minimiser_hash_19_19_example1.out)
add_cli_test (minions_counts_test.cpp FILES example1.fasta)
add_cli_test (minions_coverage_test.cpp FILES example1.fasta)
add_cli_test (minions_speed_test.cpp FILES example1.fasta)
