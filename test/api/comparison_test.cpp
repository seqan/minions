#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>

#include "compare.h"

TEST(minions, small_example)
{
    range_arguments args{};
    args.name = kmer;
    args.k_size = 19;
    args.shape = seqan3::ungapped{19};
    std::string expected{"kmer_hash_19         	159493	159493	0	159493\n"};
    args.path_out = std::filesystem::temp_directory_path();
    do_comparison({DATADIR"example1.fasta"}, args);

    std::ifstream infile;
    std::string line;
    infile.open(std::string{args.path_out} + std::string{args.path_out}  + "kmer_hash_19_speed_compression.out");
    if(infile.is_open())
    {
        while(std::getline(infile, line))
        {
            std::istringstream iss(line);
            std::vector<std::string> results((std::istream_iterator<std::string>(iss)),
                                     std::istream_iterator<std::string>());
            EXPECT_EQ(expected.substr(0,12), results[0]);
            EXPECT_EQ("0", results[3]);
        }
    }
    infile.close();
}

TEST(minions, accuracy_binary_file)
{
    range_arguments args{};
    args.name = minimiser;
    args.k_size = 19;
    args.w_size = seqan3::window_size{19};
    args.shape = seqan3::ungapped{19};
    args.path_out = std::filesystem::path{std::string{std::filesystem::temp_directory_path()} + "/"};
    do_accuracy({DATADIR"minimiser_hash_19_19_example1.out"}, args, 1000000, 1);

    seqan3::interleaved_bloom_filter ibf{};
    load_ibf(ibf, std::string{args.path_out} + "minimiser_hash_19_19.ibf");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(1, 1);
    auto & res = agent.bulk_contains(39030638997);
    EXPECT_RANGE_EQ(expected_result,  res);
}

TEST(minions, accuracy_existing_ibf)
{
    range_arguments args{};
    args.name = minimiser;
    args.k_size = 19;
    args.w_size = seqan3::window_size{19};
    args.shape = seqan3::ungapped{19};
    args.path_out = std::filesystem::path{std::string{std::filesystem::temp_directory_path()} + "/"};
    do_accuracy({DATADIR"example.ibf"}, args, 1000000, 1);

    seqan3::interleaved_bloom_filter ibf{};
    load_ibf(ibf, DATADIR"example.ibf");
    auto agent = ibf.membership_agent();

    std::vector<bool> expected_result(1, 1);
    auto & res = agent.bulk_contains(39030638997);
    EXPECT_RANGE_EQ(expected_result,  res);
}
