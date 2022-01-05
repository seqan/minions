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
