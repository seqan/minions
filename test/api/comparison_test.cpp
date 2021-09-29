#include <gtest/gtest.h>

#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "compare.h"

TEST(minions, small_example)
{
    range_arguments args{};
    args.kmers = true;
    args.k_size = 19;
    std::string expected{"kmer_hash (19)         	159493	159493	0	159493\n"};
    testing::internal::CaptureStdout();
    std::filesystem::path path_out = std::filesystem::temp_directory_path();
    do_comparison({DATADIR"example1.fasta"}, args, path_out);
    std::string std_cout = testing::internal::GetCapturedStdout();
    std::istringstream iss(std_cout);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>());


    EXPECT_EQ(expected.substr(0,9), results[0]);
    EXPECT_EQ("0", results[4]);
}
