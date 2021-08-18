#include <gtest/gtest.h>

#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "compare.h"

TEST(minions, small_example)
{
    testing::internal::CaptureStdout();
    do_comparison({DATADIR"example1.fasta"});
    std::string std_cout = testing::internal::GetCapturedStdout();

    std::string expected_compression
    {
        "COMPRESSION            	Minimum	Mean	StdDev	Maximum\n"
        "kmer_hash (19)         	1233143	1233143	0	1233143\n"
        "minimiser_hash (19, 19)	903366	903366	0	903366\n"
        "minimiser_hash (19, 23)	303050	303050	0	303050\n"

    };

    std::istringstream iss(std_cout);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>());

    EXPECT_EQ(expected_compression.substr(0,53), std_cout.substr(0,53));
    EXPECT_EQ(expected_compression.substr(52,9), results[5]);
    EXPECT_EQ("0", results[9]);
    EXPECT_EQ(expected_compression.substr(102,14), results[11]);
    EXPECT_EQ("0", results[16]);
    EXPECT_EQ(expected_compression.substr(149,14), results[18]);
    EXPECT_EQ("0", results[23]);

    std::string expected_speed{"SPEED                  	Minimum	Mean	StdDev	Maximum\n"};
    std_cout = std_cout.substr(196, 52);

    EXPECT_EQ(expected_speed, std_cout);
}
