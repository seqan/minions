#include <gtest/gtest.h>

#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "compare.h"

TEST(minions, small_example)
{
    std::string expected
    {
        "SPEED                  	Minimum	Mean	StdDev	Maximum\n"
        "kmer_hash (19)         	159493	159493	0	159493\n"
        "minimiser_hash (19, 19)	241282	241282	0	241282\n"
        "minimiser_hash (19, 23)	157535	157535	0	157535\n"

    };
    testing::internal::CaptureStdout();
    do_comparison({DATADIR"example1.fasta"});
    std::string std_cout = testing::internal::GetCapturedStdout();
    std::istringstream iss(std_cout);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>());

    EXPECT_EQ(expected.substr(0,53), std_cout.substr(0,53));
    EXPECT_EQ(expected.substr(52,9), results[5]);
    EXPECT_EQ("0", results[9]);
    EXPECT_EQ(expected.substr(99,14), results[11]);
    EXPECT_EQ("0", results[16]);
    EXPECT_EQ(expected.substr(146,14), results[18]);
    EXPECT_EQ("0", results[23]);
}
