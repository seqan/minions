#include <gtest/gtest.h>

#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "compare.h"

TEST(comparison, small_example)
{
    std::string expected
    {
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

    EXPECT_EQ(expected.substr(0,9), results[0]);
    EXPECT_EQ("0", results[4]);
    EXPECT_EQ(expected.substr(47,14), results[6]);
    EXPECT_EQ("0", results[11]);
    EXPECT_EQ(expected.substr(94,14), results[13]);
    EXPECT_EQ("0", results[18]);
}
