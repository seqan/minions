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
    EXPECT_EQ(expected.substr(0,14), std_cout.substr(0,14));
    EXPECT_EQ(expected[38], std_cout[38]);
    EXPECT_EQ(expected.substr(46,24), std_cout.substr(46,24));
    EXPECT_EQ(expected[84], std_cout[84]);
    seqan3::debug_stream << expected.substr(74,24) << " "<<expected[84] << " " <<expected[76] << " " << expected[78] << "\n";
    EXPECT_EQ(expected.substr(93,24), std_cout.substr(93,24));
    EXPECT_EQ(expected[132], std_cout[132]);
    seqan3::debug_stream <<expected.substr(93,24) << " "<< expected[132] << "\n";
}
