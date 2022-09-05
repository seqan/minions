#include <forward_list>
#include <list>
#include <type_traits>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

#include <seqan3/test/expect_range_eq.hpp>

#include <gtest/gtest.h>

#include "minstrobe_hash.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<uint64_t>;

static constexpr seqan3::shape ungapped_shape = seqan3::ungapped{4};
static constexpr seqan3::shape gapped_shape = 0b1001_shape;
static constexpr auto ungapped_view = minstrobe2_hash(ungapped_shape, 2, 7, seqan3::seed{0});
static constexpr auto gapped_view = minstrobe2_hash(gapped_shape, 2, 7, seqan3::seed{0});

static constexpr auto ungapped3_view = minstrobe3_hash(ungapped_shape, 1, 6, seqan3::seed{0});
static constexpr auto gapped3_view = minstrobe3_hash(gapped_shape, 1, 6, seqan3::seed{0});

template <typename T>
class minstrobe_hash_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;

TYPED_TEST_SUITE(minstrobe_hash_view_properties_test, underlying_range_types, );

class minstrobe_hash_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAA"_dna4};
    result_t result1{0,0,0,0}; // Same result for ungapped and gapped

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    //            ungapped minstrobes: ACGGCGAC, CGGCACGT, GGCGACGT, GCGAACGT, CGACACGT, GACGCGTT
    //              gapped minstrobes: A--GC--C, C--CA--T, G--GA--T, G--AA--T, C--CA--T, G--GC--T
    //  stop at T ungapped minstrobes: ACGGCGAC
    //    stop at T gapped minstrobes: A--GC--C
    // start at A ungapped minstrobes:                               GCGAACGT, CGACACGT, GACGCGTT
    result_t result3_ungapped{6753, 26907, 42523, 38939, 24859, 34415};
    result_t result3_gapped{37, 83, 163, 131, 83, 167};
    result_t result3_ungapped_stop{6753};
    result_t result3_gapped_stop{37};
    result_t result3_ungapped_start{38939, 24859, 34415};
    result_t result3_gapped_start{131, 83, 167};

    result_t result3_1{0,0,0}; // Same result for ungapped and gapped

    // ACGGCGACGTTTAG
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    //            ungapped minstrobes: ACGGCGGCACGT, CGGCCGACACGT, GGCGCGACACGT, GCGAACGTCGTT, CGACACGTGTTT
    //              gapped minstrobes: A--GC--CA--T, C--CC--CA--T, G--GC--CA--T, G--AA--TC--T, C--CA--TG--T
    // start at A ungapped minstrobes:                               GGCGCGACACGT, GCGAACGTCGTT, CGACACGTCGTT
    //   start at A gapped minstrobes:                               G--GC--CA--T, G--AA--TC--T, C--CA--TG--T
    result_t result3_3_ungapped{1730843,6906139,10903835,9968495,6364095};
    result_t result3_3_gapped{595,1363,2643,2103,1339};
    result_t result3_3_ungapped_start{9968495,6364095};
    result_t result3_3_gapped_start{2103,1339};
};

template <typename adaptor_t>
void compare_types(adaptor_t v)
{
    EXPECT_TRUE(std::ranges::input_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::view<decltype(v)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v), size_t>));
}

TYPED_TEST(minstrobe_hash_view_properties_test, different_input_ranges)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{7009, 27931, 46619, 55323, 24859, 34415};
    result_t gapped{53, 83, 163, 195, 83, 167};
    EXPECT_RANGE_EQ(ungapped, text | ungapped_view);
    EXPECT_RANGE_EQ(gapped, text | gapped_view);

    result_t ungapped3{7009, 27931, 46619, 55323, 24859, 34415};
    result_t gapped3{773, 1283, 2563, 3075, 1283, 2567};
    //EXPECT_RANGE_EQ(ungapped3, text | ungapped3_view);
    //EXPECT_RANGE_EQ(gapped3, text | gapped3_view);
}

TEST_F(minstrobe_hash_test, ungapped)
{
    EXPECT_RANGE_EQ(result1, text1 | ungapped_view);
    EXPECT_RANGE_EQ(result3_ungapped, text3 | ungapped_view);
    EXPECT_NO_THROW(text1 | minstrobe2_hash(ungapped_shape,3,6));
    EXPECT_THROW((text3 | minstrobe2_hash(ungapped_shape,3,2)), std::invalid_argument);

    EXPECT_RANGE_EQ(result3_1, text1 | ungapped3_view);
    EXPECT_RANGE_EQ(result3_3_ungapped, text3 | ungapped3_view);
    EXPECT_NO_THROW(text1 | minstrobe3_hash(ungapped_shape,3,6));
    EXPECT_THROW((text3 | minstrobe3_hash(ungapped_shape,3,2)), std::invalid_argument);
}

TEST_F(minstrobe_hash_test, gapped)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_view);
    EXPECT_RANGE_EQ(result3_gapped, text3 | gapped_view);
    EXPECT_NO_THROW(text1 | minstrobe2_hash(gapped_shape, 2,5));
    EXPECT_THROW((text3 | minstrobe2_hash(gapped_shape, 2,1)), std::invalid_argument);

    EXPECT_RANGE_EQ(result3_1, text1 | gapped3_view);
    EXPECT_RANGE_EQ(result3_3_gapped, text3 | gapped3_view);
    EXPECT_NO_THROW(text1 | minstrobe3_hash(gapped_shape, 2,5));
    EXPECT_THROW((text3 | minstrobe3_hash(gapped_shape, 2,1)), std::invalid_argument);
}

TEST_F(minstrobe_hash_test, combinability)
{
    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_stop, text3 | stop_at_t | ungapped_view);
    EXPECT_RANGE_EQ(result3_gapped_stop, text3 | stop_at_t | gapped_view);

    auto start_at_a = std::views::drop(3);
    EXPECT_RANGE_EQ(result3_ungapped_start, text3 | start_at_a | ungapped_view);
    EXPECT_RANGE_EQ(result3_gapped_start, text3 | start_at_a | gapped_view);
    EXPECT_RANGE_EQ(result3_3_ungapped_start, text3 | start_at_a | ungapped3_view);
    EXPECT_RANGE_EQ(result3_3_gapped_start, text3 | start_at_a | gapped3_view);
}
