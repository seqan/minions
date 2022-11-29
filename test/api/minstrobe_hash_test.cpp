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
static constexpr auto ungapped_view = minstrobe2_hash(ungapped_shape, 1, 6, seqan3::seed{0});
static constexpr auto gapped_view = minstrobe2_hash(gapped_shape, 1, 6, seqan3::seed{0});

static constexpr auto ungapped3_view = minstrobe3_hash(ungapped_shape, 1, 6, seqan3::seed{0});
static constexpr auto gapped3_view = minstrobe3_hash(gapped_shape, 1, 6, seqan3::seed{0});

template <typename T>
class minstrobe_hash_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                //seqan3::bitpacked_sequence<seqan3::dna4>,
                                                //seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;

TYPED_TEST_SUITE(minstrobe_hash_view_properties_test, underlying_range_types, );

class minstrobe_hash_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAA"_dna4};
    result_t result1{0,0,0}; // Same result for ungapped and gapped

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4}; // rev complement: CTAA ACGTCGCCGT
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //        reverse ungapped Hashes: 91,       150,      101,      217,      182,      109,      27,   6,     1,   192,  176
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    //          reverse gapped Hashes: 7,        10,        5,       13,       10,       5,        3,    2,     1,   12,    4
    //            ungapped minstrobes: ACGGACGT, aacgcgcc, aaaccgcc, GCGACGTT, CGACGTTT
    //              gapped minstrobes: A--GA--T, a--gc--c, a--cc--c, G--AC--T, c--aa--t
    // start at A ungapped minstrobes:                               GCGACGTT, CGACGTTT
    result_t result3_ungapped{6683, 1637, 357, 39023, 25023};
    result_t result3_gapped{35, 37, 21, 135, 67};
    result_t result3_ungapped_start{39023, 25023};
    result_t result3_gapped_start{135, 67};

    std::vector<seqan3::dna4> text1_3{"AAAAAAAAAAAAAAAA"_dna4};
    result_t result3_1{0}; // Same result for ungapped and gapped
    std::vector<seqan3::dna4> text3_3{"ACGGCGACGTTTAGGC"_dna4};
    // ACGGCGACGTTTAG
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG, TAGG, AGGC
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242,  202,    41
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14,   16,      1
    //            ungapped minstrobes: ACGGACGTAGGC
    //              gapped minstrobes: A--GA--TA--C
    result_t result3_3_ungapped{1710889};
    result_t result3_3_gapped{561};
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
    result_t ungapped{6939,1633,353,49261,25023};
    result_t gapped{51, 37, 21, 197, 67};
    EXPECT_RANGE_EQ(ungapped, text | ungapped_view);
    EXPECT_RANGE_EQ(gapped, text | gapped_view);
}

TEST_F(minstrobe_hash_test, ungapped)
{
    EXPECT_RANGE_EQ(result1, text1 | ungapped_view);
    EXPECT_RANGE_EQ(result3_ungapped, text3 | ungapped_view);
    EXPECT_NO_THROW(text1 | minstrobe2_hash(ungapped_shape,3,6));
    EXPECT_THROW((text3 | minstrobe2_hash(ungapped_shape,3,2)), std::invalid_argument);

    EXPECT_RANGE_EQ(result3_1, text1_3 | ungapped3_view);
    EXPECT_RANGE_EQ(result3_3_ungapped, text3_3 | ungapped3_view);
    //EXPECT_NO_THROW((text1 | minstrobe3_hash(ungapped_shape,3,6))); // Todo: Fix, Do I want to throw when sequence not long enough? 
    EXPECT_THROW((text3 | minstrobe3_hash(ungapped_shape,3,2)), std::invalid_argument);
}

TEST_F(minstrobe_hash_test, gapped)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_view);
    EXPECT_RANGE_EQ(result3_gapped, text3 | gapped_view);
    EXPECT_NO_THROW(text1 | minstrobe2_hash(gapped_shape, 2,5));
    EXPECT_THROW((text3 | minstrobe2_hash(gapped_shape, 2,1)), std::invalid_argument);

    EXPECT_RANGE_EQ(result3_1, text1_3 | gapped3_view);
    EXPECT_RANGE_EQ(result3_3_gapped, text3_3 | gapped3_view);
    //EXPECT_NO_THROW((text1 | minstrobe3_hash(gapped_shape, 2,5)), std::invalid_argument);
    EXPECT_THROW((text3 | minstrobe3_hash(gapped_shape, 2,1)), std::invalid_argument);
}

TEST_F(minstrobe_hash_test, combinability)
{
    auto start_at_a = std::views::drop(3);
    EXPECT_RANGE_EQ(result3_ungapped_start, text3 | start_at_a | ungapped_view);
    EXPECT_RANGE_EQ(result3_gapped_start, text3 | start_at_a | gapped_view);
}
