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

#include "hybridstrobe_hash.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<uint64_t>;

static constexpr seqan3::shape ungapped_shape = seqan3::ungapped{4};
static constexpr seqan3::shape gapped_shape = 0b1001_shape;
static constexpr auto ungapped_view = hybridstrobe2_hash(ungapped_shape, 1, 6, seqan3::seed{0});
static constexpr auto gapped_view = hybridstrobe2_hash(gapped_shape, 1, 6, seqan3::seed{0});

static constexpr auto ungapped3_view = hybridstrobe3_hash(ungapped_shape, 1, 6, seqan3::seed{0});
static constexpr auto gapped3_view = hybridstrobe3_hash(gapped_shape, 1, 6, seqan3::seed{0});

template <typename T>
class hybridstrobe_hash_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;

TYPED_TEST_SUITE(hybridstrobe_hash_view_properties_test, underlying_range_types, );

class hybridstrobe_hash_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAA"_dna4};
    result_t result1{0,0,0}; // Same result for ungapped and gapped

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //        reverse ungapped Hashes: 91,       150,      101,      217,      182,      109,      27,   6,     1,   192,  176
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    //          reverse gapped Hashes: 7,        10,        5,       13,       10,       5,        3,    2,     1,   12,    4
    //            ungapped hybridstrobes: ACGGACGT, aacgtcgc, aaactcgc, GCGATTTA, CGACTTTA
    //              gapped hybridstrobes: A--GA--T, a--gt--c, a--ct--c, G--AT--A, c--ac--c
    // start at A ungapped hybridstrobes:                               G--AT--A, c--ac--c
    result_t result3_ungapped{6683, 1753, 473, 39164, 25084};
    result_t result3_gapped{35, 45, 29, 140, 69};
    result_t result3_ungapped_start{39164, 25084};
    result_t result3_gapped_start{140, 69};

    std::vector<seqan3::dna4> text1_3{"AAAAAAAAAAAAAAAA"_dna4};
    result_t result3_1{0}; // Same result for ungapped and gapped
    std::vector<seqan3::dna4> text3_3{"ACGGCGACGTTTAGGC"_dna4};
    // ACGGCGACGTTTAGGC
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG, TAGG, AGGC
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242,   202,   41
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14,    14,     1
    //            ungapped hybridstrobes: ACGGACGTAGGC
    //              gapped hybridstrobes: A--GA--TA--C
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

TYPED_TEST(hybridstrobe_hash_view_properties_test, different_input_ranges)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{7009,1752,472,49261,25084};
    result_t gapped{53,44,28,74,69};
    EXPECT_RANGE_EQ(ungapped, text | ungapped_view);
    //EXPECT_RANGE_EQ(gapped, text | gapped_view);
}

TEST_F(hybridstrobe_hash_test, ungapped)
{
    EXPECT_RANGE_EQ(result1, text1 | ungapped_view);
    EXPECT_RANGE_EQ(result3_ungapped, text3 | ungapped_view);
    EXPECT_NO_THROW(text1 | hybridstrobe2_hash(ungapped_shape,3,6));
    EXPECT_THROW((text3 | hybridstrobe2_hash(ungapped_shape,3,2)), std::invalid_argument);

    EXPECT_RANGE_EQ(result3_1, text1_3 | ungapped3_view);
    EXPECT_RANGE_EQ(result3_3_ungapped, text3_3 | ungapped3_view);
    EXPECT_NO_THROW(text1 | hybridstrobe3_hash(ungapped_shape,3,6));
    EXPECT_THROW((text3 | hybridstrobe3_hash(ungapped_shape,3,2)), std::invalid_argument);
}

TEST_F(hybridstrobe_hash_test, gapped)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_view);
    //EXPECT_RANGE_EQ(result3_gapped, text3 | gapped_view);
    EXPECT_NO_THROW(text1 | hybridstrobe2_hash(gapped_shape, 2,5));
    EXPECT_THROW((text3 | hybridstrobe2_hash(gapped_shape, 2,1)), std::invalid_argument);

    EXPECT_RANGE_EQ(result3_1, text1_3 | gapped3_view);
    EXPECT_RANGE_EQ(result3_3_gapped, text3_3 | gapped3_view);
    EXPECT_NO_THROW(text1 | hybridstrobe3_hash(gapped_shape, 2,5));
    EXPECT_THROW((text3 | hybridstrobe3_hash(gapped_shape, 2,1)), std::invalid_argument);
}

TEST_F(hybridstrobe_hash_test, combinability)
{
    auto start_at_a = std::views::drop(3);
    EXPECT_RANGE_EQ(result3_ungapped_start, text3 | start_at_a | ungapped_view);
    EXPECT_RANGE_EQ(result3_gapped_start, text3 | start_at_a | gapped_view);
}
