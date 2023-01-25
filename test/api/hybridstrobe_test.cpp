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

#include "../../lib/seqan3/test/unit/range/iterator_test_template.hpp"

#include "hybridstrobe.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

inline static constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{4});
inline static constexpr auto gapped_kmer_view = seqan3::views::kmer_hash(0b1001_shape);

inline static constexpr auto hybridstrobe_view = seqan3::views::hybridstrobe(1,5,4);
inline static constexpr auto hybridstrobe_view_gapped = seqan3::views::hybridstrobe(1,5,2);

using iterator_type = std::ranges::iterator_t< decltype(std::declval<seqan3::dna4_vector&>()
                                               | kmer_view
                                               | hybridstrobe_view)>;

using order3_iterator_type = std::ranges::iterator_t<decltype(seqan3::detail::hybridstrobe_view<
                             decltype(std::declval<seqan3::dna4_vector &>() | kmer_view),3>
                             {std::declval<seqan3::dna4_vector &>() | kmer_view, 1, 3, 4})>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{6790, 27032, 42523, 39103, 24943, 34546};

    decltype(seqan3::views::hybridstrobe(seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 1, 5, 4)) test_range =
    seqan3::views::hybridstrobe(vec, 1, 5, 4);
};

template <>
struct iterator_fixture<order3_iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{1742875, 6923910, 10903919, 9968636, 6364156};

    decltype(seqan3::detail::hybridstrobe_view<decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})),3>
    (seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 1, 3, 4)) test_range =
    seqan3::detail::hybridstrobe_view<decltype(vec),3>(vec, 1, 3, 4);
};

using test_types = ::testing::Types<iterator_type,order3_iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_types, );

template <typename T>
class hybridstrobe_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const, //>;//,
                                              //  seqan3::bitpacked_sequence<seqan3::dna4>,
                                              //  seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;
TYPED_TEST_SUITE(hybridstrobe_view_properties_test, underlying_range_types, );

class hybridstrobe_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAAA"_dna4};
    result_t result1{0,0,0,0,0}; // Same result for ungapped and gapped

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    //            ungapped hybridstrobes: ACGGCGAC, CGGCACGT, GGCGACGT, GCGAACGT, CGACACGT, GACGCGTT
    //              gapped hybridstrobes: A--GC--C, C--CA--T, G--GA--T, G--AA--T, C--CA--T, G--GC--T
    //  stop at T ungapped hybridstrobes: ACGGCGAC
    //    stop at T gapped hybridstrobes: A--GC--C
    // start at A ungapped hybridstrobes:                               GCGAACGT, CGACACGT, GACGCGTT
    //   start at A gapped hybridstrobes:                               G--AA--T, C--CA--T, G--GC--T
    result_t result3_ungapped{6790, 27032, 42523, 39103, 24943, 34546};
    result_t result3_gapped{42, 83, 163, 139, 92, 171};
    result_t result3_ungapped_stop{6790};
    result_t result3_gapped_stop{42};
    result_t result3_ungapped_start{39103, 24943, 34546};
    result_t result3_gapped_start{139, 92, 171};

    result_t result3_1{0,0,0,0}; // Same result for ungapped and gapped

    // ACGGCGACGTTTAG
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    //            ungapped hybridstrobes: ACGGCGGCACGT, CGGCCGACACGT, GGCGCGACACGT, GCGAACGTCGTT, CGACACGTCGTT
    //              gapped hybridstrobes: A--GC--CA--T, C--CC--CA--T, G--GC--CA--T, G--AA--TC--T, C--CA--TG--T
    // start at A ungapped hybridstrobes:                               GGCGCGACACGT, GCGAACGTCGTT, CGACACGTCGTT
    //   start at A gapped hybridstrobes:                               G--GC--CA--T, G--AA--TC--T, C--CA--TG--T
    result_t order_3_ungapped{1742875, 6923910, 10903919, 9968636, 6364156};
    result_t order_3_gapped{643, 1367, 2647, 2108, 1406};
};

template <typename adaptor_t>
void compare_types(adaptor_t v)
{
    EXPECT_TRUE(std::ranges::input_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::view<decltype(v)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v), size_t>));
}

TYPED_TEST(hybridstrobe_view_properties_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG

    auto v = text | kmer_view | hybridstrobe_view;
    compare_types(v);

    auto v2 = seqan3::detail::hybridstrobe_view<decltype(text | kmer_view),3>(text | kmer_view,1,3,4);
    compare_types(v2);
}

TYPED_TEST(hybridstrobe_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{7021, 28001, 46703, 55393, 24943, 34546};
    result_t gapped{53, 83, 163, 197, 92, 171};
    EXPECT_RANGE_EQ(ungapped, text | kmer_view | hybridstrobe_view);
    EXPECT_RANGE_EQ(gapped, text | gapped_kmer_view | hybridstrobe_view_gapped);

    result_t ungapped3{1797473, 7198747, 11962047, 14180719, 6364156};
    result_t gapped3{853, 1367, 2647, 3159, 1406};
    EXPECT_RANGE_EQ(ungapped3, (seqan3::detail::hybridstrobe_view<decltype(text | kmer_view),3>(text | kmer_view,1,3,4)));
    EXPECT_RANGE_EQ(gapped3, (seqan3::detail::hybridstrobe_view<decltype(text | gapped_kmer_view),3>(text | gapped_kmer_view,1,3,2)));
}

TEST_F(hybridstrobe_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | hybridstrobe_view);
    EXPECT_RANGE_EQ(result3_ungapped, text3 | kmer_view | hybridstrobe_view);

    EXPECT_RANGE_EQ(result3_1, (seqan3::detail::hybridstrobe_view<decltype(text1 | kmer_view),3>(text1 | kmer_view,1,3, 4)));
    EXPECT_RANGE_EQ(order_3_ungapped, (seqan3::detail::hybridstrobe_view<decltype(text3 | kmer_view),3>(text3 | kmer_view,1,3, 4)));
}

TEST_F(hybridstrobe_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_kmer_view | hybridstrobe_view_gapped);
    EXPECT_RANGE_EQ(result3_gapped, text3 | gapped_kmer_view | hybridstrobe_view_gapped);

    EXPECT_RANGE_EQ(result3_1, (seqan3::detail::hybridstrobe_view<decltype(text1 | gapped_kmer_view),3>(text1 | gapped_kmer_view,1,3,2)));
    EXPECT_RANGE_EQ(order_3_gapped, (seqan3::detail::hybridstrobe_view<decltype(text3 | gapped_kmer_view),3>(text3 | gapped_kmer_view,1,3,2)));
}

TEST_F(hybridstrobe_test, combinability)
{
    auto start_at_a = std::views::drop(3);
    EXPECT_RANGE_EQ(result3_ungapped_start, text3 | start_at_a | kmer_view | hybridstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped_start, text3 | start_at_a | gapped_kmer_view | hybridstrobe_view_gapped);

    // This test leads to a compile error, I believe because underlying range is not sized, as I am not planing to use take_while, I leave it as it is. #Todo
    /*auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_stop, text3 | stop_at_t | kmer_view | hybridstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped_stop, text3 | stop_at_t | gapped_kmer_view | hybridstrobe_view);*/
}
