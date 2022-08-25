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

#include "minstrobe.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<std::vector<size_t>>;

inline static constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{4});
inline static constexpr auto gapped_kmer_view = seqan3::views::kmer_hash(0b1001_shape);

inline static constexpr auto minstrobe_view = seqan3::views::minstrobe(2,4);

using iterator_type = std::ranges::iterator_t< decltype(std::declval<seqan3::dna4_vector&>()
                                               | kmer_view
                                               | minstrobe_view)>;

using order3_iterator_type = std::ranges::iterator_t<decltype(seqan3::detail::minstrobe_view<
                             decltype(std::declval<seqan3::dna4_vector &>() | kmer_view),3>
                             {std::declval<seqan3::dna4_vector &>() | kmer_view, 1, 3})>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{{26,97},{105,27},{166,27},{152,27},{97,27},{134,111}};

    decltype(seqan3::views::minstrobe(seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 2, 4)) test_range =
    seqan3::views::minstrobe(vec, 2, 4);
};

template <>
struct iterator_fixture<order3_iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{{26,105,27},{105,97,27},{166,97,27},{152,27,111},{97,27,191}};

    decltype(seqan3::detail::minstrobe_view<decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})),3>
    (seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 1, 3)) test_range =
    seqan3::detail::minstrobe_view<decltype(vec),3>(vec, 1, 3);
};

using test_types = ::testing::Types<iterator_type,order3_iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_types, );

template <typename T>
class minstrobe_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;
TYPED_TEST_SUITE(minstrobe_view_properties_test, underlying_range_types, );

class minstrobe_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAA"_dna4};
    result_t result1{{0,0},{0,0},{0,0},{0,0}}; // Same result for ungapped and gapped

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    //            ungapped minstrobes: ACGGCGAC, CGGCACGT, GGCGACGT, GCGAACGT, CGACACGT, GACGCGTT
    //              gapped minstrobes: A--GC--C, C--CA--T, G--GA--T, G--AA--T, C--CA--T, G--GC--T
    //  stop at T ungapped minstrobes: ACGGCGAC
    //    stop at T gapped minstrobes: A--GC--C
    // start at A ungapped minstrobes:                               GCGAACGT, CGACACGT, GACGCGTT
    //   start at A gapped minstrobes:                               G--AA--T, C--CA--T, G--GC--T
    result_t result3_ungapped{{26,97},{105,27},{166,27},{152,27},{97,27},{134,111}};
    result_t result3_gapped{{2,5},{5,3},{10,3},{8,3},{5,3},{10,7}};
    result_t result3_ungapped_stop{{26,97}};
    result_t result3_gapped_stop{{2,5}};
    result_t result3_ungapped_start{{152,27},{97,27},{134,111}};
    result_t result3_gapped_start{{8,3},{5,3},{10,7}};

    result_t result3_1{{0,0,0},{0,0,0},{0,0,0}}; // Same result for ungapped and gapped

    // ACGGCGACGTTTAG
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    //            ungapped minstrobes: ACGGCGGCACGT, CGGCCGACACGT, GGCGCGACACGT, GCGAACGTCGTT, CGACACGTCGTT
    //              gapped minstrobes: A--GC--CA--T, C--CC--CA--T, G--GC--CA--T, G--AA--TC--T, C--CA--TG--T
    // start at A ungapped minstrobes:                               GGCGCGACACGT, GCGAACGTCGTT, CGACACGTCGTT
    //   start at A gapped minstrobes:                               G--GC--CA--T, G--AA--TC--T, C--CA--TG--T
    result_t result3_3_ungapped{{26,105,27},{105,97,27},{166,97,27},{152,27,111},{97,27,191}};
    result_t result3_3_gapped{{2,5,3},{5,5,3},{10,5,3},{8,3,7},{5,3,11}};
    result_t result3_3_ungapped_start{{152,27,111},{97,27,191}};
    result_t result3_3_gapped_start{{8,3,7},{5,3,11}};
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

TYPED_TEST(minstrobe_view_properties_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG

    auto v = text | kmer_view | minstrobe_view;
    compare_types(v);

    auto v2 = seqan3::detail::minstrobe_view<decltype(text | kmer_view),3>(text | kmer_view,1,3);
    compare_types(v2);
}

TYPED_TEST(minstrobe_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{{27,97},{109,27},{182,27},{216,27},{97,27},{134,111}};
    result_t gapped{{3,5},{5,3},{10,3},{12,3},{5,3},{10,7}};
    EXPECT_RANGE_EQ(ungapped, text | kmer_view | minstrobe_view);
    EXPECT_RANGE_EQ(gapped, text | gapped_kmer_view | minstrobe_view);

    result_t ungapped3{{27,109,27},{109,97,27},{182,97,27},{216,27,111},{97,27,191}};
    result_t gapped3{{3,5,3},{5,5,3},{10,5,3},{12,3,7},{5,3,11}};
    EXPECT_RANGE_EQ(ungapped3, (seqan3::detail::minstrobe_view<decltype(text | kmer_view),3>(text | kmer_view,1,3)));
    EXPECT_RANGE_EQ(gapped3, (seqan3::detail::minstrobe_view<decltype(text | gapped_kmer_view),3>(text | gapped_kmer_view,1,3)));
}

TEST_F(minstrobe_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | minstrobe_view);
    EXPECT_RANGE_EQ(result3_ungapped, text3 | kmer_view | minstrobe_view);

    EXPECT_RANGE_EQ(result3_1, (seqan3::detail::minstrobe_view<decltype(text1 | kmer_view),3>(text1 | kmer_view,1,3)));
    EXPECT_RANGE_EQ(result3_3_ungapped, (seqan3::detail::minstrobe_view<decltype(text3 | kmer_view),3>(text3 | kmer_view,1,3)));
}

TEST_F(minstrobe_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_kmer_view | minstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped, text3 | gapped_kmer_view | minstrobe_view);

    EXPECT_RANGE_EQ(result3_1, (seqan3::detail::minstrobe_view<decltype(text1 | gapped_kmer_view),3>(text1 | gapped_kmer_view,1,3)));
    EXPECT_RANGE_EQ(result3_3_gapped, (seqan3::detail::minstrobe_view<decltype(text3 | gapped_kmer_view),3>(text3 | gapped_kmer_view,1,3)));
}

TEST_F(minstrobe_test, combinability)
{
    auto start_at_a = std::views::drop(3);
    EXPECT_RANGE_EQ(result3_ungapped_start, text3 | start_at_a | kmer_view | minstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped_start, text3 | start_at_a | gapped_kmer_view | minstrobe_view);

    EXPECT_RANGE_EQ(result3_3_ungapped_start, (seqan3::detail::minstrobe_view<decltype(text3 | start_at_a | kmer_view),3>(text3 | start_at_a | kmer_view,1,3)));
    EXPECT_RANGE_EQ(result3_3_gapped_start, (seqan3::detail::minstrobe_view<decltype(text3 | start_at_a| gapped_kmer_view),3>(text3 | start_at_a | gapped_kmer_view,1,3)));

    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_stop, text3 | stop_at_t | kmer_view | minstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped_stop, text3 | stop_at_t | gapped_kmer_view | minstrobe_view);
}
