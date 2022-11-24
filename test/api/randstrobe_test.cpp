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

#include "randstrobe.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<std::vector<size_t>>;

inline static constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{4});
inline static constexpr auto gapped_kmer_view = seqan3::views::kmer_hash(0b1001_shape);

inline static constexpr auto randstrobe_view = seqan3::views::randstrobe(2,4);

using iterator_type = std::ranges::iterator_t< decltype(std::declval<seqan3::dna4_vector&>()
                                               | kmer_view
                                               | randstrobe_view)>;

using order3_iterator_type = std::ranges::iterator_t<decltype(seqan3::detail::randstrobe_view<
                             decltype(std::declval<seqan3::dna4_vector &>() | kmer_view),3>
                             {std::declval<seqan3::dna4_vector &>() | kmer_view, 1, 3})>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{{26,97},{105,152},{166,111},{152,191},{97,252},{134,191}};

    decltype(seqan3::views::randstrobe(seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 2, 4)) test_range =
    seqan3::views::randstrobe(vec, 2, 4);
};

template <>
struct iterator_fixture<order3_iterator_type> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{{26,166,134},{105,152,27},{166,97,27},{152,134,252},{97,27,252}};

    decltype(seqan3::detail::randstrobe_view<decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})),3>
    (seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 1, 3)) test_range =
    seqan3::detail::randstrobe_view<decltype(vec),3>(vec, 1, 3);
};

using test_types = ::testing::Types<iterator_type,order3_iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_types, );

template <typename T>
class randstrobe_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                              //  seqan3::bitpacked_sequence<seqan3::dna4>,
                                              //  seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;
TYPED_TEST_SUITE(randstrobe_view_properties_test, underlying_range_types, );

class randstrobe_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAA"_dna4};
    result_t result1{{0,0},{0,0},{0,0},{0,0}}; // Same result for ungapped and gapped

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    result_t result3_ungapped{{26,97},{105,152},{166,111},{152,191},{97,252},{134,191}};
    result_t result3_gapped{{2,5},{5,3},{10,7},{8,11},{5,12},{10,7}};
    result_t result3_ungapped_stop{{26,97}};
    result_t result3_gapped_stop{{2,5}};
    result_t result3_ungapped_start{{152,191},{97,252},{134,191}};
    result_t result3_gapped_start{{8,11},{5,12},{10,7}};

    // Reverse complement: ctaaacgtcgccgt
    //                          kmers: ctaa,     taaa,     aaac,     aacg,     acgt,     cgtc,     gtcg, tcgc, cgcc, gccg, ccgt
    //                ungapped Hashes:  112,      192,        1,        6,       27,      109,      182,  217,  101,  150,  91
    //                  gapped Hashes:    4,       12,        1,        2,        3,        5,       10,   13,    5,   10,   7
    result_t result3_rev_comp_ungapped{{112,1},{192,182},{1,27},{6,109},{27,101},{109,150}};
    result_t result3_rev_comp_gapped{{4,1},{12,5},{1,3},{2,5},{3,13},{5,13}};

    result_t result3_1{{0,0,0},{0,0,0},{0,0,0}}; // Same result for ungapped and gapped

    // ACGGCGACGTTTAG
    //                          kmers: ACGG,     CGGC,     GGCG,     GCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
    //                ungapped Hashes: 26,       105,      166,      152,      97,       134,      27,   111,  191,  252,  242
    //                  gapped Hashes: 2,        5,        10,       8,        5,        10,       3,    7,    11,   12,   14
    result_t result3_3_ungapped{{26,166,134},{105,152,27},{166,97,27},{152,134,252},{97,27,252}};
    result_t result3_3_gapped{{2,5,10},{5,5,7},{10,8,3},{8,10,7},{5,3,11}};
    result_t result3_3_ungapped_start{{152,134,252},{97,27,252}};
    result_t result3_3_gapped_start{{8,10,7},{5,3,11}};
};

template <typename adaptor_t>
void compare_types(adaptor_t v)
{
    EXPECT_TRUE(std::ranges::input_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::view<decltype(v)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v), size_t>));
}

TYPED_TEST(randstrobe_view_properties_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG

    auto v = text | kmer_view | randstrobe_view;
    compare_types(v);
    EXPECT_EQ(std::ranges::random_access_range<decltype(text)>, std::ranges::random_access_range<decltype(v)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(text)>, std::ranges::common_range<decltype(v)>);

    auto v2 = seqan3::detail::randstrobe_view<decltype(text | kmer_view),3>(text | kmer_view,1,3);
    compare_types(v2);
    EXPECT_EQ(std::ranges::random_access_range<decltype(text)>, std::ranges::random_access_range<decltype(v2)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(text)>, std::ranges::common_range<decltype(v2)>);
}

TYPED_TEST(randstrobe_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
   //                          kmers: ACGT,     CGTC,     GTCG,     TCGA,     CGAC,     GACG,     ACGT, CGTT, GTTT, TTTA, TTAG
   //                ungapped Hashes: 27,       109,      182,      216,      97,       134,      27,   111,  191,  252,  242
   //                  gapped Hashes: 3,        5,        10,       12,        5,        10,       3,    7,    11,   12,   14
    result_t ungapped{{27,97},{109,216},{182,97},{216,111},{97,252},{134,191}};
    result_t gapped{{3,5},{5,12},{10,7},{12,7},{5,12},{10,7}};
    EXPECT_RANGE_EQ(ungapped, text | kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(gapped, text | gapped_kmer_view | randstrobe_view);

    result_t ungapped3{{27,182,134},{109,216,27},{182,97,27},{216,134,252},{97,27,252}};
    result_t gapped3{{3,5,10},{5,12,3},{10,10,3},{12,5,7},{5,3,11}};
    EXPECT_RANGE_EQ(ungapped3, (seqan3::detail::randstrobe_view<decltype(text | kmer_view),3>(text | kmer_view,1,3)));
    EXPECT_RANGE_EQ(gapped3, (seqan3::detail::randstrobe_view<decltype(text | gapped_kmer_view),3>(text | gapped_kmer_view,1,3)));
}

TEST_F(randstrobe_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_ungapped, text3 | kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_rev_comp_ungapped, text3 | seqan3::views::complement  | std::views::reverse | kmer_view | randstrobe_view);

    EXPECT_RANGE_EQ(result3_1, (seqan3::detail::randstrobe_view<decltype(text1 | kmer_view),3>(text1 | kmer_view,1,3)));
    EXPECT_RANGE_EQ(result3_3_ungapped, (seqan3::detail::randstrobe_view<decltype(text3 | kmer_view),3>(text3 | kmer_view,1,3)));
}

TEST_F(randstrobe_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped, text3 | gapped_kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_rev_comp_gapped, text3 | seqan3::views::complement  | std::views::reverse | gapped_kmer_view | randstrobe_view);

    EXPECT_RANGE_EQ(result3_1, (seqan3::detail::randstrobe_view<decltype(text1 | gapped_kmer_view),3>(text1 | gapped_kmer_view,1,3)));
    EXPECT_RANGE_EQ(result3_3_gapped, (seqan3::detail::randstrobe_view<decltype(text3 | gapped_kmer_view),3>(text3 | gapped_kmer_view,1,3)));
}

TEST_F(randstrobe_test, combinability)
{
    EXPECT_RANGE_EQ(std::views::reverse(result3_rev_comp_ungapped), std::views::reverse(text3 | seqan3::views::complement  | std::views::reverse | kmer_view | randstrobe_view));
    EXPECT_RANGE_EQ(std::views::reverse(result3_rev_comp_gapped), std::views::reverse(text3 | seqan3::views::complement  | std::views::reverse | gapped_kmer_view | randstrobe_view));

    auto start_at_a = std::views::drop(3);
    EXPECT_RANGE_EQ(result3_ungapped_start, text3 | start_at_a | kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped_start, text3 | start_at_a | gapped_kmer_view | randstrobe_view);

    EXPECT_RANGE_EQ(result3_3_ungapped_start, (seqan3::detail::randstrobe_view<decltype(text3 | start_at_a | kmer_view),3>(text3 | start_at_a | kmer_view,1,3)));
    EXPECT_RANGE_EQ(result3_3_gapped_start, (seqan3::detail::randstrobe_view<decltype(text3 | start_at_a| gapped_kmer_view),3>(text3 | start_at_a | gapped_kmer_view,1,3)));

    /*auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_stop, text3 | stop_at_t | kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped_stop, text3 | stop_at_t | gapped_kmer_view | randstrobe_view);*/
}