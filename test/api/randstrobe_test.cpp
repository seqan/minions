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

inline static constexpr auto randstrobe_view = seqan3::views::randstrobe(2,5);  //namen aendern

using iterator_type = std::ranges::iterator_t< decltype(std::declval<seqan3::dna4_vector&>()
                                               | kmer_view
                                               | randstrobe_view)>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{{26,97},{105,27},{166,27},{152,27},{97,27},{134,111}};

    decltype(seqan3::views::randstrobe(seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 2, 5)) test_range =
    seqan3::views::randstrobe(vec, 2, 5);
};

using test_types = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_types, );

template <typename T>
class randstrobe_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
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
    //  old hash function von randstrobes: g(s,t) = (h(s) - h(t)) % 3
    //  new hash function von randstrobes: g(s,t) = h(t) / h(t) - h(s)
    //final hash function von randstrobes:  Sahlin  (h(k-1) + h(k')) % 5
    //            ungapped randstrobes: ACGGCGAC, CGGCACGT, GGCGACGT, GCGAACGT, CGACACGT, GACGCGTT
    //    old    ungapped randstrobes: ACGGGACG, CGGCACGT, GGCGACGT, GCGAGTTT, CGACGTTT, GACGTTAG
    //    new    ungapped randstrobes: ACGGGGCG, CGGCACGT, GGCGACGT, GCGAACGT, CGACACGT, GACGACGT
    //   final   ungapped randstrobes: ACGGGACG, CGGCGGCG, GGCGGACG, GCGAGACG, CGACGACG, GACGCGTT
    //              gapped randstrobes: A--GC--C, C--CA--T, G--GA--T, G--AA--T, C--CA--T, G--GC--T
    //    old      gapped randstrobes: A--GG--G, C--CA--T, G--GA--T, G--AG--T, C--CG--T, G--GT--G
    //    new      gapped randstrobes: A--GG--G, C--CA--T, G--GA--T, G--AA--T, C--CA--T, G--GA--T
    //    final    gapped randstrobes: A--GG--G, C--CG--G, G--GG--G, G--AG--G, C--CG--G, G--GC--T
    //  stop at T ungapped randstrobes: ACGGCGAC
    //                     randstrobe
    //  stop at ...       randstrobes: ACGGGACG
    //    stop at T gapped randstrobes: A--GC--C
    //
    // start at A ungapped randstrobes:                               GCGAACGT, CGACACGT, GACGCGTT  (erste 3 buchstaben geloescht)
    //   start at A gapped randstrobes:                               G--AA--T, C--CA--T, G--GC--T
   
   //old result_t result3_ungapped{{26,134},{105,27},{166,27},{152,191},{97,191},{134,242}};
   //old2 result_t result3_ungapped{{26,166},{105,27},{166,27},{152,27},{97,27},{134,27}};
    result_t result3_ungapped{{26,134},{105,166},{166,134},{152,134},{97,134},{134,111}};

   //old result_t result3_gapped{{2,10},{5,3},{10,3},{8,11},{5,11},{10,14}};
   //old2 result_t result3_gapped{{2,10},{5,3},{10,3},{8,3},{5,3},{10,3}};
    result_t result3_gapped{{2,10},{5,10},{10,10},{8,10},{5,10},{10,7}};
    result_t result3_ungapped_stop{{26,97}};  //?
    result_t result3_gapped_stop{{2,5}};      //?
    result_t result3_ungapped_start{{152,27},{97,27},{134,111}};  //?
    result_t result3_gapped_start{{8,3},{5,3},{10,7}};    //?
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

TYPED_TEST(randstrobe_view_properties_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG

    auto v = text | kmer_view | randstrobe_view;
    compare_types(v);
}

TYPED_TEST(randstrobe_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    //                k-mers:   ACGT,      CGTC,      GTCG,      TCGA,      CGAC,      GACG,      ACGT,      CGTT,      GTTT,      TTTA,      TTAG
    //       ungapped Hashes:    27,       109,       182,       216,        97,       134,        27,        111,      191,        252,      242
    //        gapped  Hashes:    3,         5,         10,        12,         5,       10,         3,          7,        11,        12,        14
    // old ungapped randstrobes:   ACGTACGT,  CGTCGACG,  GTCGGTTT,  TCGACGTT,  CGACGTTT,  GACGTTAG
    // new ungapped randstrobes:   ACGTTCGA,  CGTCACGT,  GTCGACGT,  TCGAACGT,  CGACACGT,  GACGACGT
    //finalungapped randstrobes:   ACGTCGTC,  CGTCTCGA,  GTCGGACG,  TCGAGACG,  CGACGACG,  GACGCGTT
    // old result_t ungapped{{27,27},{109,134},{182,191},{216,111},{97,191},{134,242}};  //changed according to randstrobes, was also originally till GACG
    // old2 result_t ungapped{{27,216},{109,27},{182,27},{216,27},{97,27},{134,27}};
    result_t ungapped{{27,109},{109,216},{182,134},{216,134},{97,134},{134,111}};
    // old result_t gapped{{3,3},{5,10},{10,11},{12,7},{5,11},{10,14}};                  //changed according to randstrobes
    // old2 result_t gapped{{3,12},{5,3},{10,3},{12,3},{5,3},{10,3}};
    result_t gapped{{3,5},{5,12},{10,10},{12,10},{5,10},{10,7}};
    EXPECT_RANGE_EQ(ungapped, text | kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(gapped, text | gapped_kmer_view | randstrobe_view);
}

TEST_F(randstrobe_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_ungapped, text3 | kmer_view | randstrobe_view);
}

TEST_F(randstrobe_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped, text3 | gapped_kmer_view | randstrobe_view);
}

TEST_F(randstrobe_test, combinability)
{
    auto start_at_a = std::views::drop(3);
    EXPECT_RANGE_EQ(result3_ungapped_start, text3 | start_at_a | kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped_start, text3 | start_at_a | gapped_kmer_view | randstrobe_view);

    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_stop, text3 | stop_at_t | kmer_view | randstrobe_view);
    EXPECT_RANGE_EQ(result3_gapped_stop, text3 | stop_at_t | gapped_kmer_view | randstrobe_view);
}
