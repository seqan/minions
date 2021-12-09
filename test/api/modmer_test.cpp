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

#include "modmer.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

inline static constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{4});
inline static constexpr auto gapped_kmer_view = seqan3::views::kmer_hash(0b1001_shape);

inline static constexpr auto modmer_view = modmer(2);

using iterator_type = std::ranges::iterator_t< decltype(std::declval<seqan3::dna4_vector&>()
                                               | kmer_view
                                               | modmer_view)>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{26, 166, 152, 134, 252, 242};

    decltype(modmer(seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 2)) test_range =
    modmer(vec, 2);
};

using test_types = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_types, );

template <typename T>
class modmer_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const,
                                                std::forward_list<seqan3::dna4>,
                                                std::forward_list<seqan3::dna4> const>;
TYPED_TEST_SUITE(modmer_view_properties_test, underlying_range_types, );

class modmer_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAA"_dna4};
    result_t result1{0, 0, 0}; // Same result for ungapped and gapped
    result_t result1_distance{0, 0, 0};

    std::vector<seqan3::dna4> too_short_text{"AC"_dna4};

    // ACGG CGGC, GGCG, GCGA, CGAC, GACG, ACGT, CGTT, GTTT, TTTA, TTAG
    // CCGT GCCG  CGCC  TCGC  GTCG  CGTC  ACGT  AACG  AAAC  TAAA  CTAA
    // ACGG CGGC cgcc GCGA CGAC cgtc ACGT aacg aaac taaa ctaa
    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    result_t result3_ungapped{26, 166, 152, 134, 252, 242}; // ACGG, GGCG, GCGA, GACG, TTTA, TTAG
    result_t result3_gapped{2, 10, 8, 10, 12, 14};      // A--G, G--G, G--A, G--G, T--A, T--G "-" for gap
    result_t result3_distance{0, 1, 0, 1, 3, 0};
    result_t result3_ungapped_stop{26, 166, 152, 134}; // ACGG, GGCG, GCGA, GACG
    result_t result3_gapped_stop{2, 10, 8, 10};
    result_t result3_distance_stop{0, 1, 0, 1};
    result_t result3_ungapped_start{252, 242};        // For start at second A, TTTA, TTAG
    result_t result3_gapped_start{14};      // For start at second A, T--G
    result_t result3_distance_start{3, 0};
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

TYPED_TEST(modmer_view_properties_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG

    auto v = text | kmer_view | modmer_view;
    compare_types(v);
}

TYPED_TEST(modmer_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{182, 216, 134, 252, 242};      // GTCG, TCGA, GACG, TTTA, TTAG
    result_t gapped{10, 12, 10, 12, 14};             // G--G, T--G, T--A, G--G, T--A, T--G "-" for gap
    EXPECT_RANGE_EQ(ungapped, text | kmer_view | modmer_view);
    EXPECT_RANGE_EQ(gapped, text | gapped_kmer_view | modmer_view);
}

TEST_F(modmer_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | modmer_view);
    auto empty_view = too_short_text | kmer_view | modmer_view;
    EXPECT_TRUE(std::ranges::empty(empty_view));
    EXPECT_RANGE_EQ(result3_ungapped, text3 | kmer_view | modmer_view);

    auto v1 = text1 | kmer_view;
    EXPECT_RANGE_EQ(result1_distance, (seqan3::detail::modmer_view<decltype(v1), true>(v1, 2)));
    auto v2 = text3 | kmer_view;
    EXPECT_RANGE_EQ(result3_distance, (seqan3::detail::modmer_view<decltype(v2), true>(v2, 2)));
}

TEST_F(modmer_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_kmer_view | modmer_view);
    auto empty_view = too_short_text | gapped_kmer_view | modmer_view;
    EXPECT_TRUE(std::ranges::empty(empty_view));
    EXPECT_RANGE_EQ(result3_gapped, text3 | gapped_kmer_view | modmer_view);

    auto v1 = text1 | gapped_kmer_view;
    EXPECT_RANGE_EQ(result1_distance, (seqan3::detail::modmer_view<decltype(v1), true>(v1, 2)));
    auto v2 = text3 | gapped_kmer_view;
    EXPECT_RANGE_EQ(result3_distance, (seqan3::detail::modmer_view<decltype(v1), true>(v2, 2)));
}

TEST_F(modmer_test, combinability)
{
    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_stop, text3 | stop_at_t | kmer_view | modmer_view);
    EXPECT_RANGE_EQ(result3_gapped_stop, text3 | stop_at_t | gapped_kmer_view | modmer_view);

    auto v1 = text3 | stop_at_t | kmer_view;
    auto v2 = text3 | stop_at_t | kmer_view;
    EXPECT_RANGE_EQ(result3_distance_stop, (seqan3::detail::modmer_view<decltype(v1), true>(v1, 2)));
    EXPECT_RANGE_EQ(result3_distance_stop, (seqan3::detail::modmer_view<decltype(v2), true>(v2, 2)));

    auto start_at_a = std::views::drop(6);
    EXPECT_RANGE_EQ(result3_ungapped_start, (seqan3::detail::modmer_view{text3 | start_at_a | kmer_view, 2}));

    auto v3 = text3 | start_at_a | kmer_view;
    auto v4 = text3 | start_at_a | gapped_kmer_view;
    EXPECT_RANGE_EQ(result3_distance_start, (seqan3::detail::modmer_view<decltype(v3), true>(v3, 2)));
    EXPECT_RANGE_EQ(result3_distance_start, (seqan3::detail::modmer_view<decltype(v4), true>(v4, 2)));
}
