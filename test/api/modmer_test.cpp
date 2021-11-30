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
inline static constexpr auto rev_kmer_view = seqan3::views::complement | std::views::reverse
                                                                       | kmer_view
                                                                       | std::views::reverse;
inline static constexpr auto gapped_kmer_view = seqan3::views::kmer_hash(0b1001_shape);
inline static constexpr auto rev_gapped_kmer_view = seqan3::views::complement | std::views::reverse
                                                                              | seqan3::views::kmer_hash(0b1001_shape)
                                                                              | std::views::reverse;
inline static constexpr auto modmer_no_rev_view = modmer(4, 2);

using iterator_type = std::ranges::iterator_t< decltype(std::declval<seqan3::dna4_vector&>()
                                               | kmer_view
                                               | modmer_no_rev_view)>;
using two_ranges_iterator_type = std::ranges::iterator_t< decltype(seqan3::detail::modmer_view{
                                                          std::declval<seqan3::dna4_vector&>()
                                                          | kmer_view,
                                                          std::declval<seqan3::dna4_vector&>()
                                                          | rev_kmer_view,
                                                          4, 2})>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{26, 166, 152, 134, 252, 242};

    decltype(modmer(seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 4, 2)) test_range =
    modmer(vec, 4, 2);
};

template <>
struct iterator_fixture<two_ranges_iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    using kmer_hash_view_t = decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4}));

    kmer_hash_view_t vec = kmer_view(text);
    result_t expected_range{26, 152, 6, 192, 112};

    using reverse_kmer_hash_view_t = decltype(rev_kmer_view(text));

    using test_range_t = decltype(seqan3::detail::modmer_view{kmer_hash_view_t{}, reverse_kmer_hash_view_t{}, 4, 2});
    test_range_t test_range = seqan3::detail::modmer_view{vec, rev_kmer_view(text), 4, 2};
};

using test_types = ::testing::Types<iterator_type, two_ranges_iterator_type>;
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

    std::vector<seqan3::dna4> too_short_text{"AC"_dna4};

    // ACGG CGGC, GGCG, GCGA, CGAC, GACG, ACGT, CGTT, GTTT, TTTA, TTAG
    // ACGG CGGC cgcc GCGA CGAC cgtc ACGT aacg aaac taaa ctaa
    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    result_t result3_ungapped{26, 152, 6, 192, 112};  // ACGG, GCGA, aacg, taaa, ctaa
    result_t result3_gapped{2, 8, 2, 12, 4};       // A--G, G--A, a--g, t--a, c--a - "-" for gap
    result_t result3_ungapped_no_rev{26, 166, 152, 134, 252, 242}; // ACGG, GGCG, GCGA, GACG, TTTA, TTAG
    result_t result3_gapped_no_rev{2, 10, 8, 10, 12, 14};      // A--G, G--G, G--A, G--G, T--A, T--G "-" for gap
    result_t result3_ungapped_stop{26, 152};       // ACGG, GCGA
    result_t result3_ungapped_no_rev_stop{26, 166, 152, 134}; // ACGG, GGCG, GCGA, GACG
    result_t result3_gapped_stop{2, 8};           // A--G, G--A
    result_t result3_gapped_no_rev_stop{2, 10, 8, 10};
    result_t result3_ungapped_start{6, 192, 112};            // For start at second A, aacg, taaa, ctaa
    result_t result3_gapped_start{2, 12, 4};       // a--g, t--a, c--a - "-" for gap
    result_t result3_ungapped_no_rev_start{242};   // For start at second A, TTAG
    result_t result3_gapped_no_rev_start{14};      // For start at second A, T--G
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

    auto v = text | kmer_view | modmer_no_rev_view;
    compare_types(v);
    auto v2 = seqan3::detail::modmer_view{text | kmer_view, text | kmer_view, 4, 2};

    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        auto v3 = seqan3::detail::modmer_view{text | kmer_view, text | rev_kmer_view, 4, 2};
        compare_types(v3);
    }
}

TYPED_TEST(modmer_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{216, 6, 192, 112};               // TCGA, aacg, taaa, ctaa
    result_t gapped{12, 2, 12, 4};                          // T--A, a--g, t--a, c--a - "-" for gap
    result_t ungapped_no_rev{182, 216, 134, 252, 242};      // GTCG, TCGA, GACG, TTTA, TTAG
    result_t gapped_no_rev{10, 12, 10, 12, 14};             // G--G, T--G, T--A, G--G, T--A, T--G "-" for gap
    EXPECT_RANGE_EQ(ungapped_no_rev, text | kmer_view | modmer_no_rev_view);
    EXPECT_RANGE_EQ(gapped_no_rev, text | gapped_kmer_view | modmer_no_rev_view);

    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        EXPECT_RANGE_EQ(ungapped, (seqan3::detail::modmer_view{text | kmer_view, text | rev_kmer_view, 4, 2})) ;
        EXPECT_RANGE_EQ(gapped, (seqan3::detail::modmer_view{text | gapped_kmer_view, text | rev_gapped_kmer_view, 4, 2}));
    }
}

TEST_F(modmer_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, (seqan3::detail::modmer_view{text1 | kmer_view, text1 | rev_kmer_view, 4, 2}));
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | modmer_no_rev_view);
    auto empty_view = seqan3::detail::modmer_view{too_short_text | kmer_view, too_short_text | rev_kmer_view, 4, 2};
    EXPECT_TRUE(std::ranges::empty(empty_view));
    auto empty_view2 = too_short_text | kmer_view | modmer_no_rev_view;
    EXPECT_TRUE(std::ranges::empty(empty_view2));
    EXPECT_RANGE_EQ(result3_ungapped, (seqan3::detail::modmer_view{text3 | kmer_view, text3 | rev_kmer_view, 4, 2}));
    EXPECT_RANGE_EQ(result3_ungapped_no_rev, text3 | kmer_view | modmer_no_rev_view);

}

TEST_F(modmer_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, (seqan3::detail::modmer_view{text1 | gapped_kmer_view,
                                                             text1 | rev_gapped_kmer_view,
                                                             4, 2}));
    EXPECT_RANGE_EQ(result1, text1 | gapped_kmer_view | modmer_no_rev_view);
    auto empty_view = seqan3::detail::modmer_view{too_short_text | gapped_kmer_view,
                                                     too_short_text | rev_gapped_kmer_view,
                                                     4, 2};
    EXPECT_TRUE(std::ranges::empty(empty_view));
    auto empty_view2 = too_short_text | gapped_kmer_view | modmer_no_rev_view;
    EXPECT_TRUE(std::ranges::empty(empty_view2));
    EXPECT_RANGE_EQ(result3_gapped, (seqan3::detail::modmer_view{text3 | gapped_kmer_view,
                                                                    text3 | rev_gapped_kmer_view,
                                                                    4, 2}));
    EXPECT_RANGE_EQ(result3_gapped_no_rev, text3 | gapped_kmer_view | modmer_no_rev_view);
}

TEST_F(modmer_test, combinability)
{
    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_no_rev_stop, text3 | stop_at_t | kmer_view | modmer_no_rev_view);
    EXPECT_RANGE_EQ(result3_gapped_no_rev_stop, text3 | stop_at_t | gapped_kmer_view | modmer_no_rev_view);

    EXPECT_RANGE_EQ(result3_ungapped_stop, (seqan3::detail::modmer_view{text3 | stop_at_t | kmer_view,
                                                                           text3 | stop_at_t | rev_kmer_view,
                                                                           4, 2}));
    EXPECT_RANGE_EQ(result3_gapped_stop, (seqan3::detail::modmer_view{text3 | stop_at_t | gapped_kmer_view,
                                                                         text3 | stop_at_t | rev_gapped_kmer_view,
                                                                         4, 2}));

    auto start_at_a = std::views::drop(6);
    EXPECT_RANGE_EQ(result3_ungapped_start, (seqan3::detail::modmer_view{text3 | start_at_a | kmer_view,
                                                                   text3 | start_at_a | rev_kmer_view,
                                                                   4, 2}));
    EXPECT_RANGE_EQ(result3_gapped_start, (seqan3::detail::modmer_view{text3 | start_at_a | gapped_kmer_view,
                                                                   text3 | start_at_a | rev_gapped_kmer_view,
                                                                   4, 2}));
}

TEST_F(modmer_test, two_ranges_unequal_size)
{
    EXPECT_THROW((seqan3::detail::modmer_view{text1 | kmer_view, text3 | rev_kmer_view, 4, 2}), std::invalid_argument);
}
