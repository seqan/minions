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

#include "opensyncmer_hash.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

using iterator_type = std::ranges::iterator_t<decltype(std::declval<seqan3::dna4_vector&>()
                                                       | opensyncmer_hash(2, 5, seqan3::seed{0}))>;

static constexpr auto opensyncmer_hash_view = opensyncmer_hash(2, 5, seqan3::seed{0});
static constexpr auto opensyncmer_hash_view_14 = opensyncmer_hash(1, 4, seqan3::seed{0});

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = false;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    result_t expected_range{105,422,111,447,764};

    using test_range_t = decltype(text | opensyncmer_hash_view);
    test_range_t test_range = text | opensyncmer_hash_view;
};

using test_type = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

template <typename T>
class opensyncmer_hash_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;

TYPED_TEST_SUITE(opensyncmer_hash_view_properties_test, underlying_range_types, );

class opensyncmer_hash_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAA"_dna4};
    result_t result1_open{0, 0};

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
//             Kmers(2,5):    ACGGC CGGCG GGCGA GCGAC CGACG GACGT ACGTT CGTTT GTTTA TTTAG
//            Hashed(2,5):    105,  422,  664,  609,  390,  539,  111,  447,  764,  1010
//       Openyncmers(2,5): ACGGC CGGCG                         ACGTT CGTTT GTTTA
//  Opensyncmer stop(2,5): ACGGC CGGCG
// Opensyncmer start(2,5):                                     ACGTT CGTTT GTTTA
    result_t result3_open{105,422,111,447,764};
    result_t result3_stop_open{105,422};
    result_t result3_start_open{111,447,764};
//             Kmers(1,4):    ACGG  CGGC  GGCG  GCGA  CGAC  GACG  ACGT  CGTT  GTTT  TTTA  TTAG
//            Hashed(1,4):    26,   105,  166,  152,  97,   134,  27,   111,  191,  252,  242
//       Openyncmers(1,4):    ACGG  CGGC                          ACGT  CGTT  GTTT
//  Openyncmers stop(1,4):    ACGG  CGGC
// Openyncmers start(1,4):                                        ACGT  CGTT  GTTT
    result_t result3_14_open{26,105,27,111,191};
    result_t result3_14_stop_open{26,105};
    result_t result3_14_start_open{27,111,191};
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

TYPED_TEST(opensyncmer_hash_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped{109,438,111,447,764};
    EXPECT_RANGE_EQ(ungapped, text | opensyncmer_hash_view);
}

TEST_F(opensyncmer_hash_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1_open, text1 | opensyncmer_hash_view);
    EXPECT_RANGE_EQ(result3_open, text3 | opensyncmer_hash_view);
    EXPECT_RANGE_EQ(result3_14_open, text3 | opensyncmer_hash_view_14);
    EXPECT_THROW((text3 | opensyncmer_hash(6, 5, seqan3::seed{0})), std::invalid_argument);
}

TEST_F(opensyncmer_hash_test, combinability)
{
    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_stop_open, text3 | stop_at_t | opensyncmer_hash_view);
    EXPECT_RANGE_EQ(result3_14_stop_open, text3 | stop_at_t | opensyncmer_hash_view_14);

    auto start_at_a = std::views::drop(6);
    EXPECT_RANGE_EQ(result3_start_open, text3 | start_at_a | opensyncmer_hash_view);
    EXPECT_RANGE_EQ(result3_14_start_open, text3 | start_at_a | opensyncmer_hash_view_14);
}
