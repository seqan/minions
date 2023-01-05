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

#include "syncmer_hash.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

using iterator_type = std::ranges::iterator_t<decltype(std::declval<seqan3::dna4_vector&>()
                                                       | syncmer_hash(2, 5, {0}, seqan3::seed{0}))>;

static auto open_view = syncmer_hash(2, 5, {0}, seqan3::seed{0});
static auto open_view_t1 = syncmer_hash(2, 5, {1}, seqan3::seed{0});
static auto open_view_t2 = syncmer_hash(2, 5, {2}, seqan3::seed{0});
static auto closed_view = syncmer_hash(2, 5, {0,3}, seqan3::seed{0});
static auto open_view_14 = syncmer_hash(1, 4, {0}, seqan3::seed{0});
static auto open_view_14_t1 = syncmer_hash(1, 4, {1}, seqan3::seed{0});
static auto open_view_14_t2 = syncmer_hash(1, 4, {2}, seqan3::seed{0});
static auto closed_view_14 = syncmer_hash(1, 4, {0,3}, seqan3::seed{0});

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = false;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    result_t expected_range{105,109,27,6,764};

    using test_range_t = decltype(text | open_view);
    test_range_t test_range = text | open_view;
};

using test_type = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

template <typename T>
class open_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;

TYPED_TEST_SUITE(open_view_properties_test, underlying_range_types, );

class syncmer_hash_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAA"_dna4};
    result_t result1_open{0, 0};

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
//             Kmers(2,5):    ACGGC CGGCG GGCGA GCGAC CGACG GACGT ACGTT CGTTT GTTTA TTTAG
//   Canonical Kmers(2,5):    ACGGC cgccg GGCGA GCGAC CGACG acgtc aacgt aaacg GTTTA ctaaa
//            Hashed(2,5):    105,  406,  664,  609,  390,  109,    27,     6,  764,  448
//          Syncmers(2,5):    ACGGC                         acgtc aacgt aaacg GTTTA
//      syncmer stop(2,5):    ACGGC
//     syncmer start(2,5):                                        aacgt aaacg GTTTA ctaaa
//        c-Syncmers(2,5):    ACGGC CGGCG       GCGAC             ACGTT CGTTT GTTTA TTTAG
//    c-syncmer stop(2,5):    ACGGC CGGCG       GCGAC
//   c-syncmer start(2,5):                                        ACGTT CGTTT GTTTA TTTAG
    result_t result3_open{105,109,27,6,764};
    result_t result3_open_t1{};
    result_t result3_open_t2{406,664,390,448};
    result_t result3_stop_open{105};
    result_t result3_start_open{27,6,764};
    result_t result3_closed{105,609,109,27,6,764};
    result_t result3_stop_closed{105,609};
    result_t result3_start_closed{27,6,764};
//             Kmers(1,4):    ACGG  CGGC  GGCG  GCGA  CGAC  GACG  ACGT  CGTT  GTTT  TTTA  TTAG
//   Canonical Kmers(1,4):    ACGG  CGGC  cgcc  GCGA  CGAC  cgtc  ACGT  aacg  aaac  taaa  ctaa
//            Hashed(1,4):    26,   105,  101,  152,  97,   109,  27,     6,     1,  192,  112
//          Syncmers(1,4):    ACGG  CGGC  cgcc              cgtc  ACGT  aacg  aaac
//     Syncmers stop(1,4):    ACGG  CGGC
//    Syncmers start(1,4):                                        ACGT  CGTT  GTTT
//        c-Syncmers(1,4):    ACGG  CGGC        GCGA              ACGT  CGTT  GTTT TTTA
//   c-Syncmers stop(1,4):    ACGG  CGGC        GCGA
//  c-Syncmers start(1,4):                                        ACGT  CGTT  GTTT TTTA
    result_t result3_14_open{26,105,101,109,27,6,1};
    result_t result3_14_open_t1{192};
    result_t result3_14_open_t2{97,112};
    result_t result3_14_stop_open{26,105,101,109};
    result_t result3_14_start_open{27,6,1};
    result_t result3_14_closed{26,105,101,152,109,27,6,1};
    result_t result3_14_stop_closed{26,105,101,152,109};
    result_t result3_14_start_closed{27,6,1};
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

TYPED_TEST(open_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t open{109,109,27,6,764};
    result_t closed{109,109,27,6,764};
    EXPECT_RANGE_EQ(open, text | open_view);
    EXPECT_RANGE_EQ(closed, text | closed_view);
}

TEST_F(syncmer_hash_test, open)
{
    EXPECT_RANGE_EQ(result1_open, text1 | open_view);
    EXPECT_RANGE_EQ(result3_open, text3 | open_view);
    EXPECT_RANGE_EQ(result3_open_t1, text3 | open_view_t1);
    EXPECT_RANGE_EQ(result3_open_t2, text3 | open_view_t2);
    EXPECT_RANGE_EQ(result3_14_open, text3 | open_view_14);
    EXPECT_RANGE_EQ(result3_14_open_t1, text3 | open_view_14_t1);
    EXPECT_RANGE_EQ(result3_14_open_t2, text3 | open_view_14_t2);
    EXPECT_THROW((text3 | syncmer_hash(6, 5, {0}, seqan3::seed{0})), std::invalid_argument);

}

TEST_F(syncmer_hash_test, closed)
{
    EXPECT_RANGE_EQ(result1_open, text1 | closed_view);
    EXPECT_RANGE_EQ(result3_closed, text3 | closed_view);
    EXPECT_RANGE_EQ(result3_14_closed, text3 | closed_view_14);
    EXPECT_THROW((text3 | syncmer_hash(6, 5, {0,1}, seqan3::seed{0})), std::invalid_argument);
}

TEST_F(syncmer_hash_test, combinability_open)
{
    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_stop_open, text3 | stop_at_t | open_view);
    EXPECT_RANGE_EQ(result3_14_stop_open, text3 | stop_at_t | open_view_14);

    auto start_at_a = std::views::drop(6);
    EXPECT_RANGE_EQ(result3_start_open, text3 | start_at_a | open_view);
    EXPECT_RANGE_EQ(result3_14_start_open, text3 | start_at_a | open_view_14);
}

TEST_F(syncmer_hash_test, combinability_closed)
{
    auto stop_at_t = std::views::take_while([] (seqan3::dna4 const x) { return x != 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_stop_closed, text3 | stop_at_t | closed_view);
    EXPECT_RANGE_EQ(result3_14_stop_closed, text3 | stop_at_t | closed_view_14);

    auto start_at_a = std::views::drop(6);
    EXPECT_RANGE_EQ(result3_start_closed, text3 | start_at_a | closed_view);
    EXPECT_RANGE_EQ(result3_14_start_closed, text3 | start_at_a | closed_view_14);
}
