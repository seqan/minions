// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitradarvish AT fu-berlin.de>
 * \brief Provides hybridstrobe_hash.
 */

#pragma once

#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "hybridstrobe.hpp"
#include "shared.hpp"

namespace seqan3::detail
{
//!\brief seqan3::views::hybridstrobe2_hash's range adaptor object type (non-closure).
//!\ingroup search_views
struct hybridstrobe2_hash_fn
{
    /*!\brief Store the shape and the window min and max offsets and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_min  The lower offset for the position of the next window from the previous one.
    * \param[in] window_len  The upper offset for the position of the next window from the previous one.
    * \throws std::invalid_argument if window_min is greater than window_len or smaller than 1.
    * \returns               A range of converted elements in vectors of size 2.
    */
    constexpr auto operator()(shape const & shape, uint32_t const window_min, uint32_t const window_len) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_min, window_len};
    }

    /*!\brief Store the shape, the window size and the seed and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_min  The lower offset for the position of the next window from the previous one.
    * \param[in] window_len  The upper offset for the position of the next window from the previous one.
    * \param[in] seed        The seed to use.
    * \throws std::invalid_argument if window_min is greater than window_len or smaller than 1.
    * \returns               A range of converted elements in vectors of size 2.
    */
    constexpr auto operator()(shape const & shape, uint32_t const window_min, uint32_t const window_len, seed const seed) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_min, window_len, seed};
    }

    /*!\brief Call the view's constructor with the underlying view, a seqan3::shape and a window size as argument.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and the reference type
     *                        of the range must model seqan3::semialphabet.
     * \param[in] shape       The seqan3::shape to use for hashing.
     * \param[in] window_min  The lower offset for the position of the next window from the previous one.
     * \param[in] window_len  The length of a window.
     * \param[in] seed        The seed to use.
     * \throws std::invalid_argument if window_min is greater than window_len or smaller than 1.
     * \returns               A range of converted elements in vectors of size 2.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange,
                              shape const & shape,
                              uint32_t const window_min,
                              uint32_t const window_len,
                              seed const seed = seqan3::seed{0x8F3F73B5CF1C9ADE}) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::hybridstrobe_hash cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::hybridstrobe_hash must model std::ranges::forward_range.");
        static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
            "The range parameter to views::hybridstrobe_hash must be over elements of seqan3::semialphabet.");

        if (window_len < window_min)
            throw std::invalid_argument{"The chosen parameters are not valid. "
                                        "Please choose a window_len greater than window_min."};

        auto hashed_values = std::forward<urng_t>(urange) | seqan3::views::kmer_hash(shape)
                                                          | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();});


        auto hybridstrobes = seqan3::detail::hybridstrobe_view(hashed_values, window_min + shape.size() - 1, window_len - shape.size() + 1);
        uint64_t multiplicator = my_pow(4, shape.count());
        auto forward = std::views::transform(hybridstrobes, [multiplicator] (std::vector<uint64_t> i)
                           {return combine_strobes(multiplicator, i[0], i[1]);});

        auto rev_hashed_values = std::forward<urng_t>(urange)  | seqan3::views::complement
                                                               | std::views::reverse
                                                               | seqan3::views::kmer_hash(shape)
                                                               | std::views::transform([seed] (uint64_t i)
                                                                            {return i ^ seed.get();});


        auto rev_hybridstrobes = seqan3::detail::hybridstrobe_view(rev_hashed_values, window_min + shape.size() - 1, window_len - shape.size() + 1);
        auto reverse = std::views::transform(rev_hybridstrobes, [multiplicator] (std::vector<uint64_t> i)
                         {return combine_strobes(multiplicator, i[0], i[1]);});

        std::vector<uint64_t> rev{};
        for(auto && h : reverse)
            rev.push_back(h);
        std::reverse(rev.begin(), rev.end());

        int i{0};
        std::vector<uint64_t> result{};
        for(auto && h : forward)
        {
            result.push_back(std::min(h, rev[i]));
            ++i;
        }

        return result;
    }
};

//!\brief seqan3::views::hybridstrobe3_hash's range adaptor object type (non-closure).
//!\ingroup search_views
struct hybridstrobe3_hash_fn
{
    /*!\brief Store the shape and the window min and max offsets and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_min  The lower offset for the position of the next window from the previous one.
    * \param[in] window_len  The upper offset for the position of the next window from the previous one.
    * \throws std::invalid_argument if window_min is greater than window_len or smaller than 1.
    * \returns               A range of converted elements in vectors of size 2.
    */
    constexpr auto operator()(shape const & shape, uint32_t const window_min, uint32_t const window_len) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_min, window_len};
    }

    /*!\brief Store the shape, the window size and the seed and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_min  The lower offset for the position of the next window from the previous one.
    * \param[in] window_len  The upper offset for the position of the next window from the previous one.
    * \param[in] seed        The seed to use.
    * \throws std::invalid_argument if window_min is greater than window_len or smaller than 1.
    * \returns               A range of converted elements in vectors of size 2.
    */
    constexpr auto operator()(shape const & shape, uint32_t const window_min, uint32_t const window_len, seed const seed) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_min, window_len, seed};
    }

    /*!\brief Call the view's constructor with the underlying view, a seqan3::shape and a window size as argument.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and the reference type
     *                        of the range must model seqan3::semialphabet.
     * \param[in] shape       The seqan3::shape to use for hashing.
     * \param[in] window_min  The lower offset for the position of the next window from the previous one.
     * \param[in] window_len  The length of a window.
     * \param[in] seed        The seed to use.
     * \throws std::invalid_argument if window_min is greater than window_len or smaller than 1.
     * \returns               A range of converted elements in vectors of size 2.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange,
                              shape const & shape,
                              uint32_t const window_min,
                              uint32_t const window_len,
                              seed const seed = seqan3::seed{0x8F3F73B5CF1C9ADE}) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::hybridstrobe_hash cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::hybridstrobe_hash must model std::ranges::forward_range.");
        static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
            "The range parameter to views::hybridstrobe_hash must be over elements of seqan3::semialphabet.");

        if (window_len < window_min)
            throw std::invalid_argument{"The chosen parameters are not valid. "
                                        "Please choose values greater than 0 and a window_len greater than window_min."};

        auto hashed_values = std::forward<urng_t>(urange) | seqan3::views::kmer_hash(shape)
                                                          | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();});


        auto hybridstrobes = seqan3::detail::hybridstrobe_view<decltype(hashed_values), 3>(hashed_values, window_min + shape.size() - 1, window_len - shape.size() + 1);
        uint64_t multiplicator = my_pow(4, shape.count()*2);
        uint64_t multiplicator2 = my_pow(4, shape.count());
        auto forward =  std::views::transform(hybridstrobes, [multiplicator, multiplicator2] (std::vector<uint64_t> i)
                           {return combine_strobes(multiplicator, multiplicator2, i[0], i[1], i[2]);});


        auto rev_hashed_values = std::forward<urng_t>(urange)  | seqan3::views::complement
                                                               | std::views::reverse
                                                               | seqan3::views::kmer_hash(shape)
                                                               | std::views::transform([seed] (uint64_t i)
                                                                           {return i ^ seed.get();});


        auto rev_hybridstrobes = seqan3::detail::hybridstrobe_view<decltype(rev_hashed_values), 3>(rev_hashed_values, window_min + shape.size() - 1, window_len - shape.size() + 1);
        auto reverse = std::views::transform(rev_hybridstrobes, [multiplicator, multiplicator2] (std::vector<uint64_t> i)
                        {return combine_strobes(multiplicator, multiplicator2, i[0], i[1], i[2]);});

        std::vector<uint64_t> rev{};
        for(auto && h : reverse)
            rev.push_back(h);

        std::reverse(rev.begin(), rev.end());

        int i{0};
        std::vector<uint64_t> result{};
        for(auto && h : forward)
        {
            result.push_back(std::min(h, rev[i]));
            ++i;
        }

        return result;
    }
};

} // namespace seqan3::detail

/*!\name Alphabet related views
 * \{
 */

/*!\brief                    Computes hybridstrobes for a range with a given shape, min and max window offsets and seed.
 * \tparam urng_t            The type of the range being processed. See below for requirements. [template parameter is
 *                           omitted in pipe notation]
 * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
 * \param[in] shape          The seqan3::shape that determines how to compute the hash value.
 * \param[in] window_min     The lower offset for the position of the next window from the previous one.
 * \param[in] window_len     The upper offset for the position of the next window from the previous one.
 * \param[in] seed           The seed used to skew the hash values. Default: 0x8F3F73B5CF1C9ADE.
 * \returns                  A range of `size_t` where each value is a vector composed of two hybridstrobes.
 *                           See below for the properties of the returned range.
 * \ingroup search_views
 *
 * \attention
 * Be aware of the requirements of the seqan3::views::kmer_hash view.
 *
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       | *required*                         | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *lost*                           |
 * | std::ranges::random_access_range |                                    | *lost*                           |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *lost*                           |
 * | std::ranges::common_range        |                                    | *lost*                           |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | seqan3::semialphabet               | std::size_t                      |
 *
 * See the views views submodule documentation for detailed descriptions of the view properties.
 *
 * \hideinitializer
 *
 */
inline constexpr auto hybridstrobe2_hash = seqan3::detail::hybridstrobe2_hash_fn{};

//!\}

/*!\name Alphabet related views
 * \{
 */

/*!\brief                    Computes hybridstrobes for a range with a given shape, min and max window offsets and seed.
 * \tparam urng_t            The type of the range being processed. See below for requirements. [template parameter is
 *                           omitted in pipe notation]
 * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
 * \param[in] shape          The seqan3::shape that determines how to compute the hash value.
 * \param[in] window_min     The lower offset for the position of the next window from the previous one.
 * \param[in] window_len     The upper offset for the position of the next window from the previous one.
 * \param[in] seed           The seed used to skew the hash values. Default: 0x8F3F73B5CF1C9ADE.
 * \returns                  A range of `size_t` where each value is a vector composed of two hybridstrobes.
 *                           See below for the properties of the returned range.
 * \ingroup search_views
 *
 * \attention
 * Be aware of the requirements of the seqan3::views::kmer_hash view.
 *
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       | *required*                         | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *lost*                           |
 * | std::ranges::random_access_range |                                    | *lost*                           |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *lost*                           |
 * | std::ranges::common_range        |                                    | *lost*                           |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | seqan3::semialphabet               | std::size_t                      |
 *
 * See the views views submodule documentation for detailed descriptions of the view properties.
 *
 * \hideinitializer
 *
 */
inline constexpr auto hybridstrobe3_hash = seqan3::detail::hybridstrobe3_hash_fn{};

//!\}
