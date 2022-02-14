// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hossein Eizadi Moghadam <hosseinem AT fu-berlin.de>
 * \brief Provides minstrobe_hash.
 */

#pragma once

#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/zip.hpp>
#include "minstrobe.hpp"
#include "shared.hpp"

namespace seqan3::detail
{
//!\brief seqan3::views::minstrobe_hash's range adaptor object type (non-closure).
//!\ingroup search_views
struct minstrobe_hash_fn
{
    /*!\brief Store the shape and the window min and max offsets and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_min  The lower offset for the position of the next window from the previous one.
    * \param[in] window_max  The upper offset for the position of the next window from the previous one.
    * \throws std::invalid_argument if window_min is greater than window_max or smaller than 1.
    * \returns               A range of converted elements in vectors of size 2.
    */
    constexpr auto operator()(shape const & shape, uint32_t const window_min, uint32_t const window_max) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_min, window_max};
    }

    /*!\brief Store the shape, the window size and the seed and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_min  The lower offset for the position of the next window from the previous one.
    * \param[in] window_max  The upper offset for the position of the next window from the previous one.
    * \param[in] seed        The seed to use.
    * \throws std::invalid_argument if window_min is greater than window_max or smaller than 1.
    * \returns               A range of converted elements in vectors of size 2.
    */
    constexpr auto operator()(shape const & shape, uint32_t const window_min, uint32_t const window_max, seed const seed) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_min, window_max, seed};
    }

    /*!\brief Call the view's constructor with the underlying view, a seqan3::shape and a window size as argument.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and the reference type
     *                        of the range must model seqan3::semialphabet.
     * \param[in] shape       The seqan3::shape to use for hashing.
     * \param[in] window_min  The lower offset for the position of the next window from the previous one.
     * \param[in] window_max  The upper offset for the position of the next window from the previous one.
     * \param[in] seed        The seed to use.
     * \throws std::invalid_argument if window_min is greater than window_max or smaller than 1.
     * \returns               A range of converted elements in vectors of size 2.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange,
                              shape const & shape,
                              uint32_t const window_min,
                              uint32_t const window_max,
                              seed const seed = seqan3::seed{0x8F3F73B5CF1C9ADE}) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::minstrobe_hash cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::minstrobe_hash must model std::ranges::forward_range.");
        static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
            "The range parameter to views::minstrobe_hash must be over elements of seqan3::semialphabet.");

        if (window_min <= 1 || window_max < window_min)
            throw std::invalid_argument{"The chosen parameters are not valid. "
                                        "Please choose values greater than 1 and a window_max greater than window_min."};

        auto hashed_values = std::forward<urng_t>(urange) | seqan3::views::kmer_hash(shape)
                                                           | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();});

        return seqan3::detail::minstrobe_view(hashed_values, window_min, window_max);
    }
};

} // namespace seqan3::detail

/*!\name Alphabet related views
 * \{
 */

/*!\brief                    Computes minstrobes for a range with a given shape, min and max window offsets and seed.
 * \tparam urng_t            The type of the range being processed. See below for requirements. [template parameter is
 *                           omitted in pipe notation]
 * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
 * \param[in] shape          The seqan3::shape that determines how to compute the hash value.
 * \param[in] window_min  The lower offset for the position of the next window from the previous one.
 * \param[in] window_max  The upper offset for the position of the next window from the previous one.
 * \param[in] seed           The seed used to skew the hash values. Default: 0x8F3F73B5CF1C9ADE.
 * \returns                  A range of `size_t` where each value is a vector composed of two minstrobes.
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
inline constexpr auto minstrobe_hash = seqan3::detail::minstrobe_hash_fn{};

//!\}
