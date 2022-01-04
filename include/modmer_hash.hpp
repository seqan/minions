// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides modmer_hash.
 */

#pragma once

#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "modmer.hpp"
#include "shared.hpp"

namespace seqan3::detail
{
//!\brief seqan3::views::modmer_hash's range adaptor object type (non-closure).
//!\ingroup search_views
struct modmer_hash_fn
{
    /*!\brief Store the shape and the window size and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] mod_used    The mod value to use.
    * \throws std::invalid_argument if the size of the shape is greater than the `mod_used`.
    * \returns               A range of converted elements.
    */
    constexpr auto operator()(shape const & shape, uint32_t const mod_used) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, mod_used};
    }

    /*!\brief Store the shape, the window size and the seed and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] mod_used    The mod value to use.
    * \param[in] seed        The seed to use.
    * \throws std::invalid_argument if the size of the shape is greater than the `mod_used`.
    * \returns               A range of converted elements.
    */
    constexpr auto operator()(shape const & shape, uint32_t const mod_used, seed const seed) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, mod_used, seed};
    }

    /*!\brief Call the view's constructor with the underlying view, a seqan3::shape and a window size as argument.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and the reference type
     *                        of the range must model seqan3::semialphabet.
     * \param[in] shape       The seqan3::shape to use for hashing.
     * \param[in] mod_used    The mod value to use.
     * \param[in] seed        The seed to use.
     * \throws std::invalid_argument if the size of the shape is greater than the `mod_used`.
     * \returns               A range of converted elements.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange,
                              shape const & shape,
                              uint32_t const mod_used,
                              seed const seed = seqan3::seed{0x8F3F73B5CF1C9ADE}) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::modmer_hash cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::modmer_hash must model std::ranges::forward_range.");
        static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
            "The range parameter to views::modmer_hash must be over elements of seqan3::semialphabet.");

        if (mod_used == 1) // Would just return urange1 without any changes
            throw std::invalid_argument{"The chosen mod_used is not valid. "
                                        "Please choose a value greater than 1."};

        auto forward_strand = std::forward<urng_t>(urange) | seqan3::views::kmer_hash(shape)
                                                           | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();});

        auto reverse_strand = std::forward<urng_t>(urange) | seqan3::views::complement
                                                           | std::views::reverse
                                                           | seqan3::views::kmer_hash(shape)
                                                           | std::views::transform([seed] (uint64_t i)
                                                                                  {return i ^ seed.get();})
                                                           | std::views::reverse;
        // fnv_hash ensures actual randomness.
        auto combined_strand = seqan3::views::zip(forward_strand, reverse_strand) | std::views::transform([seed](std::tuple<uint64_t, uint64_t> i){return fnv_hash(std::get<0>(i) + std::get<1>(i), seed.get());});
        return seqan3::detail::modmer_view(combined_strand, mod_used);
    }
};

} // namespace seqan3::detail

/*!\name Alphabet related views
 * \{
 */

/*!\brief                    Computes modmers for a range with a given shape, mod_used and seed.
 * \tparam urng_t            The type of the range being processed. See below for requirements. [template parameter is
 *                           omitted in pipe notation]
 * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
 * \param[in] shape          The seqan3::shape that determines how to compute the hash value.
 * \param[in] mod_used       The mod value to use.
 * \param[in] seed           The seed used to skew the hash values. Default: 0x8F3F73B5CF1C9ADE.
 * \returns                  A range of `size_t` where each value is the modmer of the resp. window.
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
inline constexpr auto modmer_hash = seqan3::detail::modmer_hash_fn{};

//!\}
