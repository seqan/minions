// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hossein Eizadi Moghadam <hosseinem AT fu-berlin.de>
 * \brief Provides opensyncmer_hash.
 */

#pragma once

#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/zip.hpp>
#include "syncmer.hpp"
#include "shared.hpp"

namespace seqan3::detail
{
//!\brief seqan3::views::opensyncmer_hash's range adaptor object type (non-closure).
//!\ingroup search_views
struct opensyncmer_hash_fn
{
    /*!\brief Store the kmers and the smers and return a range adaptor closure object.
    * \param[in] kmers       The K-mer size to be used.
    * \param[in] smers       The S-mer size (s<k) to be used.
    * \throws std::invalid_argument if the window size is smaller than 1.
    * \returns               A range of converted elements.
    */
    constexpr auto operator()(size_t const smers, size_t const kmers) const
    {
        return seqan3::detail::adaptor_from_functor{*this, smers, kmers};
    }

    /*!\brief Store the window size, the subwindow size and the seed and return a range adaptor closure object.
    * \param[in] kmers       The K-mer size to be used.
    * \param[in] smers       The S-mer size (s<k) to be used.
    * \param[in] seed        The seed to use.
    * \throws std::invalid_argument if the window size is smaller than 1.
    * \returns               A range of converted elements.
    */
    constexpr auto operator()(size_t const smers, size_t const kmers, seed const seed) const
    {
        return seqan3::detail::adaptor_from_functor{*this, smers, kmers, seed};
    }

    /*!\brief Call the view's constructor with the underlying view, a window size and a subwindow size as argument.
     * \param[in] urange     The input range to process. Must model std::ranges::viewable_range and
     *                       the reference type of the range must model seqan3::semialphabet.
     * \param[in] kmers      The K-mer size to be used.
     * \param[in] smers      The S-mer size (s<k) to be used.
     * \param[in] seed       The seed to use.
     * \throws std::invalid_argument if the subwindow size is smaller than 1 or kmers is smaller than smers.
     * \returns              A range of converted elements.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange,
                              size_t const smers,
			      size_t const kmers,
                              seed const seed = seqan3::seed{0x8F3F73B5CF1C9ADE}) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::opensyncmer_hash cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::opensyncmer_hash must model std::ranges::forward_range.");
        static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
            "The range parameter to views::opensyncmer_hash must be over elements of seqan3::semialphabet.");

        if (smers < 1 || kmers <= smers)
            throw std::invalid_argument{"The chosen kmers and smers are not valid."
                                        "Please choose values greater than 1 and a smer size smaller than the kmer size."};

        auto kmers_hashes = std::forward<urng_t>(urange)
                                                 | seqan3::views::kmer_hash(seqan3::shape(seqan3::ungapped(kmers)))
                                                 | std::views::transform([seed] (uint64_t i)
                                                          {return i ^ seed.get();});

        auto smers_hashes = std::forward<urng_t>(urange)
                                                 | seqan3::views::kmer_hash(seqan3::shape(seqan3::ungapped(smers)))
                                                 | std::views::transform([seed] (uint64_t i)
                                                          {return i ^ seed.get();});

        return seqan3::detail::syncmer_view<decltype(smers_hashes), decltype(kmers_hashes), true>
                                                 (smers_hashes, kmers_hashes, kmers - smers + 1);
    }
};

} // namespace seqan3::detail

/*!\name Alphabet related views
 * \{
 */

/*!\brief                    Computes opensyncmers for a range with given window and subwindow sizes, and seed.
 * \tparam urng_t            The type of the range being processed. See below for requirements.
 *                           [template parameter is omitted in pipe notation]
 * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
 * \param[in] kmers          The K-mer size to be used.
 * \param[in] smers          The S-mer size (s<k) to be used.
 * \param[in] seed           The seed used to skew the hash values. Default: 0x8F3F73B5CF1C9ADE.
 * \returns                  A range of `size_t` where each value is the opensyncmer of the resp. window.
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
inline constexpr auto opensyncmer_hash = seqan3::detail::opensyncmer_hash_fn{};

//!\}
