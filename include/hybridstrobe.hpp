// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra darvish AT fu-berlin.de>
 * \brief Provides hybridstrobe.
 */

#pragma once

#include <seqan3/std/algorithm>
#include <deque>

#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

#include "shared.hpp"

namespace seqan3::detail
{
// ---------------------------------------------------------------------------------------------------------------------
// hybridstrobe_view class
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The type returned by hybridstrobe.
 * \tparam urng_t The type of the underlying range, must model std::ranges::forward_range, the reference type must
 *                 model std::totally_ordered. The typical use case is that the reference type is the result of
 *                 seqan3::kmer_hash.
 * \tparam order  The order of strobemers, at the moment the order 2 and 3 are supported. Default is 2.
 * \implements std::ranges::view
 * \ingroup search_views
 *
 * \note Most members of this class are generated by std::ranges::view_interface which is not yet documented here.
 */
template <std::ranges::view urng_t, std::uint16_t order = 2>
class hybridstrobe_view : public std::ranges::view_interface<hybridstrobe_view<urng_t>>
{
private:
    static_assert(std::ranges::forward_range<urng_t>, "The hybridstrobe_view only works on forward_ranges.");
    static_assert(std::totally_ordered<std::ranges::range_reference_t<urng_t>>,
                  "The reference type of the underlying range must model std::totally_ordered.");

    //!\brief Whether the given ranges are const_iterable.
    static constexpr bool const_iterable = seqan3::const_iterable_range<urng_t>;

    //!\brief Whether the given order is of range 3.
    static constexpr bool order_3 = (order == 3);

    //!\brief The underlying range.
    urng_t urange{};

    //!\brief The distance of the second strobe to the first one.
    size_t window_dist{};

    //!\brief The number of elements in a window.
    size_t window_size{};

    //!\brief The multiplicator.
    uint64_t multi{};

    template <bool const_range>
    class basic_iterator;

    //!\brief The sentinel type of the hybridstrobe_view.
    using sentinel = std::default_sentinel_t;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
     /// \cond Workaround_Doxygen
    hybridstrobe_view() requires std::default_initializable<urng_t> = default; //!< Defaulted.
    /// \endcond
    hybridstrobe_view(hybridstrobe_view const & rhs) = default; //!< Defaulted.
    hybridstrobe_view(hybridstrobe_view && rhs) = default; //!< Defaulted.
    hybridstrobe_view & operator=(hybridstrobe_view const & rhs) = default; //!< Defaulted.
    hybridstrobe_view & operator=(hybridstrobe_view && rhs) = default; //!< Defaulted.
    ~hybridstrobe_view() = default; //!< Defaulted.

    /*!\brief Construct from a view and the two (lower and upper) offsets of the second window.
    * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and
    *                        std::ranges::forward_range.
    * \param[in] window_dist  The lower offset for the position of the next window from the previous one.
    * \param[in] window_size The number of elements in a window.
    * \param[in] multi       The multiplicator.
    */
    hybridstrobe_view(urng_t urange, size_t const window_dist, size_t const window_size, uint64_t const multi) :
        urange{std::move(urange)},
        window_dist{window_dist},
        window_size{window_size},
        multi{multi}
    {}

    /*!\brief Construct from a non-view that can be view-wrapped and the two (lower and upper) offsets
    *        of the second window.
    * \tparam other_urng_t   The type of another urange. Must model std::ranges::viewable_range and be
                             constructible from urng_t.
    * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and
    *                        std::ranges::forward_range.
    * \param[in] window_dist  The lower offset for the position of the next window from the previous one.
    * \param[in] window_size The number of elements in a window.
    * \param[in] multi       The multiplicator.
    */
    template <typename other_urng_t>
    //!\cond
        requires (std::ranges::viewable_range<other_urng_t> &&
                  std::constructible_from<urng_t, ranges::ref_view<std::remove_reference_t<other_urng_t>>>)
    //!\endcond
    hybridstrobe_view(other_urng_t && urange, size_t const window_dist, size_t const window_size, uint64_t const multi) :
        urange{std::views::all(std::forward<other_urng_t>(urange))},
        window_dist{window_dist},
        window_size{window_size},
        multi{multi}
    {}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the range.
     * \returns Iterator to the first element.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * Strong exception guarantee.
     */
    basic_iterator<false> begin()
    {
        return {std::ranges::begin(urange),
                std::ranges::end(urange),
                window_dist,
                window_size,
                multi};
    }

    //!\copydoc begin()
    basic_iterator<true> begin() const
    //!\cond
        requires const_iterable
    //!\endcond
    {
        return {std::ranges::cbegin(urange),
                std::ranges::cend(urange),
                window_dist,
                window_size,
                multi};
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    sentinel end() const
    {
        return {};
    }

    /*!\brief Returns the size of the range, if the underlying range is a std::ranges::sized_range.
     * \returns Size of range.
     */
    auto size() const
    //!\cond
        requires std::ranges::sized_range<urng_t>
    //!\endcond
    {
        if constexpr(order_3)
            return std::ranges::size(urange) - (2*(window_size + window_dist - 1));
        else
            return std::ranges::size(urange) - window_size - window_dist + 1;
    }
};

//!\brief Iterator for calculating hybridstrobes.
template <std::ranges::view urng_t, std::uint16_t order>
template <bool const_range>
class hybridstrobe_view<urng_t, order>::basic_iterator
{
private:
    //!\brief The sentinel type of the underlying range.
    using urng_sentinel_t = maybe_const_sentinel_t<const_range, urng_t>;
    //!\brief The iterator type of the underlying range.
    using urng_iterator_t = maybe_const_iterator_t<const_range, urng_t>;

    template <bool>
    friend class basic_iterator;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = typename std::iter_difference_t<urng_iterator_t>; //typename std::ranges::range_difference_t<urng_t>;
    //!\brief Value type of the iterator.
    using value_type = std::ranges::range_value_t<urng_t>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Reference to `value_type`.
    using reference = value_type;
    //!\brief Tag this class as a bidirectional iterator.
    using iterator_category = std::forward_iterator_tag;
    //!\brief Tag this class as a bidirectional iterator.
    using iterator_concept = iterator_category;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    basic_iterator() = default; //!< Defaulted.
    basic_iterator(basic_iterator const &) = default; //!< Defaulted.
    basic_iterator(basic_iterator &&) = default; //!< Defaulted.
    basic_iterator & operator=(basic_iterator const &) = default; //!< Defaulted.
    basic_iterator & operator=(basic_iterator &&) = default; //!< Defaulted.
    ~basic_iterator() = default; //!< Defaulted.

    //!\brief Allow iterator on a const range to be constructible from an iterator over a non-const range.
    basic_iterator(basic_iterator<!const_range> const & it)
    //!\cond
        requires const_range
    //!\endcond
        : hybridstrobe_value{std::move(it.hybridstrobe_value)},
          first_iterator{std::move(it.first_iterator)},
          second_iterator{std::move(it.second_iterator)},
          third_iterator{std::move(it.third_iterator)},
          urng_sentinel{std::move(it.urng_sentinel)},
          window_dist{std::move(it.window_dist)},
          window_size{std::move(it.window_size)},
          multiplicator{std::move(it.multiplicator)},
          multiplicator3{std::move(it.multiplicator3)}

    {}

    /*!\brief Construct from begin iterator and end iterator of a given range over std::totally_ordered values, and the
    *         lower offset of the second window and the number of elements in a window.
    * \param[in] first_iterator  Iterator pointing to the first position of the underlying range.
    * \param[in] urng_sentinel   Iterator pointing to the last position of the underlying range.
    * \param[in] window_dist     The lower offset for the position of the next window from the previous one.
    * \param[in] window_size     The number of elements in a window.
    * \param[in] power_multi           The multiplicator.
    *
    */
    basic_iterator(urng_iterator_t first_iterator,
                   urng_sentinel_t urng_sentinel,
                   size_t window_dist,
                   size_t window_size,
                   uint64_t power_multi) :
        first_iterator{first_iterator},
        second_iterator{first_iterator},
        third_iterator{first_iterator},
        urng_sentinel{std::move(urng_sentinel)},
        window_dist{window_dist},
        window_size{window_size}
    {
        size_t size = std::ranges::distance(first_iterator, urng_sentinel);

        if (window_size + 1 > size)
            throw std::invalid_argument{"The given range is too short to satisfy the given parameters.\n"
                                        "Please choose a smaller window size or pick a longer underlying range."};
        // Throws, if the second group is not as big as the first group
        if (std::ceil((window_size - std::ceil(window_size/3.0))/2.0) != std::ceil(window_size/3.0))
            throw std::invalid_argument{"The given window size is too small.\n"
                                        "Please choose a bigger window size."};
        if (window_size == 0u)
            throw std::invalid_argument{"The given window size is too small.\n"
                                        "Please choose a bigger window size greater than 0."};

        if constexpr (order_3)
        {
            multiplicator = my_pow(4, power_multi*2);
            multiplicator3 = my_pow(4, power_multi);
        }
        else
        {
            multiplicator = my_pow(4, power_multi);
        }

        elem_r = std::ceil(window_size/3.0);
        fill_window();
        determine_value();
    }
    //!\}

    //!\anchor basic_iterator_comparison_hybridstrobe
    //!\name Comparison operators
    //!\{

    //!\brief Compare to another basic_iterator.
    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        if constexpr(order_3)
            return (lhs.first_iterator == rhs.first_iterator) &&
                   (lhs.second_iterator == rhs.second_iterator) &&
                   (lhs.third_iterator == rhs.third_iterator);
        else
            return (lhs.first_iterator == rhs.first_iterator) &&
                   (lhs.second_iterator == rhs.second_iterator);
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator!=(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to the sentinel of the hybridstrobe_view.
    friend bool operator==(basic_iterator const & lhs, sentinel const &)
    {
        if constexpr(order_3)
            return lhs.third_iterator == lhs.urng_sentinel;
        else
            return lhs.second_iterator == lhs.urng_sentinel;
    }

    //!\brief Compare to the sentinel of the hybridstrobe_view.
    friend bool operator==(sentinel const & lhs, basic_iterator const & rhs)
    {
        return rhs == lhs;
    }

    //!\brief Compare to the sentinel of the hybridstrobe_view.
    friend bool operator!=(sentinel const & lhs, basic_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to the sentinel of the hybridstrobe_view.
    friend bool operator!=(basic_iterator const & lhs, sentinel const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator<(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        if constexpr(order_3)
            return (lhs.first_iterator < rhs.first_iterator) &&
                   (lhs.second_iterator < rhs.second_iterator) &&
                   (lhs.third_iterator < rhs.third_iterator);
        else
            return (lhs.first_iterator < rhs.first_iterator) &&
                   (lhs.second_iterator < rhs.second_iterator);
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator>(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        if constexpr(order_3)
            return (lhs.first_iterator > rhs.first_iterator) &&
                   (lhs.second_iterator > rhs.second_iterator) &&
                   (lhs.third_iterator > rhs.third_iterator);
        else
            return (lhs.first_iterator > rhs.first_iterator) &&
                   (lhs.second_iterator > rhs.second_iterator);
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator<=(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        if constexpr(order_3)
            return (lhs.first_iterator <= rhs.first_iterator) &&
                   (lhs.second_iterator <= rhs.second_iterator) &&
                   (lhs.third_iterator <= rhs.third_iterator);
        else
            return (lhs.first_iterator <= rhs.first_iterator) &&
                   (lhs.second_iterator <= rhs.second_iterator);
    }

    //!\brief Compare to another basic_iterator.
    friend bool operator>=(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        if constexpr(order_3)
            return (lhs.first_iterator >= rhs.first_iterator) &&
                   (lhs.second_iterator >= rhs.second_iterator) &&
                   (lhs.third_iterator >= rhs.third_iterator);
        else
            return (lhs.first_iterator >= rhs.first_iterator) &&
                   (lhs.second_iterator >= rhs.second_iterator);
    }
    //!\}

    //!\brief Pre-increment.
    basic_iterator & operator++() noexcept
    {
        next_hybridstrobe();
        return *this;
    }

    //!\brief Post-increment.
    basic_iterator operator++(int) noexcept
    {
        basic_iterator tmp{*this};
        next_hybridstrobe();
        return tmp;
    }

    /*!\brief Return offset between remote sentinel's position and this.
     * \attention This function is only available if sentinel_t and urng_t model std::sized_sentinel_for.
     */
    friend difference_type operator-(urng_sentinel_t const & lhs, basic_iterator const & rhs) noexcept
    //!\cond
        requires std::sized_sentinel_for<urng_sentinel_t, urng_iterator_t>
    //!\endcond
    {
        if constexpr (order_3)
            return static_cast<difference_type>(lhs - rhs.third_iterator);
        else
            return static_cast<difference_type>(lhs - rhs.second_iterator);
    }

    /*!\brief Return offset this and remote sentinel's position.
     * \attention This function is only available if urng_t and sentinel_t model std::sized_sentinel_for.
     */
    friend difference_type operator-(basic_iterator const & lhs, urng_sentinel_t const & rhs) noexcept
    //!\cond
        requires std::sized_sentinel_for<urng_iterator_t, urng_sentinel_t>
    //!\endcond
    {
        if constexpr (order_3)
            return static_cast<difference_type>(lhs.third_iterator - rhs);
        else
            return static_cast<difference_type>(lhs.second_iterator - rhs);
    }

    //!\brief Return the hybridstrobe.
    value_type operator*() const noexcept
    {
        return hybridstrobe_value;
    }

private:
    //!\brief The hybridstrobe value.
    value_type hybridstrobe_value{};

    //!\brief Iterator to the first strobe of hybridstrobe.
    urng_iterator_t first_iterator{};

    //!\brief Iterator to the second strobe of hybridstrobe.
    urng_iterator_t second_iterator{};

    //!\brief Iterator to the third strobe of hybridstrobe, if order is 3.
    urng_iterator_t third_iterator{};

    //!\brief Iterator to last element in range.
    urng_sentinel_t urng_sentinel{};

    //!\brief Stored values per window. It is necessary to store them, because a shift can remove the current hybridstrobe.
    std::deque<value_type> window_values{};

    //!\brief Stored values per window for order 3.
    std::deque<value_type> window_values3{};

    //!\brief The result of the mod operation, indicates, which part of a window to consider (sub-window).
    int r_pos{};

    //!\brief Measures how big the sub-windows are.
    int elem_r{};

    //!\brief The distance between the first strobe and the second.
    size_t window_dist{};

    //!\brief The number of elements in a window.
    size_t window_size{};

    //!\brief The multiplicator.
    uint64_t multiplicator{};

    //!\brief The multiplicator for order 3.
    uint64_t multiplicator3{};

    //!\brief Advances the window of the iterators to the next position.
    void advance_windows()
    {
        ++first_iterator;
        ++second_iterator;

        if constexpr(order_3)
        {
            ++third_iterator;
        }
    }

    //!\brief Fills window.
    void fill_window()
    {
        second_iterator = first_iterator;
        std::ranges::advance(second_iterator, window_dist);

        if constexpr(order_3)
        {
            third_iterator = second_iterator;
            std::ranges::advance(third_iterator, window_size + window_dist - 1);
        }

        for (int i = 1u; i < window_size; ++i)
        {
            window_values.push_back(*second_iterator);
            ++second_iterator;

            if constexpr(order_3)
            {
                window_values3.push_back(*third_iterator);
                ++third_iterator;
            }
        }
        window_values.push_back(*second_iterator);

        if constexpr(order_3)
        {
            window_values3.push_back(*third_iterator);
        }
    }

    //!\brief Determine hybridstrobe value based on the contents of the window and the first_iterator.
    void determine_value()
    {
        r_pos = *first_iterator % 3;
        auto hybridstrobe_it = std::ranges::min_element(window_values.begin() + (r_pos*elem_r), std::ranges::next(window_values.begin(),((1+r_pos)*elem_r), window_values.end()), std::less_equal<value_type>{});
        if constexpr(order_3)
        {
            auto hybridstrobe_it3 = std::ranges::min_element(window_values3.begin() + (r_pos*elem_r), std::ranges::next(window_values3.begin(), ((1+r_pos)*elem_r), window_values3.end()), std::less_equal<value_type>{});
            hybridstrobe_value = *first_iterator*multiplicator +*hybridstrobe_it*multiplicator3 + *hybridstrobe_it3;
        }
        else
        {
            hybridstrobe_value = *first_iterator*multiplicator + *hybridstrobe_it;
        }
    }

    /*!\brief Calculates the next hybridstrobe value.
     * \details
     * For the following windows, we remove the first window value (is now not in window_values) and add the new
     * value that results from the window shifting.
     */
    void next_hybridstrobe()
    {
        advance_windows();
        if (second_iterator == urng_sentinel)
            return;
        if constexpr(order_3)
        {
            if (third_iterator == urng_sentinel)
                return;
        }

        window_values.pop_front();
        window_values.push_back(*second_iterator);

        if constexpr(order_3)
        {
            window_values3.pop_front();
            window_values3.push_back(*third_iterator);
        }

        determine_value();
    }
};

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng_t>
hybridstrobe_view(rng_t &&, size_t const window_dist, size_t const window_size, uint64_t const multi) -> hybridstrobe_view<std::views::all_t<rng_t>>;

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng_t, std::uint16_t ord>
hybridstrobe_view(rng_t &&, size_t const window_dist, size_t const window_size, uint64_t const multi) -> hybridstrobe_view<std::views::all_t<rng_t>, ord>;

// ---------------------------------------------------------------------------------------------------------------------
// hybridstrobe_fn (adaptor definition)
// ---------------------------------------------------------------------------------------------------------------------

//![adaptor_def]
//!\brief hybridstrobe's range adaptor object type (non-closure).
//!\ingroup search_views
struct hybridstrobe_fn
{
    //!\brief Store the number of values in two windows and return a range adaptor closure object.
    constexpr auto operator()(bool order3, const size_t window_dist, const size_t window_size, uint64_t const multi) const
    {
        return adaptor_from_functor{*this, window_dist, window_size, multi, order3};
    }

    //!\brief Store the number of values in two windows and return a range adaptor closure object.
    constexpr auto operator()(const size_t window_dist, const size_t window_size, uint64_t const multi) const
    {
        return adaptor_from_functor{*this, window_dist, window_size, multi};
    }

    /*!\brief Call the view's constructor with three arguments: the underlying view and an integer indicating a lower
     *        offset and another integer indicating the upper offset of the second window.
     * \tparam urng_t         The type of the input range to process. Must model std::ranges::viewable_range.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and
     *                        std::ranges::forward_range.
     * \param[in] window_dist The offset for the position of the next window from the previous one.
     * \param[in] window_size The number of elements in a window.
     * \param[in] multi       The multiplicator.
     * \returns  A range of the converted values.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, size_t const window_dist, size_t const window_size, uint64_t const multi) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The range parameter to views::hybridstrobe cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
                      "The range parameter to views::hybridstrobe must model std::ranges::forward_range.");

        return hybridstrobe_view{urange, window_dist, window_size, multi};
    }

    /*!\brief Call the view's constructor with three arguments: the underlying view and an integer indicating a lower
     *        offset and another integer indicating the upper offset of the second window.
     * \tparam urng_t         The type of the input range to process. Must model std::ranges::viewable_range.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and
     *                        std::ranges::forward_range.
     * \param[in] window_dist The offset for the position of the next window from the previous one.
     * \param[in] window_size The number of elements in a window.
     * \param[in] multi       The multiplicator.
     * \param[in] order3      Use, if order 3 is wanted. TODO: The actual value does not matter. but make distinction between orders so much easier.
     * \returns  A range of the converted values.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, size_t const window_dist, size_t const window_size, uint64_t const multi, bool order3) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The range parameter to views::hybridstrobe cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
                      "The range parameter to views::hybridstrobe must model std::ranges::forward_range.");

        return hybridstrobe_view<urng_t, 3>{urange, window_dist, window_size, multi};
    }
};
//![adaptor_def]

} // namespace seqan3::detail

namespace seqan3::views
{
/*!\brief Computes hybridstrobes for a range of comparable values. A hybridstrobe consists of a starting strobe
 * concatenated with n−1 consecutively concatenated strobes.
 * \tparam urng_t         The type of the range being processed. See below for requirements. [template
 *                        parameter is omitted in pipe notation]
 * \param[in] urange      The range being processed. [parameter is omitted in pipe notation]
 * \param[in] window_dist The lower offset for the position of the next window from the previous one.
 * \param[in] window_size The number of elements in a window.
 * \param[in] multi       The multiplicator used to combine strobes. Should be the shape.count().
 * \returns A range of std::totally_ordered where each value is a vector of size 2. See below for the
 *          properties of the returned range.
 * \ingroup search_views
 *
 * \details
 *
 * A hybridstrobe defined by [Sahlin](https://genome.cshlp.org/content/31/11/2080.full.pdf) consists of
 * a starting strobe concatenated with n−1 consecutively concatenated strobes in their respective windows. The second
 * or third strobe is the minimizer in the sub-window of the window, the sub-window is picked according to the first strobe
 * by taking the hash value of the first strobe modulo 3. There are 3 sub-windows for a window. It is therefore advisable
 * that window size is divided by 3, but it is possible to have a smaller third sub-window.
 * For example for the following list of hash values `[6, 26, 41, 38, 24, 33, 6, 27, 47]` and 3 as `window_dist`,
 * 5 as `window_size`, the hybridstrobe values are `[(6,24),(26,47)]`.
 *
 * At the moment, only hybridstrobes of order (n) 2 or 3 is supported.
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       | *required*                         | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *preserved*                      |
 * | std::ranges::random_access_range |                                    | *lost*                           |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *preserved*                      |
 * | std::ranges::common_range        |                                    | *preserved, if random access*    |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | std::totally_ordered               | std::totally_ordered             |
 *
 * See the views views submodule documentation for detailed descriptions of the view properties.
 */
inline constexpr auto hybridstrobe = detail::hybridstrobe_fn{};

} // namespace seqan3::views
