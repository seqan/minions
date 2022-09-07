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

    //!\brief The lower offset for the position of the next window.
    size_t window_min{};

    //!\brief The number of elements in a window.
    size_t window_size{};

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
    * \param[in] window_min  The lower offset for the position of the next window from the previous one.
    * \param[in] window_size The number of elements in a window.
    */
    hybridstrobe_view(urng_t urange, size_t const window_min, size_t const window_size) :
        urange{std::move(urange)},
        window_min{window_min},
        window_size{window_size}
    {}

    /*!\brief Construct from a non-view that can be view-wrapped and the two (lower and upper) offsets
    *        of the second window.
    * \tparam other_urng_t   The type of another urange. Must model std::ranges::viewable_range and be
                             constructible from urng_t.
    * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and
    *                        std::ranges::forward_range.
    * \param[in] window_min  The lower offset for the position of the next window from the previous one.
    * \param[in] window_size The number of elements in a window.
    */
    template <typename other_urng_t>
    //!\cond
        requires (std::ranges::viewable_range<other_urng_t> &&
                  std::constructible_from<urng_t, ranges::ref_view<std::remove_reference_t<other_urng_t>>>)
    //!\endcond
    hybridstrobe_view(other_urng_t && urange, size_t const window_min, size_t const window_size) :
        urange{std::views::all(std::forward<other_urng_t>(urange))},
        window_min{window_min},
        window_size{window_size}
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
                std::ranges::begin(urange),
                std::ranges::end(urange),
                window_min,
                window_size};
    }

    //!\copydoc begin()
    basic_iterator<true> begin() const
    //!\cond
        requires const_iterable
    //!\endcond
    {
        return {std::ranges::cbegin(urange),
                std::ranges::cbegin(urange),
                std::ranges::cend(urange),
                window_min,
                window_size};
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * \details
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
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
    //!\}
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
    using difference_type = std::ranges::range_difference_t<urng_t>;
    //!\brief Value type of the iterator.
    using value_t = std::ranges::range_value_t<urng_t>;
    //!\brief Value type of the output.
    using value_type = std::vector<value_t>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Reference to `value_type`.
    using reference = value_type;
    //!\brief Tag this class as a forward iterator.
    using iterator_category = std::forward_iterator_tag;
    //!\brief Tag this class as a forward iterator.
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
          urng_sentinel{std::move(it.urng_sentinel)}

    {}

    /*!\brief Construct from two begin and one end iterators of a given range over std::totally_ordered values, and the two
    *         (lower and upper) offsets of the second window.
    * \param[in] second_iterator Iterator pointing to the first position of the second window of the std::totally_ordered range.
    * \param[in] third_iterator  Iterator pointing to the first position of the third window of the std::totally_ordered range.
    * \param[in] urng_sentinel   Iterator pointing to the last position of the second window of the std::totally_ordered range.
    * \param[in] window_min      The lower offset for the position of the next window from the previous one.
    * \param[in] window_size     The number of elements in a window.
    *
    * \details
    *
    * Looks at the number of values per two windows with three iterators. First iterator adds the next value in the vector as
    * the first strobe. The second iterator adds the minimum value of the second window to the second position of the vector.
    * The third iterator adds the minimum value of the third window to the third position of the vector.
    *
    */
    basic_iterator(urng_iterator_t second_iterator,
                   urng_iterator_t third_iterator,
                   urng_sentinel_t urng_sentinel,
                   size_t window_min,
                   size_t window_size) :
        first_iterator{std::move(first_iterator)},
        second_iterator{std::move(second_iterator)},
        third_iterator{std::move(third_iterator)},
        urng_sentinel{std::move(urng_sentinel)}
    {
        size_t size{};
        if constexpr(order_3)
        {
            size = std::ranges::distance(third_iterator, urng_sentinel);
        }
        else
        {
            size = std::ranges::distance(second_iterator, urng_sentinel);
        }

        if (window_size + 1 > size)
            throw std::invalid_argument{"The given sequence is too short to satisfy the given parameters.\n"
                                        "Please choose a smaller window min and size."};
        // Throws, if the second group is not as big as the first group
        if (std::ceil((window_size - std::ceil(window_size/3.0))/2.0) != std::ceil(window_size/3.0))
            throw std::invalid_argument{"The given window size is too short.\n"
                                        "Please choose a bigger window size."};
        window_first(window_min, window_size);
    }
    //!\}

    //!\anchor basic_iterator_comparison_hybridstrobe
    //!\name Comparison operators
    //!\{

    //!\brief Compare to another basic_iterator.
    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs)
    {
        return (lhs.first_iterator == rhs.first_iterator) &&
               (lhs.second_iterator == rhs.second_iterator) &&
               (lhs.third_iterator == rhs.third_iterator);
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

    //!\brief Return the hybridstrobe.
    value_type operator*() const noexcept
    {
        return hybridstrobe_value;
    }

private:
    //!\brief The hybridstrobe value.
    value_type hybridstrobe_value{};

    //!\brief The offset relative to the beginning of the window where the hybridstrobe value is found.
    size_t hybridstrobe_position_offset{};

    //!\brief The offset relative to the beginning of the window where the hybridstrobe value is found.
    size_t hybridstrobe_position_offset3{};

    //!\brief Iterator to the first strobe of hybridstrobe.
    urng_iterator_t first_iterator{};

    //!\brief Iterator to the right most value of the window and hence the second strobe of hybridstrobe.
    urng_iterator_t second_iterator{};

    //!\brief Iterator to the right most value of the window and hence the third strobe of hybridstrobe.
    urng_iterator_t third_iterator{};

    //!\brief Iterator to last element in range.
    urng_sentinel_t urng_sentinel{};

    //!\brief Stored values per window. It is necessary to store them, because a shift can remove the current hybridstrobe.
    std::deque<value_t> window_values{};

    //!\brief Stored values per window. It is necessary to store them, because a shift can remove the current hybridstrobe.
    std::deque<value_t> window_values3{};

    int r_pos{};
    int elem_r{};

    //!\brief Advances the window of the first iterator to the next position.
    void advance_windows()
    {
        ++first_iterator;
        ++second_iterator;

        if constexpr(order_3)
            ++third_iterator;
    }

    //!\brief Calculates hybridstrobes for the first window.
    void window_first(const size_t window_min, const size_t window_size)
    {
        if (window_size == 0u)
            return;

        first_iterator = second_iterator;
        r_pos = *first_iterator % 3;
        elem_r = std::ceil(window_size/3.0);
        std::ranges::advance(second_iterator, window_min);
        if constexpr(order_3)
        {
            third_iterator = second_iterator;
            std::ranges::advance(third_iterator, window_size);
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

        auto hybridstrobe_it = std::ranges::min_element(window_values.begin() + (r_pos*elem_r),std::min(window_values.begin() + ((1+r_pos)*elem_r), window_values.end()), std::less_equal<value_t>{});

        if constexpr(order_3)
        {
            window_values3.push_back(*third_iterator);
            auto hybridstrobe_it3 = std::ranges::min_element(window_values3.begin() + (r_pos*elem_r), std::min(window_values3.begin() + ((1+r_pos)*elem_r), window_values3.end()), std::less_equal<value_t>{});
            hybridstrobe_value = {*first_iterator, *hybridstrobe_it, *hybridstrobe_it3};
        }
        else
        {
            hybridstrobe_value = {*first_iterator, *hybridstrobe_it};
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
        r_pos = *first_iterator % 3;

        if (second_iterator == urng_sentinel)
            return;

        hybridstrobe_value[0]= *first_iterator;
        window_values.pop_front();
        window_values.push_back(*second_iterator);
        auto hybridstrobe_it = std::ranges::min_element(window_values.begin() + (r_pos*elem_r), std::min(window_values.begin() + ((1+r_pos)*elem_r), window_values.end()), std::less_equal<value_t>{});

        if constexpr(order_3)
        {
            window_values3.pop_front();
            window_values3.push_back(*third_iterator);
            auto hybridstrobe_it3 = std::ranges::min_element(window_values3.begin() + (r_pos*elem_r), std::min(window_values3.begin() + ((1+r_pos)*elem_r), window_values3.end()), std::less_equal<value_t>{});
            hybridstrobe_value = {*first_iterator, *hybridstrobe_it, *hybridstrobe_it3};
        }
        else
        {
            hybridstrobe_value = {*first_iterator, *hybridstrobe_it};
        }
    }
};



//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng_t>
hybridstrobe_view(rng_t &&, size_t const window_min, size_t const window_size) -> hybridstrobe_view<std::views::all_t<rng_t>>;

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng_t, std::uint16_t ord>
hybridstrobe_view(rng_t &&, size_t const window_min, size_t const window_size) -> hybridstrobe_view<std::views::all_t<rng_t>, ord>;


// ---------------------------------------------------------------------------------------------------------------------
// hybridstrobe_fn (adaptor definition)
// ---------------------------------------------------------------------------------------------------------------------

//![adaptor_def]
//!\brief hybridstrobe's range adaptor object type (non-closure).
//!\ingroup search_views
struct hybridstrobe_fn
{
    //!\brief Store the number of values in two windows and return a range adaptor closure object.
    constexpr auto operator()(const size_t window_min, const size_t window_size) const
    {
        return adaptor_from_functor{*this, window_min, window_size};
    }

    /*!\brief Call the view's constructor with three arguments: the underlying view and an integer indicating a lower
     *        offset and another integer indicating the upper offset of the second window.
     * \tparam urng_t         The type of the input range to process. Must model std::ranges::viewable_range.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and
     *                        std::ranges::forward_range.
     * \param[in] window_min  The lower offset for the position of the next window from the previous one.
     * \param[in] window_size The number of elements in a window.
     * \returns  A range of the converted values in vectors of size 2.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, size_t const window_min, size_t const window_size) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The range parameter to views::hybridstrobe cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
                      "The range parameter to views::hybridstrobe must model std::ranges::forward_range.");

        if (window_size <= window_min)
            throw std::invalid_argument{"The chosen min and max windows are not valid."
                                        "Please choose a window_size greater than window_min."};

        return hybridstrobe_view{urange, window_min, window_size};
    }
};
//![adaptor_def]

} // namespace seqan3::detail

namespace seqan3::views
{
/*!\brief Computes hybridstrobes for a range of comparable values. A hybridstrobe consists of a starting strobe
 * concatenated with n−1 consecutively concatenated minimizers.
 * \tparam urng_t         The type of the range being processed. See below for requirements. [template
 *                        parameter is omitted in pipe notation]
 * \param[in] urange      The range being processed. [parameter is omitted in pipe notation]
 * \param[in] window_min  The lower offset for the position of the next window from the previous one.
 * \param[in] window_size The number of elements in a window.
 * \returns A range of std::totally_ordered where each value is a vector of size 2. See below for the
 *          properties of the returned range.
 * \ingroup search_views
 *
 * \details
 *
 * A hybridstrobe defined by [Sahlin](https://genome.cshlp.org/content/31/11/2080.full.pdf) consists of
 * a starting strobe concatenated with n−1 consecutively concatenated minimizers in their respective windows.
 * For example for the following list of hash values `[6, 26, 41, 38, 24, 33, 6, 27, 47]` and 3 as `window_min`,
 * 4 as `window_size`, the hybridstrobe values are `[(6,24),(26,6),(41,6),(38,6)]`.
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
 * | std::ranges::range_reference_t   | std::totally_ordered               | std::totally_ordered             |
 *
 * See the views views submodule documentation for detailed descriptions of the view properties.
 */
inline constexpr auto hybridstrobe = detail::hybridstrobe_fn{};

} // namespace seqan3::views
