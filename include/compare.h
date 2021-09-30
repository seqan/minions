#pragma once

#include <iostream>
#include <numeric>
#include <seqan3/std/filesystem>
#include <vector>

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

enum methods {kmer = 0, minimiser, strobemer};

struct minimiser_arguments
{
    // Needed for minimisers
    seqan3::seed seed_se{0x8F3F73B5CF1C9ADEULL};
    seqan3::shape shape;
    seqan3::window_size w_size;
};

struct strobemer_arguments
{
    // Needed for strobemers
    unsigned int w_min;
    unsigned int w_max;
    bool rand;
    bool hybrid;
    bool minstrobers;
    unsigned int order;
};

struct range_arguments : minimiser_arguments, strobemer_arguments
{
   std::filesystem::path path_out{"./"};

   methods name;
   uint8_t k_size;
};

//!\brief Use dna4 instead of default dna5
struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

/*! \brief Function, comparing the methods.
 *  \param sequence_files A vector of sequence files.
 *  \param args The arguments about the view to be used.
 */
void do_comparison(std::vector<std::filesystem::path> sequence_files, range_arguments & args);
