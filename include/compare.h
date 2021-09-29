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


struct range_arguments
{
   // possible methods
   bool kmers;
   bool minimiser;

   uint8_t k_size;
   seqan3::seed seed_se{0x8F3F73B5CF1C9ADEULL};
   seqan3::shape shape;
   seqan3::window_size w_size;
};

//!\brief Use dna4 instead of default dna5
struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

/*! \brief Function, comparing the methods.
 *  \param sequence_files A vector of sequence files.
 *  \param args The arguments about the view to be used.
 *  \param path_out The output path.
 */
void do_comparison(std::vector<std::filesystem::path> sequence_files, range_arguments & args,
                   std::filesystem::path path_out);
