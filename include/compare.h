#pragma once

#include <iostream>
#include <numeric>
#include <seqan3/std/filesystem>
#include <vector>

#include <robin_hood.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

inline constexpr static uint64_t adjust_seed(uint8_t const kmer_size, uint64_t const seed = 0x8F3F73B5CF1C9ADEULL) noexcept
{
    return seed >> (64u - 2u * kmer_size);
}

enum methods {kmer = 0, minimiser, modmers, strobemer, syncmer};

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
    bool lib_implementation = false; // Set to true, if implementation of strobemer library should be used.
};

struct syncmer_arguments
{
    // Needed for syncmers
    bool closed;
    unsigned int t;
};

struct range_arguments : minimiser_arguments, strobemer_arguments, syncmer_arguments
{
   std::filesystem::path path_out{"./"};

   methods name;
   uint8_t k_size;
};

struct accuracy_arguments : range_arguments
{
   std::vector<std::filesystem::path> input_file{};
   uint64_t ibfsize{};
   size_t number_hashes{1};
   std::filesystem::path search_file{};
   std::filesystem::path solution_file{};
   float threshold{0.5};
};

//!\brief Use dna4 instead of default dna5
struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

//!\brief Use char and strings for strobemers
struct my_traits2 : seqan3::sequence_file_input_default_traits_dna
{
    using 	sequence_legal_alphabet = char;
    using sequence_alphabet = char;
    template <typename alph>
    using sequence_container = std::string;
};

/*! \brief Function, loading compressed and uncompressed ibfs
 *  \param ibf   ibf to load
 *  \param ipath Path, where the ibf can be found.
 */
template <class IBFType>
void load_ibf(IBFType & ibf, std::filesystem::path ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);
}

/*! \brief Function, which stored compressed and uncompressed ibfs
 *  \param ibf   The IBF to store.
 *  \param opath Path, where the IBF should be stored.
 */
template <class IBFType>
void store_ibf(IBFType const & ibf,
               std::filesystem::path opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(seqan3::interleaved_bloom_filter(ibf));
}

/*! \brief Function that creates the string name of the used view.
 *  \param args The arguments about the view to be used.
 */
std::string create_name(range_arguments & args);

/*! \brief Function, comparing the methods in regard of their coverage.
 *  \param args The arguments about the view to be used.
 */
void do_accuracy(accuracy_arguments & args);

/*! \brief Function, comparing the number of submers.
 *  \param sequence_files A vector of sequence files.
 *  \param args The arguments about the view to be used.
 */
void do_counts(std::vector<std::filesystem::path> sequence_files, range_arguments & args);

/*! \brief Function, comparing the methods in regard of their distance.
 *  \param sequence_file A sequence file.
 *  \param args The arguments about the view to be used.
 */
void do_distance(std::filesystem::path sequence_file, range_arguments & args);

/*! \brief Function, counting number of matches between two sequences.
 *  \param sequence_file1 The first sequence file.
 *  \param sequence_file2 The second sequence file.
 *  \param args The arguments about the view to be used.
 */
void do_match(std::filesystem::path sequence_file1, std::filesystem::path sequence_file2, range_arguments & args);

/*! \brief Function, comparing the speed.
 *  \param sequence_files A vector of sequence files.
 *  \param args The arguments about the view to be used.
 */
void do_speed(std::vector<std::filesystem::path> sequence_files, range_arguments & args);

/*! \brief Function that calculates the uniqueness of submers in given sequence files.
 *  \param sequence_files A vector of sequence files.
 *  \param method_name The name of the method.
 *  \param args The arguments about the view to be used.
 */
void unique(std::vector<std::filesystem::path> sequence_files, std::string method_name, range_arguments & args);
