#pragma once

#include <iostream>
#include <numeric>
#include <seqan3/std/filesystem>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>


//!\brief Use dna4 instead of default dna5
struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

/*! \brief Function, comparing the methods.
 *  \param sequence_files A vector of sequence files.
 *  \param input_view View that should be tested.
  *  \param method_name Name of the tested method.
 */

template <typename urng_t>
void compare(std::vector<std::filesystem::path> sequence_files, urng_t input_view, std::string method_name)
{
    std::vector<int> speed_results{};

    for (int i = 0; i < sequence_files.size(); ++i)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
        std::unordered_set<uint64_t> hashes{};
        auto start = std::chrono::high_resolution_clock::now();
        for (auto & [seq] : fin)
        {
            for (auto && hash : seq | input_view)
                hashes.insert(hash);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        speed_results.push_back(duration.count());
    }

    double sum = std::accumulate(speed_results.begin(), speed_results.end(), 0.0);
    double mean = sum / speed_results.size();
    std::vector<double> diff(speed_results.size());
    std::transform(speed_results.begin(), speed_results.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / speed_results.size());
    std::cout << method_name << "\t" << *std::min_element(speed_results.begin(), speed_results.end()) << "\t" << mean << "\t" << stdev << "\t" << *std::max_element(speed_results.begin(), speed_results.end()) << "\n";
}
