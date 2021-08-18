#include <chrono>

#include "compare.h"

struct output{
    std::string speed{""};
    std::string compression{""};
};

void mean_stdv(std::vector<int> const & results, std::string method_name, std::string & out)
{
    double sum = std::accumulate(results.begin(), results.end(), 0.0);
    int mean = sum / results.size();
    std::vector<double> diff(results.size());
    std::transform(results.begin(), results.end(), diff.begin(), [mean](double x) { return x - mean; });
    int sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    int stdev = std::sqrt(sq_sum / results.size());

    out = method_name + "\t" + std::to_string(*std::min_element(results.begin(), results.end())) + "\t" + std::to_string(mean) + "\t" + std::to_string(stdev) + "\t" + std::to_string(*std::max_element(results.begin(), results.end())) + "\n";
}

/*! \brief Function, comparing the methods.
 *  \param sequence_files A vector of sequence files.
 *  \param input_view View that should be tested.
 *  \param method_name Name of the tested method.
 */
template <typename urng_t>
void compare(std::vector<std::filesystem::path> sequence_files, urng_t input_view, std::string method_name, output & out)
{
    std::vector<int> compression_results{};
    std::vector<int> speed_results{};

    for (int i = 0; i < sequence_files.size(); ++i)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
        robin_hood::unordered_set<uint64_t> hashes{};
        // Measure time for speed analysis
        auto start = std::chrono::high_resolution_clock::now();
        for (auto & [seq] : fin)
        {
            for (auto && hash : seq | input_view)
                hashes.insert(hash);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        speed_results.push_back(duration.count());
        compression_results.push_back(hashes.size());
    }


    mean_stdv(speed_results, method_name, out.speed);
    mean_stdv(compression_results, method_name, out.compression);

}

void do_comparison(std::vector<std::filesystem::path> sequence_files)
{
    output kmer{};
    output minimiser_19{};
    output minimiser_23{};
    compare(sequence_files, seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{19}}), "kmer_hash (19)         ", kmer);
    compare(sequence_files, seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{19}},
                                                          seqan3::window_size{19}), "minimiser_hash (19, 19)", minimiser_19);
    compare(sequence_files, seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{19}},
                                                          seqan3::window_size{23}), "minimiser_hash (19, 23)", minimiser_23);

    std::cout << "COMPRESSION            \tMinimum\tMean\tStdDev\tMaximum\n";
    std::cout << kmer.compression << minimiser_19.compression << minimiser_23.compression;
    std::cout << "SPEED                  \tMinimum\tMean\tStdDev\tMaximum\n";
    std::cout << kmer.speed << minimiser_19.speed << minimiser_23.speed;
}
