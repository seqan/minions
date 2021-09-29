#include <chrono>

#include "compare.h"

/*! \brief Function, comparing the methods.
 *  \param sequence_files A vector of sequence files.
 *  \param input_view View that should be tested.
 *  \param method_name Name of the tested method.
 *  \param path_out The output path.
 */
template <typename urng_t>
void compare(std::vector<std::filesystem::path> sequence_files, urng_t input_view, std::string method_name, std::filesystem::path path_out)
{
    std::vector<int> speed_results{};
    std::ofstream outfile;

    for (int i = 0; i < sequence_files.size(); ++i)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
        robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
        auto start = std::chrono::high_resolution_clock::now();
        for (auto & [seq] : fin)
        {
            for (auto && hash : seq | input_view)
                hash_table[hash] = std::min<uint16_t>(65534u, hash_table[hash] + 1);
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        speed_results.push_back(duration.count());

        // Store representative k-mers
        outfile.open(std::string{path_out} + method_name + "_"+ std::string{sequence_files[i].stem()} + ".out", std::ios::binary);
        for (auto && hash : hash_table)
        {
            outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
            outfile.write(reinterpret_cast<const char*>(&hash.second), sizeof(hash.second));
        }
        outfile.close();
    }

    double sum = std::accumulate(speed_results.begin(), speed_results.end(), 0.0);
    double mean = sum / speed_results.size();
    std::vector<double> diff(speed_results.size());
    std::transform(speed_results.begin(), speed_results.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / speed_results.size());
    std::cout << method_name << "\t" << *std::min_element(speed_results.begin(), speed_results.end()) << "\t" << mean << "\t" << stdev << "\t" << *std::max_element(speed_results.begin(), speed_results.end()) << "\n";
}

void do_comparison(std::vector<std::filesystem::path> sequence_files, range_arguments & args, std::filesystem::path path_out)
{
    if (args.kmers)
        compare(sequence_files, seqan3::views::kmer_hash(seqan3::ungapped{args.k_size}), "kmer_hash_"+std::to_string(args.k_size), path_out);
    else if (args.minimiser)
        compare(sequence_files, seqan3::views::minimiser_hash(args.shape,
                                                          args.w_size), "minimiser_hash_" + std::to_string(args.k_size) + "_" + std::to_string(args.w_size.get()), path_out);
}
