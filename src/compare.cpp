#include <chrono>
#include <ranges>

#include <index.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include "compare.h"
#include "hybridstrobe_hash.hpp"
#include "minions_minimiser_hash.hpp"
#include "minstrobe_hash.hpp"
#include "modmer_hash.hpp"
#include "randstrobe_hash.hpp"
#include "syncmer_hash.hpp"

#include <seqan3/core/debug_stream.hpp>

/*! \brief Returns expected value of given list.
 *  \param results The vector from which mean and variance should be calculated of.
 */
template<typename urng_t>
double get_expected(urng_t & results)
{
    double sum{0.0};
    for(auto & elem : results)
        sum = sum + (elem*elem);
    return sum/results.size();
}

/*! \brief Calculate mean and variance of given list.
 *  \param results The vector from which mean and variance should be calculated of.
 *  \param mean Variable to store the mean.
 *  \param stdev Variable to store the variance.
 */
template<typename urng_t>
void get_mean_and_var(urng_t & results, double & mean, double & stdev)
{
    double sum = std::accumulate(results.begin(), results.end(), 0.0);
    mean = sum / results.size();
    std::vector<double> diff(results.size());
    std::transform(results.begin(),results.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    stdev = std::sqrt(sq_sum / results.size());
}


/*! \brief Function, that decides which strobemer to use.
 *  \param seq A std::string sequence.
 *  \param args The range_arguments.
 *  \param strobes_vector The vector for the strobemers.
 */
template <int strobemers_func>
void get_strobemers(std::string & seq, range_arguments const & args,
std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> & strobes_vector)
{
    if constexpr (strobemers_func == 1)
        strobes_vector = seq_to_randstrobes2(args.order, args.k_size, args.w_min, args.w_max, seq, 0);
    else if constexpr (strobemers_func == 2)
        strobes_vector = seq_to_randstrobes3(args.order, args.k_size, args.w_min, args.w_max, seq, 0);
    else if constexpr (strobemers_func == 3)
        strobes_vector = seq_to_hybridstrobes2(args.order, args.k_size, args.w_min, args.w_max, seq, 0);
    else if constexpr (strobemers_func == 4)
        strobes_vector = seq_to_minstrobes2(args.order, args.k_size, args.w_min, args.w_max, seq, 0);
}

template <typename urng_t>
void accuracy(urng_t input_view,
              std::string method_name,
              accuracy_arguments & args)
{
    // Loading/Creating the ibf.
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf;
    if ((std::filesystem::path{args.input_file[0]}.extension() == ".ibf") & (args.input_file.size() == 1))
    {
        load_ibf(ibf, args.input_file[0]);
    }
    else if (std::filesystem::path{args.input_file[0]}.extension() == ".out")
    {
        seqan3::interleaved_bloom_filter ibf_create{seqan3::bin_count{args.input_file.size()},
                                     seqan3::bin_size{args.ibfsize},
                                     seqan3::hash_function_count{args.number_hashes}};

        uint64_t minimiser;
        uint16_t minimiser_count;
        for(size_t i = 0; i < args.input_file.size(); i++)
        {
            std::ifstream infile{args.input_file[i], std::ios::binary};
            while(infile.read((char*)&minimiser, sizeof(minimiser)))
            {
                infile.read((char*)&minimiser_count, sizeof(minimiser_count));
                ibf_create.emplace(minimiser, seqan3::bin_index{i});
            }
        }
        store_ibf(ibf_create, std::string{args.path_out} + method_name + ".ibf");
        load_ibf(ibf, std::string{args.path_out} + method_name + ".ibf");
    }
    else // Sequence files
    {
        seqan3::interleaved_bloom_filter ibf_create{seqan3::bin_count{args.input_file.size()},
                                     seqan3::bin_size{args.ibfsize},
                                     seqan3::hash_function_count{args.number_hashes}};

        for(size_t i = 0; i < args.input_file.size(); i++)
        {
            for (auto && [seq] : seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>>{args.input_file[i]})
                for (auto && value : seq | input_view)
                    ibf_create.emplace(value, seqan3::bin_index{i});
        }
        store_ibf(ibf_create, std::string{args.path_out} + method_name + ".ibf");
        load_ibf(ibf, std::string{args.path_out} + method_name + ".ibf");
    }

    // Search through the ibf with a given threshold.

    // Load search sequence file.
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id, seqan3::field::seq>> fin{args.search_file};
    std::vector<std::string> ids;
    std::vector<seqan3::dna4_vector> seqs;
    for (auto & [id, seq] : fin)
    {
        ids.push_back(id);
        seqs.push_back(seq);
    }

    // Read in "solution", in which files the sequences should be present.
    robin_hood::unordered_node_map<std::string, std::vector<uint32_t>> solutions{};
    std::ifstream infile;
    std::string line;
    infile.open(args.solution_file);
    std::string name;
    int exp;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        iss >> name;
        solutions[name] = {};
        while (iss >> exp)
            solutions[name].push_back(exp);
    }
    infile.close();

    int tp = 0, tn = 0, fp = 0, fn = 0;
    std::ofstream outfile;
    outfile.open(std::string{args.path_out} + method_name + "_" + std::string{args.search_file.stem()} + ".search_out");
    // Go over the sequences in the search file.
    for (int i = 0; i < seqs.size(); ++i)
    {
        std::vector<uint32_t> counter;
        counter.assign(ibf.bin_count(), 0);
        uint64_t length = 0;
        auto agent = ibf.membership_agent();
        for (auto && hash : seqs[i] | input_view)
        {
            std::transform (counter.begin(), counter.end(), agent.bulk_contains(hash).begin(), counter.begin(),
                            std::plus<int>());
            ++length;
        }

        outfile << ids[i] << "\t";
        for (int j = 0; j < ibf.bin_count(); ++j)
        {
            bool found = (counter[j] >= (length * args.threshold));
            bool true_positive = std::binary_search(solutions[ids[i]].begin(), solutions[ids[i]].end(), j);
            if (counter[j] >= (length * args.threshold))
                outfile << j << ",";

            if (found && true_positive)
                tp++;
            else if(found && !true_positive)
                fp++;
            else if (!found &&true_positive)
                fn++;
            else if (!found && !true_positive)
                tn++;
        }
        outfile << "\n";
    }
    outfile.close();
    // Store tp, tn, fp, fn
    std::ofstream outfile2;
    outfile2.open(std::string{args.path_out} + method_name +  "_" + std::string{args.search_file.stem()} + "_accuracy.out");
    outfile2 << method_name << "\t" << tp << "\t" << tn << "\t" << fp << "\t" << fn << "\n";
    outfile2.close();
}

/*! \brief Function, counting the number of submers.
 *  \param sequence_files A vector of sequence files.
 *  \param input_view View that should be tested.
 *  \param method_name Name of the tested method.
 *  \param args The arguments about the view to be used, needed for strobemers.
 */
template <typename urng_t>
void counts(std::vector<std::filesystem::path> & sequence_files, urng_t input_view, std::string method_name, range_arguments & args)
{
    std::vector<uint64_t> counts_results{};
    std::ofstream outfile;
    for (int i = 0; i < sequence_files.size(); ++i)
    {
        robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
        for (auto & [seq] : fin)
        {
            for (auto && hash : seq | input_view)
                hash_table[hash] = std::min<uint16_t>(65534u, hash_table[hash] + 1);
        }
        counts_results.push_back(hash_table.size());

        // Store representative k-mers
        outfile.open(std::string{args.path_out} + method_name + "_"+ std::string{sequence_files[i].stem()} + "_counts.out", std::ios::binary);
        for (auto && hash : hash_table)
        {
            outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
            outfile.write(reinterpret_cast<const char*>(&hash.second), sizeof(hash.second));
        }
        outfile.close();
    }

    double mean_counts, stdev_counts;
    get_mean_and_var(counts_results, mean_counts, stdev_counts);

    // Store speed and counts
    outfile.open(std::string{args.path_out} + method_name + "_counts.out");
    outfile << method_name << "\t" << *std::min_element(counts_results.begin(), counts_results.end()) << "\t" << mean_counts << "\t" << stdev_counts << "\t" << *std::max_element(counts_results.begin(), counts_results.end()) << "\n";
    outfile.close();
}

template <typename urng_t, typename urng_t2>
void counts_strobemer(std::vector<std::filesystem::path> & sequence_files, urng_t input_view, urng_t2 input_view2, std::string method_name, range_arguments & args)
{
    std::vector<uint64_t> counts_results{};
    std::ofstream outfile;
    for (int i = 0; i < sequence_files.size(); ++i)
    {
        robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
        for (auto & [seq] : fin)
        {
            std::vector<uint64_t> res = seq | input_view;
            for (auto && hash : res | input_view2)
            {
                hash_table[hash] = std::min<uint16_t>(65534u, hash_table[hash] + 1);
            }
        }

        counts_results.push_back(hash_table.size());

        // Store representative k-mers
        outfile.open(std::string{args.path_out} + method_name + "_"+ std::string{sequence_files[i].stem()} + "_counts.out", std::ios::binary);
        for (auto && hash : hash_table)
        {
            outfile.write(reinterpret_cast<const char*>(&hash.first), sizeof(hash.first));
            outfile.write(reinterpret_cast<const char*>(&hash.second), sizeof(hash.second));
        }
        outfile.close();
    }

    double mean_counts, stdev_counts;
    get_mean_and_var(counts_results, mean_counts, stdev_counts);

    // Store speed and counts
    outfile.open(std::string{args.path_out} + method_name + "_counts.out");
    outfile << method_name << "\t" << *std::min_element(counts_results.begin(), counts_results.end()) << "\t" << mean_counts << "\t" << stdev_counts << "\t" << *std::max_element(counts_results.begin(), counts_results.end()) << "\n";
    outfile.close();
}

template <typename urng_t, bool syncmer = false>
std::vector<uint64_t> read_seq_file(std::filesystem::path sequence_file, urng_t input_view, range_arguments & args)
{
    std::vector<uint64_t> vector{};
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_file};
    for (auto & [seq] : fin)
    {
        if constexpr (syncmer)
        {
            for (auto && hash : seq | input_view)
            {
                if (syncmer_filter(hash, args.w_size.get(), (args.k_size *args.order),  args.positions, args.seed_se.get()))
                    vector.push_back(hash);
            }
        }
        else
        {
            for (auto && hash : seq | input_view)
                vector.push_back(hash);
        }
    }

    return vector;
}

template <typename urng_t, typename urng_t2>
std::vector<uint64_t> read_seq_file(std::filesystem::path sequence_file, urng_t input_view, urng_t2 input_view2)
{
    std::vector<uint64_t> vector{};
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_file};
    for (auto & [seq] : fin)
    {
        auto v = seq | input_view;
        for (auto && hash : v | input_view2)
            vector.push_back(hash);
    }

    return vector;
}

template <int strobemers>
std::vector<uint64_t> read_seq_file(std::filesystem::path sequence_file, range_arguments & args)
{
    std::vector<uint64_t> vector{};
    seqan3::sequence_file_input<my_traits2, seqan3::fields<seqan3::field::seq>> fin{sequence_file};
    for (auto & [seq] : fin)
    {
        std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> strobes_vector;
        get_strobemers<strobemers>(seq, args, strobes_vector);
        for (auto & t : strobes_vector) // iterate over the strobemer tuples
            vector.push_back(std::get<0>(t));
    }

    return vector;
}

template <typename urng_t, typename urng_t2>
void distance(std::filesystem::path sequence_file, urng_t input_view, urng_t2 compare_view, std::string method_name)
{
    std::vector<uint64_t> distances{};
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_file};
    for (auto & [seq] : fin)
    {
        int distance = 0;
        auto representative = seq | input_view;
        auto rep_it = representative.begin();
        auto compare = seq | compare_view;
        auto comp_it = compare.begin();
        while((rep_it != representative.end()) & (comp_it != compare.end()))
        {
            if (*rep_it == *comp_it)
            {
                if (comp_it != compare.begin())
                {
                    distances.push_back(distance);
                    distance = 0;
                }
                rep_it++;
            }
            distance++;
            comp_it++;
        }
    }

    std::ofstream outfile;
    outfile.open(method_name + "_"+ std::string{sequence_file.stem()} + "_distances.out");
    for (int i = *std::min_element(distances.begin(), distances.end()); i <= *std::max_element(distances.begin(), distances.end()); i++)
    {
        auto count = std::count(distances.begin(), distances.end(), i);
        if (count > 0)
            outfile << i << "\t" << count << "\n";
    }
    outfile.close();
}

template <typename urng_t, typename urng_t2>
void distance_strobemer(std::filesystem::path sequence_file, urng_t input_view, urng_t2 compare_view, std::string method_name)
{
    std::vector<uint64_t> distances{};
    int distance = 0;
    int it_1 = 0;
    int it_2 = 0;

    std::vector<uint64_t> vector{};
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_file};
    for (auto & [seq] : fin)
    {
        int distance = 0;

        auto rep1 = seq | compare_view;
        auto representative = rep1 | input_view;
        auto rep_it = representative.begin();
        auto compare = seq | compare_view;
        auto comp_it = compare.begin();
        do
        {
            if (*rep_it == *comp_it)
            {
                if (comp_it != compare.begin())
                {
                    distances.push_back(distance);
                    distance = 0;
                }
                rep_it++;
            }
            else
            {
                distance++;
            }
            comp_it++;
        }
        while((rep_it != representative.end()) & (comp_it != compare.end()));
    }

    std::ofstream outfile;
    outfile.open(method_name + "_"+ std::string{sequence_file.stem()} + "_distances.out");
    for (int i = *std::min_element(distances.begin(), distances.end()); i <= *std::max_element(distances.begin(), distances.end()); i++)
    {
        auto count = std::count(distances.begin(), distances.end(), i);
        if (count > 0)
            outfile << i << "\t" << count << "\n";
    }
    outfile.close();
}

template <typename urng_t>
void distance_syncmer(std::filesystem::path sequence_file, urng_t input_view, range_arguments & args, std::string method_name)
{
    std::vector<uint64_t> distances{};
    int distance = 0;
    int it_1 = 0;
    int it_2 = 0;

    std::vector<uint64_t> vector{};
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_file};
    for (auto & [seq] : fin)
    {
        int distance = 0;
        auto representative = seq | input_view;
        auto rep_it = representative.begin();
        do
        {
            if (syncmer_filter(*rep_it, args.w_size.get(), (args.k_size *args.order),  args.positions, args.seed_se.get()))
            {
                if (rep_it != representative.begin())
                {
                    distances.push_back(distance);
                    distance = 0;
                }
            }
            else
            {
                distance++;
            }
            rep_it++;
        }
        while(rep_it != representative.end());
    }

    std::ofstream outfile;
    outfile.open(method_name + "_"+ std::string{sequence_file.stem()} + "_distances.out");
    for (int i = *std::min_element(distances.begin(), distances.end()); i <= *std::max_element(distances.begin(), distances.end()); i++)
    {
        auto count = std::count(distances.begin(), distances.end(), i);
        if (count > 0)
            outfile << i << "\t" << count << "\n";
    }
    outfile.close();
}

void fill_positions(std::vector<bool> & positions, int pos, int match_length)
{
    for(int j = pos; j < pos+match_length; j++)
        positions[j] = true;
}

/*! \brief Count the number of matches and the match coverage found in two sequence files.
 *  \param sequence_file1 The first sequence file.
 *  \param sequence_file2 The second sequence file.
 *  \param input_view View that should be tested.
 *  \param method_name Name of the tested method.
 *  \param args The arguments about the view to be used, needed for strobemers.
 */
template <typename urng_t, int strobemers = 0>
void match(std::filesystem::path sequence_file1, std::filesystem::path sequence_file2, urng_t input_view, std::string method_name, range_arguments & args)
{
    std::vector<uint64_t> seq1_vector{};
    std::vector<uint64_t> seq2_vector{};
    uint64_t matches{0};
    uint64_t missed{0};
    if constexpr (strobemers > 0)
    {
        seq1_vector = read_seq_file<strobemers>(sequence_file1, args);
        seq2_vector = read_seq_file<strobemers>(sequence_file2, args);
    }
    else
    {
        seq1_vector = read_seq_file(sequence_file1, input_view, args);
        seq2_vector = read_seq_file(sequence_file2, input_view, args);
    }

    std::size_t length{0};
    switch(args.name)
    {
        case kmer: length = seq1_vector.size()+args.shape.size()-1;
                        break;
        case strobemer: length = seq1_vector.size()+(args.k_size * args.order)-1;
    }
    std::vector<bool> positions(length, false);
    std::vector<uint64_t> islands{};
    uint64_t current_island{0};
    bool new_island{false};

    for(int i = 0; i < seq1_vector.size(); ++i)
    {
        if (seq1_vector[i] == seq2_vector[i])
        {
            matches++;
            switch(args.name)
            {
                case kmer: fill_positions(positions, i, args.shape.size());
                                break;
                case strobemer: fill_positions(positions, i, (args.k_size * args.order));
            }
        }
        else
        {
            missed++;
        }
    }

    for (int i = 0; i < positions.size(); ++i)
    {
        if (!positions[i])
        {
            current_island++;
        }
        else if (i > 0)
        {
            islands.push_back(current_island);
            current_island = 0;
        }
    }
    islands.push_back(current_island);


    double mean_island, stdev_island;
    get_mean_and_var(islands, mean_island, stdev_island);

    std::cout << "Matches: " << matches << "\t" << "Missed: " << missed << "\n";
    std::cout << "Match Coverage: " << std::count(positions.begin(), positions.end(), true)*100.0/positions.size() << "\n";
    std::cout << "Islands: " << *std::min_element(islands.begin(), islands.end()) << "\t" << mean_island << "\t" << stdev_island << "\t" << *std::max_element(islands.begin(), islands.end()) << "\n";
    std::cout << "Expected Island Size: " << get_expected(islands) << "\n";
}

/*! \brief Do the actual matching for repesentative methods.
 *  \param seq1_vector The first vector based on the first file.
 *  \param seq2_vector The second vector based on the second file.
 *  \param all1_vector The first vector containig all submers (kmers or strobmers) of the first file.
 *  \param all2_vector The second vector containig all submers (kmers or strobmers) of the second file.
 *  \param method_name Name of the tested method.
 *  \param args The arguments about the view to be used, needed for strobemers.
 */
void match_vectors(std::vector<uint64_t> & seq1_vector,
                   std::vector<uint64_t> & seq2_vector,
                   std::vector<uint64_t> & all1_vector,
                   std::vector<uint64_t> & all2_vector,
                   std::string method_name,
                   range_arguments & args)
{
    uint64_t matches{0};
    uint64_t missed{0};
    int it_1{0};
    int it_2{0};
    bool changed{true};
    int i{0};

    std::size_t length{0};
    length = all1_vector.size()+args.shape.size()-1;
    std::vector<bool> positions(length, false);
    std::vector<uint64_t> islands{};
    uint64_t current_island{0};
    bool new_island{false};

    while((it_1 < seq1_vector.size()) & (it_2 < seq2_vector.size()) & (i < all1_vector.size()))
    {
        if ((seq1_vector[it_1] == seq2_vector[it_2]) & changed)
        {
            matches++;
            fill_positions(positions, i, args.k_size);
        }

        if (seq1_vector[it_1] == all1_vector[i])
        {
            it_1++;
        }
        if (seq2_vector[it_2] == all2_vector[i])
        {
            it_2++;
        }
        i++;
    }

    for (int i = 0; i < positions.size(); ++i)
    {
        if (!positions[i])
        {
            current_island++;
        }
        else if (i > 0)
        {
            islands.push_back(current_island);
            current_island = 0;
        }
    }
    islands.push_back(current_island);

    double mean_island, stdev_island;
    get_mean_and_var(islands, mean_island, stdev_island);

    missed = std::min(seq1_vector.size(), seq2_vector.size()) - matches;
    std::cout << "Matches: " << matches << "\t" << "Missed: " << missed << "\n";
    std::cout << "Match Coverage: " << std::count(positions.begin(), positions.end(), true)*100.0/positions.size() << "\n";
    std::cout << "Islands: " << *std::min_element(islands.begin(), islands.end()) << "\t" << mean_island << "\t" << stdev_island << "\t" << *std::max_element(islands.begin(), islands.end()) << "\n";
    std::cout << "Expected Island Size: " << get_expected(islands) << "\n";
}

/*! \brief Count the number of matches and match coverage found in two sequence files.
 *  \param sequence_file1 The first sequence file.
 *  \param sequence_file2 The second sequence file.
 *  \param input_view View that should be tested.
 *  \param compare_view View for comparison, should be kmer_hash view.
 *  \param method_name Name of the tested method.
 *  \param args The arguments about the view to be used, needed for strobemers.
 */
template <typename urng_t, typename urng2_t>
void match(std::filesystem::path sequence_file1, std::filesystem::path sequence_file2, urng_t input_view, urng2_t compare_view, std::string method_name, range_arguments & args)
{
    std::vector<uint64_t> seq1_vector = read_seq_file(sequence_file1, input_view, args);
    std::vector<uint64_t> seq2_vector = read_seq_file(sequence_file2, input_view, args);
    std::vector<uint64_t> all1_vector = read_seq_file(sequence_file1, compare_view, args);
    std::vector<uint64_t> all2_vector = read_seq_file(sequence_file2, compare_view, args);

    match_vectors(seq1_vector, seq2_vector, all1_vector, all2_vector, method_name, args);
}

/*! \brief Count the number of matches and match coverage found in two sequence files.
 *  \param sequence_file1 The first sequence file.
 *  \param sequence_file2 The second sequence file.
 *  \param input_view View that should be tested.
 *  \param compare_view View for comparison, should be kmer_hash view.
 *  \param method_name Name of the tested method.
 *  \param args The arguments about the view to be used, needed for strobemers.
 */
template <typename urng_t, typename urng2_t>
void match_strobemer(std::filesystem::path sequence_file1, std::filesystem::path sequence_file2, urng_t input_view, urng2_t compare_view, std::string method_name, range_arguments & args)
{
    std::vector<uint64_t> seq1_vector = read_seq_file(sequence_file1, compare_view, input_view);
    std::vector<uint64_t> seq2_vector = read_seq_file(sequence_file2, compare_view, input_view);
    std::vector<uint64_t> all1_vector = read_seq_file(sequence_file1, compare_view, args);
    std::vector<uint64_t> all2_vector = read_seq_file(sequence_file2, compare_view, args);

    match_vectors(seq1_vector, seq2_vector, all1_vector, all2_vector, method_name, args);
}

/*! \brief Count the number of matches and match coverage found in two sequence files for syncmers bases on strobemers.
 *  \param sequence_file1 The first sequence file.
 *  \param sequence_file2 The second sequence file.
 *  \param input_view View that should be tested.
 *  \param method_name Name of the tested method.
 *  \param args The arguments about the view to be used, needed for strobemers.
 */
template <typename urng_t>
void match_syncmer(std::filesystem::path sequence_file1, std::filesystem::path sequence_file2, urng_t input_view, std::string method_name, range_arguments & args)
{
    std::vector<uint64_t> seq1_vector = read_seq_file<urng_t, true>(sequence_file1, input_view, args);
    std::vector<uint64_t> seq2_vector = read_seq_file<urng_t, true>(sequence_file2, input_view, args);
    std::vector<uint64_t> all1_vector = read_seq_file(sequence_file1, input_view, args);
    std::vector<uint64_t> all2_vector = read_seq_file(sequence_file2, input_view, args);

    match_vectors(seq1_vector, seq2_vector, all1_vector, all2_vector, method_name, args);
}

/*! \brief Function, that measures the speed of a method.
 *  \param sequence_files A vector of sequence files.
 *  \param input_view View that should be tested.
 *  \param method_name Name of the tested method.
 *  \param args The arguments about the view to be used, needed for strobemers.
 */
template <typename urng_t, int strobemers = 0>
void speed(std::vector<std::filesystem::path> sequence_files, urng_t input_view, std::string method_name, range_arguments & args)
{
   std::vector<int> speed_results{};
   std::ofstream outfile;
   int count{};
   for (int i = 0; i < sequence_files.size(); ++i)
   {
       robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
       auto start = std::chrono::high_resolution_clock::now();
       auto end = std::chrono::high_resolution_clock::now();
       auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
       if constexpr (strobemers > 0)
       {
           seqan3::sequence_file_input<my_traits2, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
           for (auto & [seq] : fin)
           {
               std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> strobes_vector;
               start = std::chrono::high_resolution_clock::now();
               get_strobemers<strobemers>(seq, args, strobes_vector);
               for (auto & t : strobes_vector) // iterate over the strobemer tuples
                    count += std::get<0>(t);
               end = std::chrono::high_resolution_clock::now();
               duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
               speed_results.push_back(duration.count());
           }
       }
       else
       {
           seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
           for (auto & [seq] : fin)
           {
               start = std::chrono::high_resolution_clock::now();
               for (auto && hash : seq | input_view)
                   count += hash; // Store hash value to enforce evaluation of it
               end = std::chrono::high_resolution_clock::now();
               duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
               speed_results.push_back(duration.count());

           }
       }
   }

   double mean_speed, stdev_speed;
   get_mean_and_var(speed_results, mean_speed, stdev_speed);

   // Store speed, the count value is stored so the compiler can not optimize the speed by not calculating the hash values
   outfile.open(std::string{args.path_out} + method_name + "_speed.out");
   outfile << method_name << "\t" << *std::min_element(speed_results.begin(), speed_results.end()) << "\t" << mean_speed << "\t" << stdev_speed << "\t" << *std::max_element(speed_results.begin(), speed_results.end()) << "\t" << count << "\n";
   outfile.close();
}

// Input files should be the output files from count
void unique(std::vector<std::filesystem::path> input_files, std::filesystem::path oname)
{
   std::ifstream infile;
   std::ofstream outfile;
   outfile.open(oname);

   for (int i = 0; i < input_files.size(); ++i)
   {
       //robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
       infile.open(input_files[i], std::ios::binary);
       uint64_t submer;
       uint16_t submer_count;
       uint64_t singletons{0};
       uint64_t all_counts{0};
       while(infile.read((char*)&submer, sizeof(submer)))
       {
            infile.read((char*)&submer_count, sizeof(submer_count));
	    if (submer_count == 1)
            	++singletons;
	    all_counts++;
       }
       infile.close();

       outfile << input_files[i].stem() << "\t" << (singletons * 100.0)/all_counts << "\n";

   }
   outfile.close();
}

std::string create_name(range_arguments & args, bool underlying_strobemer)
{
    std::string prefix{""};
    if (underlying_strobemer)
    {
        prefix = "Strobemer_";
    }

    switch(args.name)
    {
        case kmer: return "kmer_hash_"+std::to_string(args.k_size);
                   break;
        case minimiser: return prefix +"minimiser_hash_" + std::to_string(args.k_size) + "_" + std::to_string(args.w_size.get());
                        break;
        case modmers: return prefix +"modmer_hash_" + std::to_string(args.k_size) + "_" + std::to_string(args.w_size.get());
                        break;
        case strobemer: {
                            std::ranges::empty_view<seqan3::detail::empty_type> empty{};
                            if (args.lib_implementation)
                            {
                                if (args.rand & (args.order == 2))
                                    return "randstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                                else if (args.rand & (args.order == 3))
                                    return "randstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                                else if (args.hybrid)
                                    return "hybridstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                                else if (args.minstrobers & (args.order == 2))
                                    return "minstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                            }
                            else
                            {
                                if (args.hybrid & (args.order == 2))
                                    return "hybridstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                                else if (args.hybrid & (args.order == 3))
                                    return "hybridstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                                else if (args.minstrobers & (args.order == 2))
                                    return "minstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                                else if (args.minstrobers & (args.order == 3))
                                    return "minstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                                else if (args.rand & (args.order == 2))
                                    return "randstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                                else if (args.rand & (args.order == 3))
                                    return "randstrobemers_" + std::to_string(args.k_size) + "_" + std::to_string(args.order) + "_" +  std::to_string(args.w_min) + "_" +  std::to_string(args.w_max);
                                else
                                    return "";
                            }
                    }
        case syncmer: return prefix + "syncmer_hash_" + std::to_string(args.k_size) + "_" + std::to_string(args.w_size.get())+ "_" + std::to_string(args.positions[0]) + "_" + std::to_string(args.positions[args.positions.size()-1]);

        default: return "";
    }
}

void do_accuracy(accuracy_arguments & args)
{
    switch(args.name)
    {
        case kmer: accuracy(seqan3::views::kmer_hash(args.shape), create_name(args), args);
                        break;
        case minimiser: accuracy(seqan3::views::minimiser_hash(args.shape,
                                args.w_size, args.seed_se), create_name(args), args);
                        break;
        case modmers: accuracy(modmer_hash(args.shape,
                                args.w_size.get(), args.seed_se), create_name(args), args);
                        break;
        case syncmer: accuracy(syncmer_hash(args.w_size.get(), args.k_size, args.positions, args.seed_se),
                               create_name(args), args);
                        break;
        case strobemer: {
                            if (args.hybrid & (args.order == 2))
                                accuracy(hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                            else if (args.hybrid & (args.order == 3))
                                accuracy(hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                            else if (args.minstrobers & (args.order == 2))
                                accuracy(minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                            else if (args.minstrobers & (args.order == 3))
                                accuracy(minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                            else if (args.rand & (args.order == 2))
                                accuracy(randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                            else if (args.rand & (args.order == 3))
                                accuracy(randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                        }
    }
}

void do_counts(std::vector<std::filesystem::path> sequence_files, range_arguments & args, bool underlying_strobemer)
{
    if(underlying_strobemer)
    {
        switch(args.name)
        {
            case minimiser: {
                                if (args.hybrid & (args.order == 2))
                                    counts_strobemer(sequence_files, hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1), create_name(args, true), args);
                                if (args.hybrid & (args.order == 3))
                                    counts_strobemer(sequence_files, hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), minions::views::minimiser(args.w_size.get()-(args.shape.size()*3)+1), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 2))
                                    counts_strobemer(sequence_files, minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 3))
                                    counts_strobemer(sequence_files, minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), minions::views::minimiser(args.w_size.get()-(args.shape.size()*3)+1), create_name(args, true), args);
                                if (args.rand & (args.order == 2))
                                    counts_strobemer(sequence_files, randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1), create_name(args, true), args);
                                if (args.rand & (args.order == 3))
                                    counts_strobemer(sequence_files, randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), minions::views::minimiser(args.w_size.get()-(args.shape.size()*3)+1), create_name(args, true), args);
                            }
                            break;
            case modmers: {
                                if (args.hybrid & (args.order == 2))
                                    counts_strobemer(sequence_files, hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), modmer(args.w_size.get()), create_name(args, true), args);
                                if (args.hybrid & (args.order == 3))
                                    counts_strobemer(sequence_files, hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), modmer(args.w_size.get()), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 2))
                                    counts_strobemer(sequence_files, minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), modmer(args.w_size.get()), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 3))
                                    counts_strobemer(sequence_files, minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), modmer(args.w_size.get()), create_name(args, true), args);
                                if (args.rand & (args.order == 2))
                                    counts_strobemer(sequence_files, randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), modmer(args.w_size.get()), create_name(args, true), args);
                                if (args.rand & (args.order == 3))
                                    counts_strobemer(sequence_files, randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), modmer(args.w_size.get()), create_name(args, true), args);
                            }
                            break;
            case syncmer:  {
                                auto syncmer = std::views::filter([args] (uint64_t i)
                                                     {return syncmer_filter(i,args.w_size.get(), (args.k_size *args.order),  args.positions, args.seed_se.get());});
                                if (args.hybrid & (args.order == 2))
                                    counts_strobemer(sequence_files, hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), syncmer, create_name(args, true), args);
                                if (args.hybrid & (args.order == 3))
                                    counts_strobemer(sequence_files, hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), syncmer, create_name(args, true), args);
                                if (args.minstrobers & (args.order == 2))
                                    counts_strobemer(sequence_files, minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), syncmer, create_name(args, true), args);
                                if (args.minstrobers & (args.order == 3))
                                    counts_strobemer(sequence_files, minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), syncmer, create_name(args, true), args);
                                if (args.rand & (args.order == 2))
                                    counts_strobemer(sequence_files, randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), syncmer, create_name(args, true), args);
                                if (args.rand & (args.order == 3))
                                    counts_strobemer(sequence_files, randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), syncmer, create_name(args, true), args);
                           }
                           break;
        }
    }
    else
    {
        switch(args.name)
        {
            case kmer: counts(sequence_files, seqan3::views::kmer_hash(args.shape), create_name(args), args);
                       break;
            case minimiser: counts(sequence_files, seqan3::views::minimiser_hash(args.shape,
                                    args.w_size, args.seed_se), create_name(args), args);
                            break;
            case modmers: counts(sequence_files, modmer_hash(args.shape,
                                    args.w_size.get(), args.seed_se), create_name(args), args);
                            break;
            case syncmer:  counts(sequence_files, syncmer_hash(args.w_size.get(), args.k_size, args.positions, args.seed_se),
                                  create_name(args), args);
                            break;
            case strobemer: {
                                if (args.hybrid & (args.order == 2))
                                    counts(sequence_files, hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                                else if (args.hybrid & (args.order == 3))
                                    counts(sequence_files, hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                                else if (args.minstrobers & (args.order == 2))
                                    counts(sequence_files, minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                                else if (args.minstrobers & (args.order == 3))
                                    counts(sequence_files, minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                                else if (args.rand & (args.order == 2))
                                    counts(sequence_files, randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                                else if (args.rand & (args.order == 3))
                                    counts(sequence_files, randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                            }
        }
    }
}

void do_distance(std::filesystem::path sequence_file, range_arguments & args, bool underlying_strobemer)
{
    if (underlying_strobemer)
    {
        switch(args.name)
        {
            case minimiser: {
                                if (args.hybrid & (args.order == 2))
                                    distance_strobemer(sequence_file, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.hybrid & (args.order == 3))
                                    distance_strobemer(sequence_file, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.minstrobers & (args.order == 2))
                                    distance_strobemer(sequence_file, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.minstrobers & (args.order == 3))
                                    distance_strobemer(sequence_file, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.rand & (args.order == 2))
                                    distance_strobemer(sequence_file, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.rand & (args.order == 3))
                                    distance_strobemer(sequence_file, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                            }
                            break;
            case modmers: {
                                if (args.hybrid & (args.order == 2))
                                    distance_strobemer(sequence_file, modmer(args.w_size.get()), hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.hybrid & (args.order == 3))
                                    distance_strobemer(sequence_file, modmer(args.w_size.get()), hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.minstrobers & (args.order == 2))
                                    distance_strobemer(sequence_file, modmer(args.w_size.get()), minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.minstrobers & (args.order == 3))
                                    distance_strobemer(sequence_file, modmer(args.w_size.get()), minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.rand & (args.order == 2))
                                    distance_strobemer(sequence_file, modmer(args.w_size.get()),randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                                if (args.rand & (args.order == 3))
                                    distance_strobemer(sequence_file, modmer(args.w_size.get()), randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), std::string{args.path_out} + create_name(args, true));
                            }
                            break;
            case syncmer:  {
                                if (args.hybrid & (args.order == 2))
                                    distance_syncmer(sequence_file, hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), args, std::string{args.path_out} + create_name(args, true));
                                if (args.hybrid & (args.order == 3))
                                    distance_syncmer(sequence_file, hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), args, std::string{args.path_out} + create_name(args, true));
                                if (args.minstrobers & (args.order == 2))
                                    distance_syncmer(sequence_file, minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), args, std::string{args.path_out} + create_name(args, true));
                                if (args.minstrobers & (args.order == 3))
                                    distance_syncmer(sequence_file, minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), args, std::string{args.path_out} + create_name(args, true));
                                if (args.rand & (args.order == 2))
                                    distance_syncmer(sequence_file, randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), args, std::string{args.path_out} + create_name(args, true));
                                if (args.rand & (args.order == 3))
                                    distance_syncmer(sequence_file, randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), args, std::string{args.path_out} + create_name(args, true));
                           }
                           break;
        }
    }
    else
    {
        switch(args.name)
        {
            case minimiser: distance(sequence_file, seqan3::views::minimiser_hash(args.shape, args.w_size, args.seed_se),
                                                    seqan3::views::minimiser_hash(args.shape, seqan3::window_size{args.shape.size()}, args.seed_se), std::string{args.path_out} + create_name(args));
                            break;
            case modmers: distance(sequence_file, modmer_hash(args.shape, args.w_size.get(), args.seed_se),
                                                  modmer_hash(args.shape, 1, args.seed_se), std::string{args.path_out} + create_name(args));
                            break;
            case syncmer: distance(sequence_file, syncmer_hash(args.w_size.get(), args.k_size, args.positions, args.seed_se),
                                                  seqan3::views::minimiser_hash(args.shape, seqan3::window_size{args.shape.size()}, args.seed_se), std::string{args.path_out} + create_name(args));
                          break;
        }
    }
}

void do_match(std::filesystem::path sequence_file1, std::filesystem::path sequence_file2, range_arguments & args, bool underlying_strobemer)
{
    seqan3::seed seed{args.seed_se};

    if(underlying_strobemer)
    {
        switch(args.name)
        {
            case minimiser: {
                                if (args.hybrid & (args.order == 2))
                                    match_strobemer(sequence_file1, sequence_file2, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.hybrid & (args.order == 3))
                                    match_strobemer(sequence_file1, sequence_file2, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 2))
                                    match_strobemer(sequence_file1, sequence_file2, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 3))
                                    match_strobemer(sequence_file1, sequence_file2, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.rand & (args.order == 2))
                                    match_strobemer(sequence_file1, sequence_file2, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.rand & (args.order == 3))
                                    match_strobemer(sequence_file1, sequence_file2, minions::views::minimiser(args.w_size.get()-(args.shape.size()*2)+1),randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                            }
                            break;
            case modmers: {
                                if (args.hybrid & (args.order == 2))
                                    match_strobemer(sequence_file1, sequence_file2, modmer(args.w_size.get()), hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.hybrid & (args.order == 3))
                                    match_strobemer(sequence_file1, sequence_file2, modmer(args.w_size.get()), hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 2))
                                    match_strobemer(sequence_file1, sequence_file2, modmer(args.w_size.get()), minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 3))
                                    match_strobemer(sequence_file1, sequence_file2, modmer(args.w_size.get()), minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.rand & (args.order == 2))
                                    match_strobemer(sequence_file1, sequence_file2, modmer(args.w_size.get()),randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.rand & (args.order == 3))
                                    match_strobemer(sequence_file1, sequence_file2, modmer(args.w_size.get()), randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                            }
                            break;
            case syncmer:  {
                                if (args.hybrid & (args.order == 2))
                                    match_syncmer(sequence_file1, sequence_file2, hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.hybrid & (args.order == 3))
                                    match_syncmer(sequence_file1, sequence_file2, hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 2))
                                    match_syncmer(sequence_file1, sequence_file2, minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.minstrobers & (args.order == 3))
                                    match_syncmer(sequence_file1, sequence_file2, minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.rand & (args.order == 2))
                                    match_syncmer(sequence_file1, sequence_file2, randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                                if (args.rand & (args.order == 3))
                                    match_syncmer(sequence_file1, sequence_file2, randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args, true), args);
                           }
                           break;
        }
    }
    else
    {
        switch(args.name)
        {
            case kmer: match(sequence_file1, sequence_file2, seqan3::views::kmer_hash(args.shape), create_name(args), args);
                       break;
            case minimiser: match(sequence_file1, sequence_file2, seqan3::views::minimiser_hash(args.shape,
                                    args.w_size, args.seed_se),seqan3::views::minimiser_hash(args.shape,
                                                            seqan3::window_size{args.shape.size()}, args.seed_se), create_name(args), args);
                            break;
            case modmers: match(sequence_file1, sequence_file2, modmer_hash(args.shape,
                                    args.w_size.get(), args.seed_se), modmer_hash(args.shape, 1, args.seed_se), create_name(args), args);
                            break;
            case syncmer: match(sequence_file1, sequence_file2, syncmer_hash(args.w_size.get(), args.k_size, args.positions, args.seed_se), seqan3::views::minimiser_hash(args.shape,
                                                            seqan3::window_size{args.shape.size()}, args.seed_se), create_name(args), args);
                            break;
            case strobemer: std::ranges::empty_view<seqan3::detail::empty_type> empty{};
                            if (args.lib_implementation)
                            {
                                if (args.rand & (args.order == 2))
                                    match<std::ranges::empty_view<seqan3::detail::empty_type>, 1>(sequence_file1, sequence_file2, empty, create_name(args), args);
                                else if (args.rand & (args.order == 3))
                                    match<std::ranges::empty_view<seqan3::detail::empty_type>, 2>(sequence_file1, sequence_file2, empty, create_name(args), args);
                                else if (args.hybrid)
                                    match<std::ranges::empty_view<seqan3::detail::empty_type>, 3>(sequence_file1, sequence_file2, empty, create_name(args), args);
                                else if (args.minstrobers & (args.order == 2))
                                    match<std::ranges::empty_view<seqan3::detail::empty_type>, 4>(sequence_file1, sequence_file2, empty, create_name(args), args);
                            }
                            else
                            {
                                if (args.hybrid & (args.order == 2))
                                    match(sequence_file1, sequence_file2, hybridstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                                else if (args.hybrid & (args.order == 3))
                                    match(sequence_file1, sequence_file2, hybridstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se),create_name(args), args);
                                else if (args.minstrobers & (args.order == 2))
                                    match(sequence_file1, sequence_file2, minstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                                else if (args.minstrobers & (args.order == 3))
                                    match(sequence_file1, sequence_file2, minstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                                else if (args.rand & (args.order == 2))
                                    match(sequence_file1, sequence_file2, randstrobe2_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                                else if (args.rand & (args.order == 3))
                                    match(sequence_file1, sequence_file2, randstrobe3_hash(args.shape, args.w_min, args.w_max, args.seed_se), create_name(args), args);
                        }
                        break;
        }
    }
}

// Note: Speed is based on non-canonical version!
void do_speed(std::vector<std::filesystem::path> sequence_files, range_arguments & args)
{
    switch(args.name)
    {
        case kmer: speed(sequence_files, seqan3::views::kmer_hash(args.shape), create_name(args), args);
                   break;
        case minimiser: speed(sequence_files, seqan3::views::kmer_hash(args.shape) | seqan3::views::minimiser(args.w_size.get()-args.shape.size()+1), create_name(args), args);
                        break;
        case modmers: speed(sequence_files, seqan3::views::kmer_hash(args.shape) | modmer(args.w_size.get()), create_name(args), args);
                        break;
        case strobemer: {std::ranges::empty_view<seqan3::detail::empty_type> empty{};
                        if (args.lib_implementation)
                        {
                            if (args.rand & (args.order == 2))
                                speed<std::ranges::empty_view<seqan3::detail::empty_type>, 1>(sequence_files, empty, create_name(args), args);
                            else if (args.rand & (args.order == 3))
                                speed<std::ranges::empty_view<seqan3::detail::empty_type>, 2>(sequence_files, empty, create_name(args), args);
                            else if (args.hybrid)
                                speed<std::ranges::empty_view<seqan3::detail::empty_type>, 3>(sequence_files, empty, create_name(args), args);
                            else if (args.minstrobers & (args.order == 2))
                                speed<std::ranges::empty_view<seqan3::detail::empty_type>, 4>(sequence_files, empty, create_name(args), args);
                        }
                        else
                        {
                            if (args.hybrid & (args.order == 2))
                                speed(sequence_files, seqan3::views::kmer_hash(args.shape) | seqan3::views::hybridstrobe(args.w_min + args.shape.size() - 1, args.w_max - args.shape.size() + 1, args.shape.count()), create_name(args), args);
                            else if (args.hybrid & (args.order == 3))
                                speed(sequence_files, seqan3::views::kmer_hash(args.shape) | seqan3::views::hybridstrobe(true, args.w_min + args.shape.size() - 1, args.w_max - args.shape.size() + 1, args.shape.count()), create_name(args), args);
                            else if (args.minstrobers & (args.order == 2))
                                speed(sequence_files, seqan3::views::kmer_hash(args.shape) | seqan3::views::minstrobe(args.w_min + args.shape.size() - 1, args.w_max - args.shape.size() + 1, args.shape.count()), create_name(args), args);
                            else if (args.minstrobers & (args.order == 3))
                                speed(sequence_files, seqan3::views::kmer_hash(args.shape) | seqan3::views::minstrobe(true, args.w_min + args.shape.size() - 1, args.w_max - args.shape.size() + 1, args.shape.count()), create_name(args), args);
                            else if (args.rand & (args.order == 2))
                                speed(sequence_files, seqan3::views::kmer_hash(args.shape) | seqan3::views::randstrobe(args.w_min + args.shape.size() - 1, args.w_max - args.shape.size() + 1, args.shape.count()), create_name(args), args);
                            else if (args.rand & (args.order == 3))
                                speed(sequence_files, seqan3::views::kmer_hash(args.shape) | seqan3::views::randstrobe(true, args.w_min + args.shape.size() - 1, args.w_max - args.shape.size() + 1, args.shape.count()), create_name(args), args);
                    }
                    break;}
        case syncmer: speed(sequence_files, syncmer_hash_no_reverse(args.w_size.get(), args.k_size, args.positions), create_name(args), args);
    }
}
