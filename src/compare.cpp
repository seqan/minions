#include <chrono>
#include <ranges>

#include <index.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>

#include "compare.h"
#include "hybridstrobe_hash.hpp"
#include "minimiser_hash_distance.hpp"
#include "minstrobe_hash.hpp"
#include "modmer_hash.hpp"
#include "modmer_hash_distance.hpp"
#include "randstrobe_hash.hpp"
#include "syncmer_hash.hpp"

/*! \brief Calculate mean and variance of given list.
 *  \param results The vector from which mean and varaince should be calculated of.
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
template <typename urng_t, int strobemers = 0>
void counts(std::vector<std::filesystem::path> sequence_files, urng_t input_view, std::string method_name, range_arguments & args)
{
    std::vector<uint64_t> counts_results{};
    std::ofstream outfile;
    for (int i = 0; i < sequence_files.size(); ++i)
    {
        robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
        if constexpr (strobemers > 0)
        {
            seqan3::sequence_file_input<my_traits2, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            for (auto & [seq] : fin)
            {
                std::vector<std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int>> strobes_vector;
                get_strobemers<strobemers>(seq, args, strobes_vector);
                for (auto & t : strobes_vector) // iterate over the strobemer tuples
                    hash_table[std::get<0>(t)] = std::min<uint16_t>(65534u, hash_table[std::get<0>(t)] + 1);
            }
        }
        else
        {
            seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_files[i]};
            for (auto & [seq] : fin)
            {
                for (auto && hash : seq | input_view)
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

/*! \brief Function, that increases how often a position is covered by one.
 *  \param covs A vector, where every entry represents a position in the sequence and how often it is covered.
 *  \param position Starting position.
 *  \param shape The positions of covered positions.
 */
void positions_covered(std::vector<uint32_t> & covs, int position, seqan3::shape shape)
{
    for(int pos = 0; pos < shape.size(); pos++)
    {
        if (shape[pos] == 1)
            covs[position+pos]++;
    }
}

/*! \brief Function, get the actual coverage of some submers.
 *  \param kmers A view that contains all possible k-mers (minimisers with window length = kmer length, if reverse strand should be considered).
 *  \param submers The submers, which coverage should be obtained.
 *  \param covs A vector, where every entry represents a position in the sequence and how often it is covered.
 *  \param shape Shape of a submer.
 */
template <typename urng_t>
void coverage(urng_t kmers, std::deque<uint64_t> & submers, std::vector<uint32_t> & covs, seqan3::shape shape)
{
    int i{};
    for (auto && elem : kmers)
    {
        if (elem == submers[0])
        {
            positions_covered(covs, i, shape);
            submers.pop_front();
        }

        if (submers.size() == 0)
            break;
        i++;
    }
}

/*! \brief Function, get the coverage of one sequence file for a method.
 *  Strobemers not supported, because strobemers should cover every position at least once because they are not
 *  a sampling method.
 *  \param sequence_file A sequence file.
 *  \param kmer_view View to compare to (Should be k-mers in most cases or minimisers with window length = kmer length).
 *  \param input_view View that should be evaluated.
 *  \param method_name Name of the tested method.
 *  \param args The arguments about the view to be used, needed for strobemers.
 */
template <typename urng_t, typename urng_t2>
void compare_cov(std::filesystem::path sequence_file, urng_t kmer_view, urng_t2 input_view, std::string method_name, range_arguments & args)
{
    std::vector<int> islands{};
    int island{};
    int covered{};
    std::vector<int> avg_islands{};
    std::vector<int> largest_islands{};
    std::vector<int> covered_percentage{};
    std::vector<int> covereage_avg{};
    std::ofstream outfile;

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::id,seqan3::field::seq>> fin{sequence_file};
    for (auto & [id, seq] : fin)
    {
        auto kmers = seq | kmer_view;
        auto submers = seq | input_view;

        std::vector<uint32_t> covs{};
        covs.assign(seq.size(), 0);
        std::deque<uint64_t> submers2{};

        for(auto && sub: submers)
            submers2.push_back(sub);
        coverage(kmers, submers2, covs, args.shape);

        int i{0};
        for(auto & elem : covs)
        {
            if (elem == 0)
            {
                island++;
            }
            else
            {
                if (island > 0)
                {
                    if (island >23)
                        std::cout << i << " " << island << ", " << id<< "\n";
                    islands.push_back(island);
                    island = 0;
                }
                covered++;
            }
            i++;
        }

        if (island > 0)
        {
            if (island >23)
                std::cout << i << " " << island << ", " << id<< "\n";
            islands.push_back(island);
            island = 0;
        }

        covered_percentage.push_back(covered);

        double mean_cov{};
        double var_cov{};
        get_mean_and_var(covs, mean_cov, var_cov);
        covereage_avg.push_back(mean_cov);

        if (islands.size() == 0)
        {
            avg_islands.push_back(0);
            largest_islands.push_back(0);
        }
        else
        {
            double mean{};
            double var{};
            get_mean_and_var(islands, mean, var);
            avg_islands.push_back(mean);
            largest_islands.push_back(*std::max_element(islands.begin(), islands.end()));
        }

        islands.clear();
        island = 0;
        covered = 0;
    }

    double mean_covered, mean_largest_island, mean_avg_island, stdev_covered, stdev_largest_island, stdev_avg_island;
    get_mean_and_var(covered_percentage, mean_covered, stdev_covered);
    get_mean_and_var(largest_islands, mean_largest_island, stdev_largest_island);
    get_mean_and_var(avg_islands, mean_avg_island, stdev_avg_island);

    std::nth_element(covereage_avg.begin(), covereage_avg.begin() + covereage_avg.size()/2, covereage_avg.end());
    int median = covereage_avg[covereage_avg.size()/2];

    // Store speed and counts
    outfile.open(std::string{args.path_out} + method_name + "_coverage.out");
    outfile << "Covered\t"<< method_name << "\t" << *std::min_element(covered_percentage.begin(), covered_percentage.end()) << "\t" << mean_covered << "\t" << stdev_covered << "\t" << *std::max_element(covered_percentage.begin(), covered_percentage.end()) << "\n";
    outfile << "Covered Median\t"<< method_name << "\t" << *std::min_element(covereage_avg.begin(), covereage_avg.end()) << "\t" << median << "\t" << *std::max_element(covereage_avg.begin(), covereage_avg.end()) << "\n";
    outfile << "Largest Island\t"<< method_name << "\t" << *std::min_element(largest_islands.begin(), largest_islands.end()) << "\t" << mean_largest_island << "\t" << stdev_largest_island << "\t" << *std::max_element(largest_islands.begin(), largest_islands.end()) << "\n";
    outfile << "Avg Island\t"<< method_name << "\t" << *std::min_element(avg_islands.begin(), avg_islands.end()) << "\t" << mean_avg_island << "\t" << stdev_avg_island << "\t" << *std::max_element(avg_islands.begin(), avg_islands.end()) << "\n";
    outfile.close();
}

template <typename urng_t>
void compare_cov2(std::filesystem::path sequence_file, urng_t distance_view, std::string method_name, range_arguments & args)
{
    std::vector<double> coverage{};
    std::vector<double> stdev{};
    std::ofstream outfile;

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_file};
    for (auto & [seq] : fin)
    {
        for (auto && hash : seq | distance_view)
            coverage.push_back(hash);
    }
    double mean_coverage, stdev_coverage;
    get_mean_and_var(coverage, mean_coverage, stdev_coverage);

    // Store coverage
    outfile.open(std::string{args.path_out} + method_name + "_coverage.out");
    outfile << method_name << "\t" << *std::min_element(coverage.begin(), coverage.end()) << "\t" << mean_coverage << "\t" << stdev_coverage << "\t" << *std::max_element(coverage.begin(), coverage.end()) << "\n";
    outfile.close();
}

template <typename urng_t>
std::vector<uint64_t> read_seq_file(std::filesystem::path sequence_file, urng_t input_view)
{
    std::vector<uint64_t> vector{};
    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{sequence_file};
    for (auto & [seq] : fin)
    {
        for (auto && hash : seq | input_view)
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
        seq1_vector = read_seq_file(sequence_file1, input_view);
        seq2_vector = read_seq_file(sequence_file2, input_view);
    }

    std::size_t length{0};
    switch(args.name)
    {
        case kmer: length = seq1_vector.size()+args.shape.size()-1;
                        break;
        case strobemer: length = seq1_vector.size()+(args.k_size * args.order)-1;
    }
    std::vector<bool> positions(length, false);

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
    std::cout << "Matches: " << matches << "\t" << "Missed: " << missed << "\n";
    std::cout << "Match Coverage: " << std::count(positions.begin(), positions.end(), true)*100.0/positions.size() << "\n";
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
    uint64_t matches{0};
    uint64_t missed{0};
    std::vector<uint64_t> seq1_vector = read_seq_file(sequence_file1, input_view);
    std::vector<uint64_t> seq2_vector = read_seq_file(sequence_file2, input_view);
    std::vector<uint64_t> all1_vector = read_seq_file(sequence_file1, compare_view);
    std::vector<uint64_t> all2_vector = read_seq_file(sequence_file2, compare_view);

    int it_1{0};
    int it_2{0};
    bool changed{true};
    int i{0};

    std::size_t length{0};
    length = all1_vector.size()+args.shape.size()-1;
    std::vector<bool> positions(length, false);

    while((it_1 < seq1_vector.size()) & (it_2 < seq2_vector.size()) & (i < all1_vector.size()))
    {
        if ((seq1_vector[it_1] == seq2_vector[it_2]) & changed)
        {
            matches++;
            changed = false;
            switch(args.name)
            {
                case minimiser: fill_positions(positions, i, args.w_size.get());
                                break;
                case modmers: fill_positions(positions, i, args.k_size);
            }
        }

        if (seq1_vector[it_1] == all1_vector[i])
        {
            it_1++;
            changed = true;
        }
        if (seq2_vector[it_2] == all2_vector[i])
        {
            it_2++;
            changed = true;
        }
        i++;
    }

    missed = std::min(seq1_vector.size(), seq2_vector.size()) - matches;
    std::cout << "Matches: " << matches << "\t" << "Missed: " << missed << "\n";
    std::cout << "Match Coverage: " << std::count(positions.begin(), positions.end(), true)*100.0/positions.size() << "\n";
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
                   count++;
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
                   count++;
               end = std::chrono::high_resolution_clock::now();
               duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
               speed_results.push_back(duration.count());
           }
       }
   }

   double mean_speed, stdev_speed;
   get_mean_and_var(speed_results, mean_speed, stdev_speed);

   // Store speed
   outfile.open(std::string{args.path_out} + method_name + "_speed.out");
   outfile << method_name << "\t" << *std::min_element(speed_results.begin(), speed_results.end()) << "\t" << mean_speed << "\t" << stdev_speed << "\t" << *std::max_element(speed_results.begin(), speed_results.end()) << "\n";
   outfile.close();
}

uint64_t count_singletons(robin_hood::unordered_node_map<uint64_t, uint16_t> & hash_table)
{
    uint64_t singletons{0};
    for (auto && hash: hash_table)
    {
        if (hash_table[hash.first] == 1)
            ++singletons;
    }
    return singletons;
}

void unique(std::vector<std::filesystem::path> sequence_files, std::string method_name, range_arguments & args)
{
   std::ifstream infile;
   std::ofstream outfile;
   outfile.open(std::string{args.path_out} + method_name + "_unique.out");

   for (int i = 0; i < sequence_files.size(); ++i)
   {
       robin_hood::unordered_node_map<uint64_t, uint16_t> hash_table{};
       infile.open(std::string{args.path_out} + method_name + "_"+ std::string{sequence_files[i].stem()} + "_counts.out", std::ios::binary);
       uint64_t submer;
       uint16_t submer_count;
       while(infile.read((char*)&submer, sizeof(submer)))
       {
            infile.read((char*)&submer_count, sizeof(submer_count));
            hash_table[submer] = submer_count;
       }
       infile.close();

       outfile << sequence_files[i].stem() << "\t" << (count_singletons(hash_table) * 100.0)/hash_table.size() << "\n";

   }
   outfile.close();
}

std::string create_name(range_arguments & args)
{
    switch(args.name)
    {
        case kmer: return "kmer_hash_"+std::to_string(args.k_size);
                   break;
        case minimiser: return "minimiser_hash_" + std::to_string(args.k_size) + "_" + std::to_string(args.w_size.get());
                        break;
        case modmers: return "modmer_hash_" + std::to_string(args.k_size) + "_" + std::to_string(args.w_size.get());
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
        case syncmer:  {
                            if (args.closed)
                                return "closedsyncmer_hash_" + std::to_string(args.k_size) + "_" + std::to_string(args.w_size.get())+ "_" + std::to_string(args.t);
                            else
                                return "opensyncmer_hash_" + std::to_string(args.k_size) + "_" + std::to_string(args.w_size.get())+ "_" + std::to_string(args.t);
                            break;
                        }
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
        case syncmer: {
                        if (args.closed)
                            accuracy(syncmer_hash<false>(args.w_size.get(), args.k_size, args.t, args.seed_se),
                                     create_name(args), args);
                        else
                            accuracy(syncmer_hash<true>(args.w_size.get(), args.k_size, args.t, args.seed_se),
                                     create_name(args), args);
                        break;
                    }
    }
}

void do_counts(std::vector<std::filesystem::path> sequence_files, range_arguments & args)
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
        case syncmer:  {
                            if (args.closed)
                                counts(sequence_files, syncmer_hash<false>(args.w_size.get(), args.k_size, args.t, args.seed_se),
                                create_name(args), args);
                            else
                                counts(sequence_files, syncmer_hash<true>(args.w_size.get(), args.k_size, args.t, args.seed_se),
                                create_name(args), args);
                            break;
                        }
        case strobemer: {
                            std::ranges::empty_view<seqan3::detail::empty_type> empty{};
                            if (args.rand & (args.order == 2))
                                counts<std::ranges::empty_view<seqan3::detail::empty_type>, 1>(sequence_files, empty, create_name(args), args);
                            else if (args.rand & (args.order == 3))
                                counts<std::ranges::empty_view<seqan3::detail::empty_type>, 2>(sequence_files, empty, create_name(args), args);
                            else if (args.hybrid)
                                counts<std::ranges::empty_view<seqan3::detail::empty_type>, 3>(sequence_files, empty, create_name(args), args);
                            else if (args.minstrobers)
                                counts<std::ranges::empty_view<seqan3::detail::empty_type>, 4>(sequence_files, empty, create_name(args), args);
                            break;
                        }
    }
}

void do_coverage(std::filesystem::path sequence_file, range_arguments & args)
{
    switch(args.name)
    {
        case minimiser: compare_cov2(sequence_file, minimiser_hash_distance(args.shape,
                                args.w_size, args.seed_se), create_name(args), args);
                        break;
        case modmers: compare_cov2(sequence_file, modmer_hash_distance(args.shape,
                                args.w_size.get(), args.seed_se), create_name(args), args);
                        break;
    }
}

void do_match(std::filesystem::path sequence_file1, std::filesystem::path sequence_file2, range_arguments & args)
{
    seqan3::seed seed{args.seed_se};
    switch(args.name)
    {
        case kmer: match(sequence_file1, sequence_file2, seqan3::views::kmer_hash(args.shape), create_name(args), args);
                   break;
        case minimiser: match(sequence_file1, sequence_file2, seqan3::views::minimiser_hash(args.shape,
                                args.w_size, args.seed_se),seqan3::views::minimiser_hash(args.shape,
                                                        seqan3::window_size{args.shape.size()}, args.seed_se), create_name(args), args);
                        break;
        case modmers: match(sequence_file1, sequence_file2, modmer_hash(args.shape,
                                args.w_size.get(), args.seed_se), modmer_hash(args.shape,
                                                        1, args.seed_se), create_name(args), args);
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
                                match(sequence_file1, sequence_file2, hybridstrobe2_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                            else if (args.hybrid & (args.order == 3))
                                match(sequence_file1, sequence_file2, hybridstrobe3_hash(args.shape, args.w_min, args.w_max),create_name(args), args);
                            else if (args.minstrobers & (args.order == 2))
                                match(sequence_file1, sequence_file2, minstrobe2_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                            else if (args.minstrobers & (args.order == 3))
                                match(sequence_file1, sequence_file2, minstrobe3_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                            else if (args.rand & (args.order == 2))
                                match(sequence_file1, sequence_file2, randstrobe2_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                            else if (args.rand & (args.order == 3))
                                match(sequence_file1, sequence_file2, randstrobe3_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                    }
    }
}


void do_speed(std::vector<std::filesystem::path> sequence_files, range_arguments & args)
{
    switch(args.name)
    {
        case kmer: speed(sequence_files, seqan3::views::kmer_hash(args.shape), create_name(args), args);
                   break;
        case minimiser: speed(sequence_files, seqan3::views::minimiser_hash(args.shape,
                                args.w_size, args.seed_se), create_name(args), args);
                        break;
        case modmers: speed(sequence_files, modmer_hash(args.shape,
                                args.w_size.get(), args.seed_se), create_name(args), args);
                        break;
        case strobemer: std::ranges::empty_view<seqan3::detail::empty_type> empty{};
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
                                speed(sequence_files, hybridstrobe2_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                            else if (args.hybrid & (args.order == 3))
                                speed(sequence_files, hybridstrobe3_hash(args.shape, args.w_min, args.w_max),create_name(args), args);
                            else if (args.minstrobers & (args.order == 2))
                                speed(sequence_files, minstrobe2_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                            else if (args.minstrobers & (args.order == 3))
                                speed(sequence_files, minstrobe3_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                            else if (args.rand & (args.order == 2))
                                speed(sequence_files, randstrobe2_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                            else if (args.rand & (args.order == 3))
                                speed(sequence_files, randstrobe3_hash(args.shape, args.w_min, args.w_max), create_name(args), args);
                    }
    }
}
