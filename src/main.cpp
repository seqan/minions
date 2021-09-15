#include <sstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "compare.h"

uint32_t w_size;
uint64_t shape{};
uint64_t se;

void read_range_arguments(seqan3::argument_parser & parser, range_arguments & args)
{
    parser.add_option(args.k_size, 'k', "kmer-size", "Define kmer size.");
    parser.add_option(w_size, 'w', "window", "Define window size. Default: 60.");
    parser.add_option(shape, '\0', "shape", "Define a shape by the decimal of a bitvector, where 0 symbolizes a "
                                           "position to be ignored, 1 a position considered. Default: ungapped.");
    parser.add_option(se, '\0', "seed", "Define seed.");
}

void parsing(range_arguments & args)
{
    args.w_size = seqan3::window_size{w_size};
    if (shape == 0)
        args.shape = seqan3::ungapped{args.k_size};
    else
        args.shape = seqan3::bin_literal{shape};
    args.seed_se = seqan3::seed{adjust_seed(args.k_size, se)};
}

int speed(seqan3::argument_parser & parser)
{
    range_arguments args{};
    std::vector<std::filesystem::path> sequence_files{};
    std::filesystem::path path_out{"./"};
    parser.add_positional_option(sequence_files, "Please provide at least one sequence file.");
    parser.add_option(path_out, 'o', "out", "Directory, where output files should be saved.");
    parser.add_flag(args.kmers, '\0', "kmers", "If k-mers should be calculated.");
    parser.add_flag(args.minimiser, '\0', "minimiser", "If minimisers should be calculated.");
    read_range_arguments(parser, args);

    try
    {
        parser.parse();
        parsing(args);
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Error. Incorrect command line input for speed. " << ext.what() << "\n";
        return -1;
    }

    try
    {
        do_comparison(sequence_files, args, path_out);
    }
    catch (const std::invalid_argument & e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser top_level_parser{"minions", argc, argv, seqan3::update_notifications::on,
                                            {"speed"}};

    // Parser
    top_level_parser.info.author = "Mitra Darvish"; // give parser some infos
    top_level_parser.info.version = "0.1.0";

    try
    {
         top_level_parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n"; // give error message
        return -1;
    }

    seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser

    if (sub_parser.info.app_name == std::string_view{"minions-speed"})
        speed(sub_parser);

    return 0;
}
