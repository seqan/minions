#include <sstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "compare.h"

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"comparison", argc, argv};

    // Declarations for argument parser
    std::vector<std::filesystem::path> sequence_files{};

    // Parser
    parser.info.author = "Mitra Darvish"; // give parser some infos
    parser.info.version = "0.1.0";
    parser.add_positional_option(sequence_files, "Please provide at least one sequence file.");

    try
    {
         parser.parse();                                                  // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n"; // give error message
        return -1;
    }

    compare(sequence_files, seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{19}}), "kmer_hash (19)         ");
    compare(sequence_files, seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{19}},
                                                          seqan3::window_size{19}), "minimiser_hash (19, 19)");
    compare(sequence_files, seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{19}},
                                                          seqan3::window_size{23}), "minimiser_hash (19, 23)");

    return 0;
}
