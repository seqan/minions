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

    do_comparison(sequence_files);

    return 0;
}
