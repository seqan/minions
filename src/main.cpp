#include <sstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "compare.h"


int speed(seqan3::argument_parser & parser)
{
    parser.add_positional_option(sequence_files, "Please provide at least one sequence file.");
    do_comparison(sequence_files);

    return 0;
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser top_level_parser{"minions", argc, argv, seqan3::update_notifications::on,
                                            {"speed"}};

    // Declarations for argument parser
    std::vector<std::filesystem::path> sequence_files{};

    // Parser
    parser.info.author = "Mitra Darvish"; // give parser some infos
    parser.info.version = "0.1.0";

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
