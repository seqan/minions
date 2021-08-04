#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("Comparison");
    std::string expected
    {
        "comparison\n"
        "==========\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("Comparison", "-v");
    std::string expected
    {
        "Parsing error. Unknown option -v. In case this is meant to be a non-option/argument/parameter, please "
        "specify the start of non-options with '--'. See -h/--help for program information.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, with_argument)
{
    cli_test_result result = execute_app("Comparison", data("example1.fasta"));
    std::string expected
    {
        "kmer_hash (19)         	159493	159493	0	159493\n"
        "minimiser_hash (19, 19)	241282	241282	0	241282\n"
        "minimiser_hash (19, 23)	157535	157535	0	157535\n"

    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(expected.substr(0,14), result.out.substr(0,14));
    EXPECT_EQ(expected[38], result.out[38]);
    EXPECT_EQ(expected.substr(46,24), result.out.substr(46,24));
    EXPECT_EQ(expected[84], result.out[84]);
    EXPECT_EQ(expected.substr(93,24), result.out.substr(93,24));
    EXPECT_EQ(expected[132], result.out[132]);
    EXPECT_EQ(result.err, std::string{});
}
