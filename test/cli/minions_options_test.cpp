#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("minions");
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
    cli_test_result result = execute_app("minions", "-v");
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
    cli_test_result result = execute_app("minions", data("example1.fasta"));
    std::string expected
    {
        "SPEED                  	Minimum	Mean	StdDev	Maximum\n"
        "kmer_hash (19)         	159493	159493	0	159493\n"
        "minimiser_hash (19, 19)	241282	241282	0	241282\n"
        "minimiser_hash (19, 23)	157535	157535	0	157535\n"

    };
    EXPECT_EQ(result.exit_code, 0);
    std::istringstream iss(result.out);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>());

    EXPECT_EQ(expected.substr(0,53), result.out.substr(0,53));
    EXPECT_EQ(expected.substr(52,9), results[5]);
    EXPECT_EQ("0", results[9]);
    EXPECT_EQ(expected.substr(99,14), results[11]);
    EXPECT_EQ("0", results[16]);
    EXPECT_EQ(expected.substr(146,14), results[18]);
    EXPECT_EQ("0", results[23]);
    EXPECT_EQ(result.err, std::string{});
}
