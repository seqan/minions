#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("minions");
    std::string expected
    {
        "minions\n"
        "=======\n"
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
        "Parsing error. You either forgot or misspelled the subcommand! Please specify which sub-program you want to "
        "use: one of [speed]. Use -h/--help for more information.\n"

    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, with_argument)
{
    cli_test_result result = execute_app("minions speed --method kmer -k 19", data("example1.fasta"));
    std::string expected
    {
        "kmer_hash_19         	159493	159493	0	159493\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    std::istringstream iss(result.out);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>());

    EXPECT_EQ(expected.substr(0,12), results[0]);
    EXPECT_EQ("0", results[3]);
}

TEST_F(cli_test, minimiser)
{
    cli_test_result result = execute_app("minions speed --method minimiser -k 19 -w 19 ", data("example1.fasta"));
    std::string expected
    {
        "minimiser_hash_19_19         	159493	159493	0	159493\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    std::istringstream iss(result.out);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>());

    EXPECT_EQ(expected.substr(0,20), results[0]);
    EXPECT_EQ("0", results[3]);
}

TEST_F(cli_test, strobemer)
{
    cli_test_result result = execute_app("minions speed --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --randstrobemers", data("example1.fasta"));
    std::string expected
    {
        "randstrobemers_19_2_16_30         	159493	159493	0	159493\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    std::istringstream iss(result.out);
    std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                     std::istream_iterator<std::string>());

    EXPECT_EQ(expected.substr(0,25), results[0]);
    EXPECT_EQ("0", results[3]);
}
