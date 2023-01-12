#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("minions unique");
    std::string expected
    {
        "minions-unique - Calculates the percentage of unique submers of a method for the given sequence files.\n"
        "======================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, with_argument)
{
    cli_test_result result = execute_app("minions unique --method kmer -k 19", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minimiser)
{
    cli_test_result result = execute_app("minions unique --method minimiser -k 19 -w 19 ", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, modmer)
{
    cli_test_result result = execute_app("minions unique --method modmer -k 19 -w 2 ", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, strobemer)
{
    cli_test_result result = execute_app("minions unique --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --rand", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, hybridstrobemer)
{
    cli_test_result result = execute_app("minions unique --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --hybrid", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minstrobemer)
{
    cli_test_result result = execute_app("minions unique --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --min", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, wrong_method)
{
    cli_test_result result = execute_app("minions unique --method submer -k 19", data("example1.fasta"));
    std::string expected
    {
        "Error. Incorrect command line input for unique. Validation failed "
        "for option --method: Value submer is not one of [kmer,minimiser,modmer,strobemer].\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
