#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("minions counts");
    std::string expected
    {
        "minions-counts - Counts the number of submers in the given sequence files.\n"
        "==========================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, with_argument)
{
    cli_test_result result = execute_app("minions counts --method kmer -k 19", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minimiser)
{
    cli_test_result result = execute_app("minions counts --method minimiser -k 19 -w 19 ", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, modmer)
{
    cli_test_result result = execute_app("minions counts --method modmer -k 19 -w 2 ", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, strobemer)
{
    cli_test_result result = execute_app("minions counts --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --randstrobemers", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, hybridstrobemer)
{
    cli_test_result result = execute_app("minions counts --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --hybrid", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minstrobers)
{
    cli_test_result result = execute_app("minions counts --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --minstrobers", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, opensyncmer)
{
    cli_test_result result = execute_app("minions counts --method opensyncmer -k 19 -w 3", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, wrong_method)
{
    cli_test_result result = execute_app("minions counts --method submer -k 19", data("example1.fasta"));
    std::string expected
    {
        "Error. Incorrect command line input for counts. Validation failed "
        "for option --method: Value submer is not one of [kmer,minimiser,modmer,strobemer,opensyncmer].\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
