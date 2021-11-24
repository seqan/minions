#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("minions coverage");
    std::string expected
    {
        "minions-coverage\n"
        "================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, kmer)
{
    cli_test_result result = execute_app("minions coverage --method kmer -k 19", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, gapped_kmer)
{
    cli_test_result result = execute_app("minions coverage --method kmer --shape 524223", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minimiser)
{
    cli_test_result result = execute_app("minions coverage --method minimiser -k 19 -w 19 ", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, modmer)
{
    cli_test_result result = execute_app("minions coverage --method modmer -k 19 -w 19 ", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, wrong_method)
{
    cli_test_result result = execute_app("minions coverage --method submer -k 19", data("example1.fasta"));
    std::string expected
    {
        "Error. Incorrect command line input for coverage. Validation failed "
        "for option --method: Value submer is not one of [kmer,minimiser,modmer].\n"
    };

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.err, expected);
    EXPECT_EQ(result.out, std::string{});
}
