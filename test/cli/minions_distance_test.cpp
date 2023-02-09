#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("minions distance");
    std::string expected
    {
        "minions-distance - Estimates the distance of the singular submers to each other for different methods.\n"
        "======================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minimiser)
{
    cli_test_result result = execute_app("minions distance --method minimiser -k 19 -w 19 ", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Distances: 0\t0\t0\t0\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, gapped_minimiser)
{
    cli_test_result result = execute_app("minions distance --method minimiser -k 19 -w 19 --shape 524223", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Distances: 0\t0\t0\t0\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, modmer)
{
    cli_test_result result = execute_app("minions distance --method modmer -k 19 -w 2", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Distances: 0\t1.03566\t1.51183\t49\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, syncmer)
{
    cli_test_result result = execute_app("minions distance --method syncmer -k 19 -w 2 -p 0", data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Distances: 0\t3.74393\t4.88286\t40\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, wrong_method)
{
    cli_test_result result = execute_app("minions distance --method submer -k 19", data("example1.fasta"));
    std::string expected
    {
        "Error. Incorrect command line input for distance. Validation failed "
        "for option --method: Value submer is not one of [minimiser,modmer,syncmer].\n"
    };

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.err, expected);
    EXPECT_EQ(result.out, std::string{});
}
