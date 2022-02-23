#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("minions accuracy");
    std::string expected
    {
        "minions-accuracy - Counts the true positive, false positive, true negatives and false negatives of a sequence file given the ground truth by a solution file.\n"
        "=============================================================================================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, kmer)
{
    cli_test_result result = execute_app("minions accuracy --method kmer -k 19 ", data("example.ibf"), "--search-file", data("search.fasta"), "--solution-file", data("expected_search_result.out"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minimiser)
{
    cli_test_result result = execute_app("minions accuracy --method minimiser -k 19 -w 19 ", data("example.ibf"), "--search-file", data("search.fasta"), "--solution-file", data("expected_search_result.out"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, gapped_minimiser)
{
    cli_test_result result = execute_app("minions accuracy --method minimiser -k 19 -w 19 --shape 524223 ", data("example.ibf"), "--search-file", data("search.fasta"), "--solution-file", data("expected_search_result.out"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, modmer)
{
    cli_test_result result = execute_app("minions accuracy --method modmer -k 19 -w 2 ", data("example.ibf"), "--search-file", data("search.fasta"), "--solution-file", data("expected_search_result.out"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, opensyncmer)
{
    cli_test_result result = execute_app("minions accuracy --method opensyncmer -k 19 -w 2 ", data("example.ibf"), "--search-file", data("search.fasta"), "--solution-file", data("expected_search_result.out"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, no_ibf_yet)
{
    cli_test_result result = execute_app("minions accuracy --method minimiser -k 19 -w 19 --ibfsize 10000 ", data("minimiser_hash_19_19_example1.out"), "--search-file", data("search.fasta"), "--solution-file", data("expected_search_result.out"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, wrong_method)
{
    cli_test_result result = execute_app("minions accuracy --method submer -k 19 --search-file ", data("search.fasta"), data("example.ibf"));
    std::string expected
    {
        "Error. Incorrect command line input for accuracy. Validation failed "
        "for option --method: Value submer is not one of [kmer,minimiser,modmer,opensyncmer].\n"
    };

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.err, expected);
    EXPECT_EQ(result.out, std::string{});
}
