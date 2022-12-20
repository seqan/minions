#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("minions match");
    std::string expected
    {
        "minions-match - Counts the number of matches for a given method between the two given files.\n"
        "============================================================================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, with_argument)
{
    cli_test_result result = execute_app("minions match --method kmer -k 19", data("example1.fasta"), data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Matches: 6425136\tMissed: 0\nMatch Coverage: 100\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minimiser)
{
    cli_test_result result = execute_app("minions match --method minimiser -k 19 -w 19 ", data("example1.fasta"), data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Matches: 6425136\tMissed: 0\nMatch Coverage: 100\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, modmer)
{
    cli_test_result result = execute_app("minions match --method modmer -k 19 -w 2 ", data("example1.fasta"), data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Matches: 3139105\tMissed: 0\nMatch Coverage: 99.9723\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, strobemer)
{
    cli_test_result result = execute_app("minions match --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --randstrobemers", data("example1.fasta"), data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Matches: 3722976\tMissed: 0\nMatch Coverage: 100\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, hybridstrobemer)
{
    cli_test_result result = execute_app("minions match --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --hybrid", data("example1.fasta"), data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Matches: 3722976\tMissed: 0\nMatch Coverage: 100\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, minstrobers)
{
    cli_test_result result = execute_app("minions match --method strobemer -k 19 --w-min 16 --w-max 30 --order 2 --minstrobers", data("example1.fasta"), data("example1.fasta"));
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{"Matches: 3722976\tMissed: 0\nMatch Coverage: 100\n"});
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, wrong_method)
{
    cli_test_result result = execute_app("minions match --method submer -k 19", data("example1.fasta"), data("example1.fasta"));
    std::string expected
    {
        "Error. Incorrect command line input for match. Validation failed "
        "for option --method: Value submer is not one of [kmer,minimiser,modmer,strobemer].\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
