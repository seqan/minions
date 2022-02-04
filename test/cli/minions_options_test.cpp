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
        "Parsing error. You either forgot or misspelled the subcommand! Please "
        "specify which sub-program you want to use: one of [accuracy,counts,"
        "coverage,speed]. Use -h/--help for more information.\n"

    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}
