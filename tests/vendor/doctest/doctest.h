
// =================================================================================================
//
// doctest.h - the lightest feature-rich C++ single-header testing framework for unit tests and TDD
//
// Copyright (c) 2016-2017 Viktor Kirilov
//
// Distributed under the MIT Software License
// See accompanying LICENSE.txt file or copy at
// https://opensource.org/licenses/MIT
//
// The documentation can be found at the library's page:
// https://github.com/doctest/doctest
//
// =================================================================================================

#ifndef TESTS_VENDOR_DOCTEST_DOCTEST_H
#define TESTS_VENDOR_DOCTEST_DOCTEST_H

#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace doctest
{

struct TestCase
{
    const char* name;
    void (*func)();
};

inline std::vector<TestCase>& registry()
{
    static std::vector<TestCase> tests;
    return tests;
}

struct Registrar
{
    Registrar(const char* name, void (*func)())
    {
        registry().push_back({name, func});
    }
};

class AssertionFailure : public std::runtime_error
{
public:
    explicit AssertionFailure(const std::string& message)
    :
        std::runtime_error(message)
    {}
};

template<class Left, class Right>
[[noreturn]] inline void failEquality
(
    const char* file,
    int line,
    const char* lhsExpr,
    const char* rhsExpr,
    const Left& lhs,
    const Right& rhs
)
{
    std::ostringstream os;
    os << file << ":" << line << ": CHECK_EQ(" << lhsExpr << ", " << rhsExpr
       << ") failed with lhs=" << lhs << " rhs=" << rhs;
    throw AssertionFailure(os.str());
}

[[noreturn]] inline void failCheck
(
    const char* file,
    int line,
    const char* expression
)
{
    std::ostringstream os;
    os << file << ":" << line << ": CHECK(" << expression << ") failed";
    throw AssertionFailure(os.str());
}

template<class Left, class Right>
inline void checkEqual
(
    const Left& lhs,
    const Right& rhs,
    const char* lhsExpr,
    const char* rhsExpr,
    const char* file,
    int line
)
{
    if (!(lhs == rhs))
    {
        failEquality(file, line, lhsExpr, rhsExpr, lhs, rhs);
    }
}

inline void check
(
    bool condition,
    const char* expression,
    const char* file,
    int line
)
{
    if (!condition)
    {
        failCheck(file, line, expression);
    }
}

inline int runTests()
{
    int failed = 0;
    int passed = 0;

    for (const auto& test : registry())
    {
        try
        {
            test.func();
            ++passed;
            std::cout << "[pass] " << test.name << '\n';
        }
        catch (const std::exception& err)
        {
            ++failed;
            std::cerr << "[fail] " << test.name << '\n'
                      << err.what() << '\n';
        }
        catch (...)
        {
            ++failed;
            std::cerr << "[fail] " << test.name << '\n'
                      << "Unknown exception\n";
        }
    }

    std::cout << "Executed " << (passed + failed) << " test case(s): "
              << passed << " passed, " << failed << " failed\n";

    return failed == 0 ? 0 : 1;
}

} // namespace doctest

#define DOCTEST_DETAIL_CONCAT_IMPL(lhs, rhs) lhs##rhs
#define DOCTEST_DETAIL_CONCAT(lhs, rhs) DOCTEST_DETAIL_CONCAT_IMPL(lhs, rhs)

#define TEST_CASE(name)                                                         \
    static void DOCTEST_DETAIL_CONCAT(doctest_case_, __LINE__)();               \
    static doctest::Registrar DOCTEST_DETAIL_CONCAT                             \
    (                                                                          \
        doctest_registrar_,                                                     \
        __LINE__                                                                \
    )(name, &DOCTEST_DETAIL_CONCAT(doctest_case_, __LINE__));                   \
    static void DOCTEST_DETAIL_CONCAT(doctest_case_, __LINE__)()

#define CHECK(expression)                                                       \
    doctest::check(static_cast<bool>(expression), #expression, __FILE__, __LINE__)

#define CHECK_EQ(lhs, rhs)                                                      \
    doctest::checkEqual((lhs), (rhs), #lhs, #rhs, __FILE__, __LINE__)

#endif
