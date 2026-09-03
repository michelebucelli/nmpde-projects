#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <utility>

// Minimal check-and-count harness. Silent on success, one line per failure,
// and a final summary; the exit code is what ctest looks at.
class TestSuite {
public:
  explicit TestSuite(std::string name) : name(std::move(name)) {}

  auto check(const bool condition, const std::string &what) -> void {
    if (condition) {
      ++n_passed;
    } else {
      ++n_failed;
      std::cout << name << ": FAILED " << what << std::endl;
    }
  }

  auto close(const double actual, const double expected, const double tolerance,
             const std::string &what) -> void {
    const bool ok = std::abs(actual - expected) <=
                    tolerance * std::max(1.0, std::abs(expected));
    if (!ok) {
      std::cout << name << ": FAILED " << what << " (got " << actual
                << ", expected " << expected << ")" << std::endl;
      ++n_failed;
    } else {
      ++n_passed;
    }
  }

  auto report() const -> int {
    std::cout << name << ": " << n_passed << "/" << n_passed + n_failed
              << " checks passed" << std::endl;
    return n_failed == 0 ? 0 : 1;
  }

private:
  std::string name;
  unsigned int n_passed = 0;
  unsigned int n_failed = 0;
};
