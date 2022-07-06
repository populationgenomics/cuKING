#pragma once

#include <absl/base/optimization.h>
#include <absl/status/status.h>
#include <absl/time/time.h>

namespace cuking {

// Returns ceil(a / b) for integers a, b.
template <typename T>
inline T CeilIntDiv(const T a, const T b) {
  return (a + b - 1) / b;
}

// Keeps track of time intervals.
class StopWatch {
 public:
  absl::Duration GetElapsedAndReset() {
    const auto now = absl::Now();
    const auto result = now - last_time_;
    last_time_ = now;
    return result;
  }

 private:
  absl::Time last_time_ = absl::Now();
};

inline absl::Status ToAbslStatus(absl::Status status) {
  return std::move(status);
}

// Unfortunately these status macros are still not part of Abseil. This is
// adapted from
// https://source.chromium.org/chromiumos/chromiumos/codesearch/+/main:src/platform2/missive/util/status_macros.h
// but allows converting from different Status implementations
// (e.g. absl::Status and arrow::Status) by providing an override for the
// ToAbslStatus converter.

#define RETURN_IF_ERROR(expr)                                                \
  do {                                                                       \
    /* Using _status below to avoid capture problems if expr is "status". */ \
    const auto _status = (expr);                                             \
    if (ABSL_PREDICT_FALSE(!_status.ok())) {                                 \
      return cuking::ToAbslStatus(_status);                                  \
    }                                                                        \
  } while (0)

// Internal helper for concatenating macro values.
#define STATUS_MACROS_CONCAT_NAME_INNER(x, y) x##y
#define STATUS_MACROS_CONCAT_NAME(x, y) STATUS_MACROS_CONCAT_NAME_INNER(x, y)

#define ASSIGN_OR_RETURN_IMPL(result, lhs, rexpr) \
  auto result = rexpr;                            \
  if (ABSL_PREDICT_FALSE(!result.ok())) {         \
    return cuking::ToAbslStatus(result.status()); \
  }                                               \
  lhs = *std::move(result);

// Executes an expression that returns a StatusOr, extracting its value
// into the variable defined by lhs (or returning on error).
//
// Example: Assigning to an existing value
//   ValueType value;
//   ASSIGN_OR_RETURN(value, MaybeGetValue(arg));
//
// Example: Creating and assigning variable in one line.
//   ASSIGN_OR_RETURN(ValueType value, MaybeGetValue(arg));
//   DoSomethingWithValueType(value);
//
// WARNING: ASSIGN_OR_RETURN expands into multiple statements; it cannot be used
//  in a single statement (e.g. as the body of an if statement without {})!
#define ASSIGN_OR_RETURN(lhs, rexpr) \
  ASSIGN_OR_RETURN_IMPL(             \
      STATUS_MACROS_CONCAT_NAME(_status_or_value, __COUNTER__), lhs, rexpr)

}  // namespace cuking
