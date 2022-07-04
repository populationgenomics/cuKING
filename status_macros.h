#pragma once

// Unfortunately these status macros are still not part of Abseil. This is
// adapted from
// https://source.chromium.org/chromiumos/chromiumos/codesearch/+/main:src/platform2/missive/util/status_macros.h
// but allows converting from different Status implementations
// (e.g. absl::Status and arrow::Status) by providing an override for the
// ToAbslStatus converter.

#include <absl/base/optimization.h>
#include <absl/status/status.h>

namespace status_macros {

inline absl::Status ToAbslStatus(absl::Status status) {
  return std::move(status);
}

}  // namespace status_macros

#define RETURN_IF_ERROR(expr)                                                \
  do {                                                                       \
    /* Using _status below to avoid capture problems if expr is "status". */ \
    const auto _status = (expr);                                             \
    if (ABSL_PREDICT_FALSE(!_status.ok())) {                                 \
      return status_macros::ToAbslStatus(_status);                           \
    }                                                                        \
  } while (0)

// Internal helper for concatenating macro values.
#define STATUS_MACROS_CONCAT_NAME_INNER(x, y) x##y
#define STATUS_MACROS_CONCAT_NAME(x, y) STATUS_MACROS_CONCAT_NAME_INNER(x, y)

#define ASSIGN_OR_RETURN_IMPL(result, lhs, rexpr)        \
  auto result = rexpr;                                   \
  if (ABSL_PREDICT_FALSE(!result.ok())) {                \
    return status_macros::ToAbslStatus(result.status()); \
  }                                                      \
  lhs = std::move(result).ValueOrDie()

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
