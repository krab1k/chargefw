#pragma once

#include <chargefw/calculation/observer.h>

namespace chargefw::calculation::detail {

// Observers are optional and purely observational. Their failures must never alter calculation
// control flow, particularly while reporting a terminal event during exception unwinding.
inline auto report_progress(const CalculationObserver& observer,
                            const CalculationProgress& progress) noexcept -> void {
    try {
        observer.on_progress(progress);
    } catch (...) {
    }
}

} // namespace chargefw::calculation::detail
