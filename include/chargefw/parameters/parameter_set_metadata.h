#pragma once

#include <cstdint>
#include <string>

namespace chargefw::parameters {

struct ParameterSetMetadata {
    std::string id;
    std::string method_id;
    std::string name;
    std::string publication;
    std::string notes;
    std::uint16_t priority = 0;
};

} // namespace chargefw::parameters