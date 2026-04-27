#pragma once

#include <string>

namespace chargefw::parameters {

struct ParameterSetMetadata {
    std::string id;
    std::string method_id;
    std::string name;
    std::string publication;
    std::string notes;
};

} // namespace chargefw::parameters