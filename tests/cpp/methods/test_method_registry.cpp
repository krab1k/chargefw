#include <chargefw/charges/atomic_charges.h>
#include <chargefw/methods/method_registry.h>

#include <cassert>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace methods = chargefw::methods;

namespace {

class StubMethod final : public methods::Method {
  public:
    explicit StubMethod(std::string id) : id_{std::move(id)},
                                          metadata_{.id = id_,
                                                    .name = id_,
                                                    .full_name = id_,
                                                    .publication = std::nullopt,
                                                    .priority = 0} {}

    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override
    {
        return metadata_;
    }

    [[nodiscard]] auto requirements() const noexcept -> methods::MethodRequirements override
    {
        return {};
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override
    {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> chargefw::charges::AtomicCharges override
    {
        return chargefw::charges::AtomicCharges{std::vector(input.prepared_molecule.molecule().atom_count(), 0.0)};
    }

  private:
    std::string id_;
    methods::MethodMetadata metadata_;
};

} // namespace

auto main() -> int
{
    auto methods_for_registry = std::vector<std::unique_ptr<methods::Method>>{};
    methods_for_registry.push_back(std::make_unique<StubMethod>("zeta"));
    methods_for_registry.push_back(std::make_unique<StubMethod>("alpha"));

    const methods::MethodRegistry registry{std::move(methods_for_registry)};

    assert(registry.methods().size() == 2);
    assert(registry.find("alpha") != nullptr);
    assert(registry.find("zeta") != nullptr);
    assert(registry.find("missing") == nullptr);
    assert((registry.names() == std::vector<std::string>{"alpha", "zeta"}));

    bool rejected_null_method = false;

    try {
        auto invalid_methods = std::vector<std::unique_ptr<methods::Method>>{};
        invalid_methods.push_back(std::unique_ptr<methods::Method>{});

        [[maybe_unused]] const methods::MethodRegistry invalid{std::move(invalid_methods)};
    } catch (const std::invalid_argument&) {
        rejected_null_method = true;
    }

    assert(rejected_null_method);

    bool rejected_empty_id = false;

    try {
        auto invalid_methods = std::vector<std::unique_ptr<methods::Method>>{};
        invalid_methods.push_back(std::make_unique<StubMethod>(""));

        [[maybe_unused]] const methods::MethodRegistry invalid{std::move(invalid_methods)};
    } catch (const std::invalid_argument&) {
        rejected_empty_id = true;
    }

    assert(rejected_empty_id);

    bool rejected_duplicate_id = false;

    try {
        auto invalid_methods = std::vector<std::unique_ptr<methods::Method>>{};
        invalid_methods.push_back(std::make_unique<StubMethod>("same"));
        invalid_methods.push_back(std::make_unique<StubMethod>("same"));

        [[maybe_unused]] const methods::MethodRegistry invalid{std::move(invalid_methods)};
    } catch (const std::invalid_argument&) {
        rejected_duplicate_id = true;
    }

    assert(rejected_duplicate_id);

    return 0;
}
