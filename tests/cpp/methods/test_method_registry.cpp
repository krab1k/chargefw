#include <chargefw/charges/atomic_charges.h>
#include <chargefw/methods/method_registry.h>

#include <memory>
#include <optional>
#include <snitch/snitch.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace methods = chargefw::methods;

namespace {

class StubMethod final : public methods::Method {
  public:
    explicit StubMethod(std::string id)
        : id_{std::move(id)}, metadata_{.id = id_,
                                        .name = id_,
                                        .full_name = id_,
                                        .publication = std::nullopt,
                                        .priority = 0} {}

    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        return metadata_;
    }

    [[nodiscard]] auto requirements() const noexcept -> methods::MethodRequirements override {
        return {};
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> chargefw::charges::AtomicCharges override {
        return chargefw::charges::AtomicCharges{std::vector(input.molecule().atom_count(), 0.0)};
    }

  private:
    std::string id_;
    methods::MethodMetadata metadata_;
};

} // namespace

TEST_CASE("method registry stores, finds, and rejects invalid methods",
          "[methods][method-registry]") {
    auto methods_for_registry = std::vector<std::unique_ptr<methods::Method>>{};
    methods_for_registry.push_back(std::make_unique<StubMethod>("zeta"));
    methods_for_registry.push_back(std::make_unique<StubMethod>("alpha"));

    const methods::MethodRegistry registry{std::move(methods_for_registry)};

    CHECK(registry.methods().size() == 2);
    CHECK(registry.find("alpha") != nullptr);
    CHECK(registry.find("zeta") != nullptr);
    CHECK(registry.find("missing") == nullptr);
    CHECK((registry.names() == std::vector<std::string>{"alpha", "zeta"}));

    const auto null_method = [] {
        auto invalid_methods = std::vector<std::unique_ptr<methods::Method>>{};
        invalid_methods.push_back(std::unique_ptr<methods::Method>{});
        [[maybe_unused]] const methods::MethodRegistry invalid{std::move(invalid_methods)};
    };
    CHECK_THROWS_AS(null_method(), std::invalid_argument);

    const auto empty_id = [] {
        auto invalid_methods = std::vector<std::unique_ptr<methods::Method>>{};
        invalid_methods.push_back(std::make_unique<StubMethod>(""));
        [[maybe_unused]] const methods::MethodRegistry invalid{std::move(invalid_methods)};
    };
    CHECK_THROWS_AS(empty_id(), std::invalid_argument);

    const auto duplicate_id = [] {
        auto invalid_methods = std::vector<std::unique_ptr<methods::Method>>{};
        invalid_methods.push_back(std::make_unique<StubMethod>("same"));
        invalid_methods.push_back(std::make_unique<StubMethod>("same"));
        [[maybe_unused]] const methods::MethodRegistry invalid{std::move(invalid_methods)};
    };
    CHECK_THROWS_AS(duplicate_id(), std::invalid_argument);
}