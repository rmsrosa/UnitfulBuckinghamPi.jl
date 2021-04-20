using Unitful
using Test
using UnitfulBuckinghamPi

@testset "Test" begin
    # Define Unitful paramaters, which can be quantities, units or dimensions
    ℓ = u"m"
    g = 9.8u"m/s^2"
    m = u"g"
    τ = u"𝐓"
    θ = u"NoDims"
    # Set parameters
    @setparameters ℓ g m τ θ
    @test UnitfulBuckinghamPi._ubp_pars == [:ℓ, :g, :m, :τ, :θ]
    @test UnitfulBuckinghamPi._ubp_vals == [u"m", 9.8u"m/s^2", u"g", u"𝐓", u"NoDims"]
    # Check parameters
    Π = UnitfulBuckinghamPi.pi_groups_str()
    @test Π[1] == "ℓ^(-1//2)*g^(1//2)*τ^(1//1)"
    @test Π[2] == "θ^(1//1)"

    τ = u"s"
    @setparameters ℓ g m τ θ
    @test UnitfulBuckinghamPi._ubp_vals == [u"m", 9.8u"m/s^2", u"g", u"s", u"NoDims"]
    Π = UnitfulBuckinghamPi.pi_groups()
    @test Π[1] == :(ℓ ^ (-1 // 2) * g ^ (1 // 2) * τ ^ (1 // 1))
    @test Π[2] == :(θ ^ (1 // 1))
    #@test eval(Π[1]) ≈ 3.1304951684997055
    #@test eval(Π[2]) == NoDims

    v = u"m/s"
    @addparameters v
    @test UnitfulBuckinghamPi._ubp_pars == [:ℓ, :g, :m, :τ, :θ, :v]
    @test UnitfulBuckinghamPi._ubp_vals == [u"m", 9.8u"m/s^2", u"g", u"s", u"NoDims", u"m/s"]
end
