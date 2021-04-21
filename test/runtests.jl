using Unitful
using Test
using UnitfulBuckinghamPi

# Define Unitful paramaters, which can be quantities, units or dimensions
ℓ = u"m"
g = 9.8u"m/s^2"
m = u"g"
T = u"𝐓"
τ = u"s"
θ = u"NoDims"
v = u"m/s"
α = 2
s = "blah"

@testset "Test" begin
    # Set and check parameters
    @setparameters ℓ g m T θ
    @test UnitfulBuckinghamPi.param_symbols == [:ℓ, :g, :m, :T, :θ]
    @test UnitfulBuckinghamPi.param_values == [u"m", 9.8u"m/s^2", u"g", u"𝐓", u"NoDims"]

    # Check adimensional groups as String eltype
    Π = pi_groups(:String)
    @test Π[1] == "ℓ^(-1//2)*g^(1//2)*T^(1//1)"
    @test Π[2] == "θ^(1//1)"

    # Set and check parameters
    @setparameters ℓ g m τ θ
    @test UnitfulBuckinghamPi.param_values == [u"m", 9.8u"m/s^2", u"g", u"s", u"NoDims"]
    # Check adimensional groups as Expr eltype
    Π = pi_groups()
    @test Π[1] == :(ℓ ^ (-1 // 2) * g ^ (1 // 2) * τ ^ (1 // 1))
    @test Π[2] == :(θ ^ (1 // 1))
    # Test evaluating expressions
    @test eval(Π[1]) ≈ 3.1304951684997055
    @test eval(Π[2]) == NoDims

    # Check add parameter and check results
    @addparameters v
    @test UnitfulBuckinghamPi.param_symbols == [:ℓ, :g, :m, :τ, :θ, :v]
    @test UnitfulBuckinghamPi.param_values == [u"m", 9.8u"m/s^2", u"g", u"s", u"NoDims", u"m/s"]
    Π = pi_groups(:Expr)
    @test length(Π) == 3

    # Test parameter of type Number
    @addparameters α
    @test length(UnitfulBuckinghamPi.param_symbols) == 7
    Π = pi_groups()
    @test length(Π) == 4

    # Avoid adding duplicates
    @addparameters α
    @test length(UnitfulBuckinghamPi.param_symbols) == 7

    # Check setting no parameters for an empty list
    @setparameters
    Π = pi_groups()
    @test size(Π) == (0,)

    # Test errors
    @test_throws ArgumentError pi_groups(:NotImplemented)
    @test_throws MethodError @setparameters 1
    @test_throws ArgumentError @setparameters s
end
