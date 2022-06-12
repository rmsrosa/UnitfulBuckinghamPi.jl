using Unitful
using Test
using UnitfulBuckinghamPi
using LinearAlgebra

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

u = u"m/s"
ρ = u"g/m^3"
μ = u"g/m/s"

@testset "UnitfulBuckinghamPi" begin
    # Set and check parameters
    @setparameters ℓ g m T θ
    @test UnitfulBuckinghamPi.param_symbols == [:ℓ, :g, :m, :T, :θ]
    @test UnitfulBuckinghamPi.param_values == [u"m", 9.8u"m/s^2", u"g", u"𝐓", u"NoDims"]

    # Check adimensional groups as String eltype
    Π = pi_groups(:String)
    @test Π[1] == "g^(1//2)*ℓ^(-1//2)*T^(1//1)"
    @test Π[2] == "θ^(1//1)"

    # Set and check parameters
    @setparameters ℓ g m τ θ
    @test UnitfulBuckinghamPi.param_values == 
        [u"m", 9.8u"m/s^2", u"g", u"s", u"NoDims"]
    # Check adimensional groups as Expr eltype
    Π = pi_groups()
    @test Π[1] == :(g ^ (1 // 2) * ℓ ^ (-1 // 2) * τ ^ (1 // 1))
    @test Π[2] == :(θ ^ (1 // 1))
    # Test evaluating expressions
    @test eval(Π[1]) ≈ 3.1304951684997055
    @test eval(Π[2]) == NoDims

    # Check add parameter and check results
    @addparameters v
    @test UnitfulBuckinghamPi.param_symbols == [:ℓ, :g, :m, :τ, :θ, :v]
    @test UnitfulBuckinghamPi.param_values == 
        [u"m", 9.8u"m/s^2", u"g", u"s", u"NoDims", u"m/s"]
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

    # Check singularity in the LU decomposition
    @setparameters u ρ μ ℓ
    @test pi_groups() == 
        [:(ρ ^ (1 // 1) * u ^ (1 // 1) * μ ^ (-1 // 1) * ℓ ^ (1 // 1))]
    @setparameters ℓ ρ μ u
    @test pi_groups() == 
        [:(ρ ^ (1 // 1) * μ ^ (-1 // 1) * ℓ ^ (1 // 1) * u ^ (1 // 1))]

    # Test errors
    @setparameters ℓ g m τ θ
    @test_throws ArgumentError pi_groups(:NotImplemented)
    @test_throws MethodError @setparameters 1
    @test_throws ArgumentError @setparameters s
end

@testset "LU with full pivoting" begin
    # Example from LinearAlgebra.lu
    A = [4 3; 6 3]
    F = UnitfulBuckinghamPi.lu_pq(A)
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test L ≈ [1.0 0.0; 0.666667 1.0] (atol = 0.00001)
    @test (F.L == L) && (F.U == U) && (F.p == p) && (F.q == q)
    @test L * U == A[p, q]

    # Singular rectangular wide matrix for which regular `LinearAlgebra.lu` fails
    A = reshape(collect(1:12),3,4)
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test L * U == A[p, q]
    @test rank(A) == rank(U) == 2

    # Same but with Rational{Int} eltype
    A = convert.(Rational, reshape(collect(1:12),3,4))
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test L * U == A[p, q]
    @test (L isa Matrix{Rational{Int64}}) && (U isa Matrix{Rational{Int64}})

    # Tall matrix version
    A = convert.(Rational, reshape(collect(1:12),4,3))
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test L * U == A[p, q]
    @test U == [
        12//1   4//1   8//1;
        0//1  -2//1  -1//1;
        0//1   0//1   0//1;
        0//1   0//1   0//1;
    ]

    # Complex{Int} version
    A = complex(reshape(collect(1:12),3,4))
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test L * U == A[p, q]
    @test (L isa Matrix{Complex{Float64}}) && (U isa Matrix{Complex{Float64}})

    # Complex{Float} version
    A = complex(reshape(collect(1.0:12.0),3,4))
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test L * U == A[p, q]
    @test (L isa Matrix{Complex{Float64}}) && (U isa Matrix{Complex{Float64}})
    
    # SVD values very far apart, for which `LinearAlgebra.rank` fails
    A = [1 1; typemax(Int64)//2 1]
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test L * U == A[p, q]

    # Singular case associated with the second Reynolds number example above
    A = convert.(Rational,[-1 0 -2 -1; 0 1 1 1; 1 -3 -1 -1])
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test U == [
        -3//1 -1//1 -1//1 1//1;
        0//1  -2//1 -1//1 -1//1;
        0//1  0//1  1//3  0//1
    ]

    # Complex{Int} version, with round-offs
    A = convert.(Complex,[-1 0 -2 -1; 0 1 1 1; 1 -3 -1 -1])
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test (L isa Matrix{Complex{Float64}}) && (U isa Matrix{Complex{Float64}})
    @test U ≈ [
        -3.0+0.0im -1.0+0.0im  -1.0+0.0im 1.0+0.0im;
        0.0+0.0im  -2.0+0.0im -1.0+0.0im -1.0+0.0im;
        0.0+0.0im   0.0+0.0im  0.333333+0.0im  0.0+0.0im] (atol = 0.00001)

    # Random wide matrix
    A = rand(7,11)
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test L * U ≈ A[p,q]

    # Random tall matrix
    A = rand(11,7)
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(A)
    @test L * U ≈ A[p,q]

    # Rational{BigInt} matrix that overflows with Rational{Int}
    A = reshape(1 .//collect(10001:10012), 4, 3)
    @test_throws OverflowError UnitfulBuckinghamPi.lu_pq(A)
    L, U, p, q = UnitfulBuckinghamPi.lu_pq(convert.(Rational{BigInt}, A))
    @test L * U == A[p, q]
end
