module UnitfulBuckinghamPi

using Unitful
using LinearAlgebra

export @setparameters, @addparameters, pi_groups

"""
    ParamTypes

Type `Union{...}` with the allowed types for the parameters: `Unitful.Quantity`,
`Unitful.Units`, `Unitful.Dimensions`, and `Number`.
"""
const ParamTypes = 
    Union{Unitful.Quantity, Unitful.FreeUnits, Unitful.Dimensions, Number}

"""
    param_symbols = Vector{Symbol}()

Vector containing the list of symbols of the registered parameters.
"""
param_symbols = Vector{Symbol}()

"""
    param_values = Vector{ParamTypes}()

Vector containing the list of values of the registered parameters.

See [`ParamTypes`](@ref) for the allowed parameter types.
"""
param_values = Vector{ParamTypes}()

"""
    clearregister()

Clear the parameter register.
"""
function clearregister()
    empty!(param_symbols)
    empty!(param_values)
end

"""
    addparameter(x, y)

Add parameter with symbol `x` and value `y` to the register.
"""
function addparameter(x,y)
    y isa ParamTypes || 
        throw(ArgumentError(
            "Parameter should be either a Unitful.Quantity, Unitful.FreeUnits,"
            * "Unitful.Dimensions or Number"))
    
    if x ∉ UnitfulBuckinghamPi.param_symbols
        push!(param_symbols, x)
        push!(param_values, y)
    end
end

"""
    @setparameters args...

Register a list of parameters.

The parameters should be of any of following types:
    * Unitful.Quantity
    * Unitful.FreeUnits
    * Unitful.Dimensions
    * Number

# Examples

```jldoctest
julia> ℓ = u"m"; g = 9.8u"m/s^2"; m = u"g"; T = u"𝐓"; θ = u"NoDims";

julia> @setparameters ℓ g m τ θ
[ Info: Parameter(s) registered:
[ Info:  ℓ = m
[ Info:  g = 9.8 m s⁻²
[ Info:  m = g
[ Info:  τ = s
[ Info:  θ = NoDims
```
"""
macro setparameters(args...)    
    quote
        $UnitfulBuckinghamPi.clearregister()
        $([:($UnitfulBuckinghamPi.addparameter($(Meta.quot(arg)), $(esc(arg))))
            for arg in args]...)
        $UnitfulBuckinghamPi.displayregister()
    end
end

"""
    @addparameters args...

Add a list of parameters to the register.

The parameters should be of any of following types:
    * Unitful.Quantity
    * Unitful.FreeUnits
    * Unitful.Dimensions
    * Number

# Examples

```jldoctest
julia> ℓ = u"m"; g = 9.8u"m/s^2"; m = u"g"; T = u"𝐓"; θ = u"NoDims";

julia> @setparameters ℓ g m
[ Info: Parameter(s) registered:
[ Info:  ℓ = m
[ Info:  g = 9.8 m s⁻²
[ Info:  m = g

julia> @addparameters T θ
[ Info: Parameter(s) registered:
[ Info:  ℓ = m
[ Info:  g = 9.8 m s⁻²
[ Info:  m = g
[ Info:  T = 𝐓
[ Info:  θ = NoDims
```
"""
macro addparameters(args...)
    quote
        $([:($UnitfulBuckinghamPi.addparameter($(Meta.quot(arg)), $(esc(arg)))) for arg in args]...)
        $UnitfulBuckinghamPi.displayregister()
    end
end

"""
    displayregister()

Display info with the list of registered parameters.
"""
function displayregister()
    @info "Parameter(s) registered:"
    for (n, v) in zip(param_symbols, param_values)
        @info " $n = $v"
    end
end

"""
    param_dimensions(x)

Return name, abbreviation and power of the dimensions of the given argument `x`.

The argument `x` can be a `Uniftul.Quantity`, a `Unitful.FreeUnits`,
a `Unitful.Dimensions` or a `Number`.
"""
function param_dimensions(::Unitful.Dimensions{D}) where {D}
    return [(Unitful.name(z), Unitful.abbr(z), Unitful.power(z)) for z in D]
end
param_dimensions(u::Unitful.FreeUnits) = param_dimensions(dimension(u))
param_dimensions(q::Unitful.Quantity) = param_dimensions(dimension(q))
param_dimensions(x::Number) = param_dimensions(dimension(x))

"""
    parameterdimensionmatrix()

Build and return the parameter-to-dimension matrix.
"""
function parameterdimensionmatrix()
    dims = Set(dims[2] for p in param_values for dims in param_dimensions(p))
    dims_dict = Dict(v => i for (i,v) in enumerate(dims))
    pdmat = fill(0//1, length(dims), length(param_values))
    for (j,p) in enumerate(param_values)
        for d in param_dimensions(p)
            pdmat[dims_dict[d[2]],j] = d[3]
        end
    end
    return pdmat
end

"""
    lu_pq(A::Matrix)

Compute the LU factorization of `A` with full pivoting.

Output is a `NamedTuple`, say `F=lu_pq(A)`, whose fields are
* `F.L`: `L` (lower triangular) part of the factorization
* `F.U`: `U` (upper triangular) part of the factorization
* `p`: row permutation `Vector`
* `q`: column permutation `Vector`

Factorization yields the identity `F.L * F.U == A[F.p, F.q]`.

More information about LU factorization with full pivoting,
see [^GolubVanLoan1996], [^TrefethenBauIII1997].

The choice of the tolerance, however, is beyond me, at the moment,
but see [^issue8859] and [^pinv]

# Examples

```jldoctest
julia> A = [4 3; 6 3]
2×2 Matrix{Int64}:
 4  3
 6  3

julia> F = lu_pq(A)
(L = [1.0 0.0; 0.6666666666666666 1.0], U = [6.0 3.0; 0.0 1.0], p = [2, 1], q = [1, 2])

julia> F.L * F.U == A[F.p, F.q]
true

julia> L, U, p, q = lu_pq(A)
(L = [1.0 0.0; 0.6666666666666666 1.0], U = [6.0 3.0; 0.0 1.0], p = [2, 1], q = [1, 2])

julia> L * U == A[p, q]
true

julia> A = reshape(collect(1:12),3,4)
3×4 Matrix{Int64}:
 1  4  7  10
 2  5  8  11
 3  6  9  12

julia> L, U, p, q = lu_pq(A)
(L = [1.0 0.0 0.0; 0.8333333333333334 1.0 0.0; 0.9166666666666666 0.5 1.0], U = [12.0 3.0 9.0 6.0; 0.0 -1.5 -0.5 -1.0; 0.0 0.0 0.0 0.0], p = [3, 1, 2], q = [4, 1, 3, 2])

julia> L * U == A[p, q]
true

julia> rank(A) == rank(U) == 2
true

julia> A = convert.(Rational, reshape(collect(1:12),3,4))
3×4 Matrix{Rational{Int64}}:
 1//1  4//1  7//1  10//1
 2//1  5//1  8//1  11//1
 3//1  6//1  9//1  12//1

julia> L, U, p, q = lu_pq(A);

julia> U
3×4 Matrix{Rational{Int64}}:
 12//1   3//1   9//1   6//1
  0//1  -3//2  -1//2  -1//1
  0//1   0//1   0//1   0//1

julia> A = complex(reshape(collect(1.0:12.0),3,4))
3×4 Matrix{ComplexF64}:
 1.0+0.0im  4.0+0.0im  7.0+0.0im  10.0+0.0im
 2.0+0.0im  5.0+0.0im  8.0+0.0im  11.0+0.0im
 3.0+0.0im  6.0+0.0im  9.0+0.0im  12.0+0.0im

julia> L, U, p, q = lu_pq(A);

julia> U
3×4 Matrix{ComplexF64}:
 12.0+0.0im   3.0+0.0im           9.0+0.0im   6.0+0.0im
  0.0+0.0im  -1.5+0.0im          -0.5+0.0im  -1.0+0.0im
  0.0+0.0im   0.0+0.0im  -4.44089e-16+0.0im   0.0+0.0im

julia> L * U == A[p, q]
true

julia> A = [1 1; typemax(Int64)//2 1]
2×2 Matrix{Rational{Int64}}:
                   1//1  1//1
 9223372036854775807//2  1//1

julia> L, U, p, q = lu_pq(A);

julia> U
2×2 Matrix{Rational{Int64}}:
 9223372036854775807//2                    1//1
                   0//1  9223372036854775805//9223372036854775807

julia> L * U == A[p, q]
true

julia> A = reshape( 1 .//collect(BigInt, 10001:10012), 4, 3)
4×3 Matrix{Rational{BigInt}}:
 1//10001  1//10005  1//10009
 1//10002  1//10006  1//10010
 1//10003  1//10007  1//10011
 1//10004  1//10008  1//10012

julia> L, U, p, q = lu_pq(A);

julia> L * U == A[p, q]
true

```

# References:

[^GolubVanLoan1996]: Gene H. Golub, Charles F. Van Loan, "Matrix Computations", Johns Hopkins 
Studies in Mathematical Sciences, 3rd Edition, 1996.

[^TrefethenBauIII1997]: Lloyd N. Trefethen, David Bau III, "Numerical Linear Algebra",
First Edition, SIAM, 1997.

[^issue8859]: Issue 8859, "Fix least squares", [https://github.com/JuliaLang/julia/pull/8859](https://github.com/JuliaLang/julia/pull/8859)

[^pinv]: [LinearAlgebra.pinv](https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/dense.jl#L1362)
"""
function lu_pq(A::Matrix{T}) where T <: Number
    tol = T <: Rational ? 0//1 : min(size(A)...)*eps(real(float(one(T))))
    U = copy(A)
    (n,m) = size(U)
    L = diagm(ones(eltype(A),n))
    p = collect(1:n)
    q = collect(1:m)
    for k =1:min(n,m)
        i, j = k-1 .+ Tuple(argmax(abs.(U[k:end,k:end])))
        abs(U[i,j]) > tol || break
        if i > k
            p[k], p[i] = p[i], p[k]
            U[k,:], U[i,:] = U[i,:], U[k,:]
            L[k,1:k-1], L[i,1:k-1] = L[i,1:k-1], L[k,1:k-1]
        end
        if j > k
            q[k], q[j] = q[j], q[k]
            U[:,k], U[:,j] = U[:,j], U[:,k]
        end
        τ = U[k+1:end,k] / U[k,k]
        U[k+1:end,k:end] -=  τ * U[k,k:end]'
        L[k+1:end,k] = τ
    end
    return (L=L, U=U, p=p, q=q)
end
lu_pq(A::Matrix{T}) where T <: Integer = lu_pq(float(A))
lu_pq(A::Matrix{Complex{T}}) where T <: Integer = lu_pq(convert.(float(Complex{Int}),A))

"""
    lu_nullspace(mat)

Return a matrix whose columns span the null space of the matrix `mat`.

If `mat` has full rank, the returned matrix has zero columns.

The algorithm uses the implementation [`lu_pq`](@ref) of the LU decomposition
with full pivoting. This decomposition not only retains the Rational
eltype of the matrix, but also handles singular matrix.
"""
function lu_nullspace(mat)
    function ivec(n,i)
        ohv = fill(zero(eltype(mat)), n)
        ohv[i] = one(eltype(mat))
        return ohv
    end
    mat_lu = lu_pq(mat)
    mat_nrows, mat_ncols = size(mat)
    mat_rank = rank(mat)
    mat_null = fill(zero(eltype(mat)), mat_ncols, mat_ncols - mat_rank)
    for j in 1:mat_ncols-mat_rank
        mat_null[1:mat_rank,j] = -(mat_lu.U[:,1:mat_rank] \ (mat_lu.U * ivec(mat_ncols, mat_rank+j)))
        mat_null[mat_rank+j,j] = one(eltype(mat))
    end
    return mat_null, mat_lu.q
end

"""
    pi_groups(type::Symbol = :Expr)

Return a vector with the Pi groups associated with the registered parameters.

The returned vector can be of one of the following eltypes depending
on the value of the argument `type`:
    * eltype `Expr`, if `type` is set to `Expr`
    * eltype `String`, if `type` is set to `:String`

# Examples

```jldoctest
julia> ℓ = u"m"; g = 9.8u"m/s^2"; m = u"g"; τ = u"s"; θ = u"NoDims";

julia> ℓ = u"m"; g = 9.8u"m/s^2"; m = u"g"; τ = u"s"; θ = u"NoDims";

julia> @setparameters ℓ g m τ θ
[ Info: Parameter(s) registered:
[ Info:  ℓ = m
[ Info:  g = 9.8 m s⁻²
[ Info:  m = g
[ Info:  τ = s
[ Info:  θ = NoDims

julia> pi_groups()
2-element Vector{Expr}:
 :(g ^ (1 // 2) * ℓ ^ (-1 // 2) * τ ^ (1 // 1))
 :(θ ^ (1 // 1))

julia> pi_groups(:Expr)
 2-element Vector{Expr}:
  :(g ^ (1 // 2) * ℓ ^ (-1 // 2) * τ ^ (1 // 1))
  :(θ ^ (1 // 1))

julia> pi_groups(:String)
2-element Vector{String}:
  "g^(1//2)*ℓ^(-1//2)*τ^(1//1)"
  "θ^(1//1)"
```
"""
function pi_groups(type::Symbol = :Expr)
    pdmat = parameterdimensionmatrix()
    pdmat_null, q = lu_nullspace(pdmat)
    groups = [join(["$p^($a)" for (p,a) in zip(param_symbols[q], col) if !iszero(a)], "*") for col in eachcol(pdmat_null)]
    type == :String && return groups
    type == :Expr && return [Meta.parse(group) for group in groups]
    throw(
        ArgumentError(
            "Type $type not implemented. Choose between `:Expr` and `:String`"
        )
    )
end

end # module
