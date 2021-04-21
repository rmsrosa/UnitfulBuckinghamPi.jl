
module UnitfulBuckinghamPi

using Unitful
using LinearAlgebra

export @setparameters, @addparameters, pi_groups

"""
    param_simbols = Vector{Symbol}()

Vector containing the list of symbols of the registered parameters.
"""
param_symbols = Vector{Symbol}()

"""
    param_values = Vector{Union{Unitful.Quantity, Unitful.FreeUnits, Unitful.Dimensions, Number}}()

Vector containing the list of values of the registered parameters.
"""
param_values = Vector{Union{Unitful.Quantity, Unitful.FreeUnits, Unitful.Dimensions, Number}}()

"""
    clearregister()

Clear the parameter register.
"""
function clearregister()
    intersect!(param_symbols, Vector{Symbol}())
    intersect!(param_values, Vector{Union{Unitful.Quantity, Unitful.FreeUnits, Unitful.Dimensions, Number}}())
end

"""
    addparameter(x, y)

Add parameter with symbol `x` and value `y` to the register.
"""
function addparameter(x,y)
    y isa Union{Unitful.Quantity, Unitful.FreeUnits, Unitful.Dimensions, Number} || 
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
    rationalnullspace(mat)

Return a matrix whose columns span the null space of the matrix `mat`.

If `mat` has full rank, the returned matrix has zero columns.

The algorithm uses the LU decomposition in LinearAlgebra since
this decomposition retains the Rational eltype of the matrix.
"""
function rationalnullspace(mat)
    function ivec(n,i)
        ohv = fill(0//1, n)
        ohv[i] = 1//1
        return ohv
    end
    mat_lu = lu(mat, check=false) 
    mat_nrows, mat_ncols = size(mat)
    mat_rank = rank(mat)
    mat_null = fill(0//1, mat_ncols, mat_ncols - mat_rank)
    for j in 1:mat_ncols-mat_rank
        mat_null[1:mat_rank,j] = -(mat_lu.U[:,1:mat_rank] \ (mat_lu.U * ivec(mat_ncols, mat_rank+j)))
        mat_null[mat_rank+j,j] = 1//1
    end
    return mat_null
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
 :(ℓ ^ (-1 // 2) * g ^ (1 // 2) * τ ^ (1 // 1))
 :(θ ^ (1 // 1))

julia> pi_groups(:Expr)
2-element Vector{Expr}:
 :(ℓ ^ (-1 // 2) * g ^ (1 // 2) * τ ^ (1 // 1))
 :(θ ^ (1 // 1))

julia> pi_groups(:String)
2-element Vector{String}:
 "ℓ^(-1//2)*g^(1//2)*τ^(1//1)"
 "θ^(1//1)"
```
"""
function pi_groups(type::Symbol = :Expr)
    pdmat = parameterdimensionmatrix()
    pdmat_null = rationalnullspace(pdmat)
    groups = [join(["$p^($a)" for (p,a) in zip(param_symbols, col) if !iszero(a)], "*") for col in eachcol(pdmat_null)]
    type == :String && return groups
    type == :Expr && return [Meta.parse(group) for group in groups]
    throw(
        ArgumentError(
            "Type $type not implemented. Choose between `:Expr` and `:String`"
        )
    )
end

end # module
