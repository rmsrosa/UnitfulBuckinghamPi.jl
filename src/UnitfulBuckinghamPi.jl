
module UnitfulBuckinghamPi

using Unitful
using LinearAlgebra

export @setparameters, @addparameters, pi_groups

param_symbols = Vector{Symbol}()
param_values = Vector{Any}()
#param_values = Vector{Unitful.Quantity, Unitful.Unitlike, Unitul.Dimensions{()}}()

"""
    clearregister()

Clear the parameter register.
"""
function clearregister()
    intersect!(param_symbols, Vector{Symbol}())
    intersect!(param_values, Vector{Any}())
end

"""
    add2register(x, y)

Add parameter symbol `x` and value `y` to the register.
"""
function add2register(x,y)
    push!(param_symbols, x)
    push!(param_values, y)
end


"""
    @setparameters args...

Register a list of parameters.
"""
macro setparameters(args...)
    quote
        $UnitfulBuckinghamPi.clearregister()
        $([:($UnitfulBuckinghamPi.add2register($(Meta.quot(arg)), $(esc(arg)))) for arg in args]...)
        $UnitfulBuckinghamPi.register()
    end
end

"""
    @addparameters args...

Add a list of parameters to the register.
"""
macro addparameters(args...)
    quote
        $([:($UnitfulBuckinghamPi.add2register($(Meta.quot(arg)), $(esc(arg)))) for arg in args]...)
        $UnitfulBuckinghamPi.register()
    end
end

"""
    register()

Display info with the list of registered parameters
"""
function register()
    @info "Parameter(s) registered:"
    for (n, v) in zip(param_symbols, param_values)
        @info " $n = $v"
    end
end

"""
    get_dims(x)

Return name, abbreviation and power of the given argument.

The argument `x` can be a Uniftul.Quantity, a Unitful.FreeUnits, a Unitful.Dimensions or a nondimenional type.
"""
function get_dims(::Unitful.Dimensions{D}) where {D}
    return [(Unitful.name(z), Unitful.abbr(z), Unitful.power(z)) for z in D]
end
get_dims(u::Unitful.FreeUnits) = get_dims(dimension(u))
get_dims(q::Unitful.Quantity) = get_dims(dimension(q))
get_dims(x::Number) = get_dims(dimension(x))

"""
    parameterdimensionmatrix()

Build and return the parameter-to-dimension matrix.
"""
function parameterdimensionmatrix()
    param_dimensions = Set(dims[2] for p in param_values for dims in get_dims(p))
    param_dims_dict = Dict(v => i for (i,v) in enumerate(param_dimensions))
    pdmat = fill(0//1, length(param_dimensions), length(param_values))
    for (j,p) in enumerate(param_values)
        for d in get_dims(p)
            pdmat[param_dims_dict[d[2]],j] = d[3]
        end
    end
    return pdmat
end

"""
    rationalnullspace(mat)

Return a matrix whose columns span the null space of the matrix `mat`.

If `mat` has full rank, the returned matrix has zero columns.

Algorithm uses the LU decomposition in LinearAlgebra since this decomposition
retains the Rational eltype of the matrix.
"""
function rationalnullspace(mat)
    function ivec(n,i)
        ohv = fill(0//1, n)
        ohv[i] = 1//1
        return ohv
    end
    mat_lu = lu(mat)
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
    pi_groups_str()

Return a vector of strings with the pi groups associated with the registered parameters.
"""
function pi_groups_str()
    pdmat = parameterdimensionmatrix()
    pdmat_null = rationalnullspace(pdmat)
    return [join(["$p^($a)" for (p,a) in zip(param_symbols, col) if !iszero(a)], "*") for col in eachcol(pdmat_null)]
end

"""
    pi_groups()

Return a vector of expressions with the pi groups associated with the registered parameters.

This only works if all the parameters are dimensions or if none of them are,
since it is not possible to multiply a Unitful.Dimensions with a Unitful.FreeUnits
or with at Unitful.Quantity.
"""
function pi_groups()
    groups = pi_groups_str()
    return [Meta.parse(group) for group in groups]
end

end # module
