# UnitfulBuckinghamPi

![Main Tests Workflow Status](https://github.com/rmsrosa/UnitfulCurrencies.jl/workflows/CI/badge.svg) [![codecov](https://codecov.io/gh/rmsrosa/UnitfulBuckinghamPi.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/rmsrosa/UnitfulBuckinghamPi.jl) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![GitHub repo size](https://img.shields.io/github/repo-size/rmsrosa/UnitfulCurrencies.jl) ![Lifecycle Experimental](https://img.shields.io/badge/lifecycle-experimental-orange) [![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package is for solving for the adimensional Pi groups (or Π groups) in a given list of parameters, according to the [Buckingham-Pi Theorem](https://en.wikipedia.org/wiki/Buckingham_π_theorem).

We use the package [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) in order to facilitate the construction of the parameters and to easily handle the dimensions associated with each parameter.

This package is inspired by a similar package written in python: [ian-r-rose/buckinghampy](https://github.com/ian-r-rose/buckinghampy).

## Example

As a simple example, let us consider the period of a simple pendulum.

We consider the *length* of the rod, the *mass* of the bob, the *acceleration of gravity*, the *angle* of the rod with the downards vertical direction, and the *period* of the swinging pendulum as the relevant parameters.

We defined these parameters as `Unitful.FreeUnits`. Except for the acceleration of gravity, which is a constant and is given a `Unitful.Quantity` value, and for the period, for which we do not associate any unit, only a dimension, just for fun.

To tell `UnitfulBuckinghamPi` that these are the parameters to consider, we use the macro `@setparameters`. Then, we find the adimensional Π groups with the function `pi_groups_str()`, which returns the groups as a vector of strings. Or we can use `pi_groups()`, which returns a vector of expressions.

```julia
julia> using Unitful

julia> using UnitfulBuckinghamPi

julia> ℓ = u"m"
m

julia> g = 9.8u"m/s^2"
9.8 m s⁻²

julia> m = u"g"
g

julia> T = u"𝐓"
𝐓

julia> θ = u"NoDims"
NoDims

julia> @setparameters ℓ g m T θ
[ Info: Parameter(s) registered:
[ Info:  ℓ = m
[ Info:  g = 9.8 m s⁻²
[ Info:  m = g
[ Info:  T = 𝐓
[ Info:  θ = NoDims

julia> Π_str = pi_groups_str()
2-element Vector{String}:
 "ℓ^(-1//2)*g^(1//2)*T^(1//1)"
 "θ^(1//1)"

julia> Π = pi_groups()
2-element Vector{Expr}:
 :(ℓ ^ (-1 // 2) * g ^ (1 // 2) * T ^ (1 // 1))
 :(θ ^ (1 // 1)) 
```

There are two adimensional groups, `Π[1]` and `Π[2]`.

One can use [korsbo/Latexify.jl](https://github.com/korsbo/Latexify.jl) to display the groups in Latex format, but be aware that Latexify doesn't properly render Rational numbers when they appear as powers of another quantity. So, one needs to replace the double backslashes with a single backslash for a proper display, like with `latexify(replace(Π_str[1], "//" => "/"))`.

With the parameters above, one cannot evaluate the adimensional group since that would amount to multiplying Unitful.FreeUnits or Unitful.Quantities like the Unitful.Dimensions parameter `T`. That i not allowed by `Unitful.jl`. One can solve that, however, by substituting `T` with a unit. Then, we can either parse each element in the vector of strings returned by `pi_groups_str()` and evaluate that or we can use `pi_groups()` to obtain directly the corresponding expressions and evaluate the expressions.

```julia
julia> τ = u"s"
s

julia> @setparameters ℓ g m τ θ
[ Info: Parameter(s) registered:
[ Info:  ℓ = m
[ Info:  g = 9.8 m s⁻²
[ Info:  m = g
[ Info:  τ = s
[ Info:  θ = NoDims

julia> Π = pi_groups()
2-element Vector{Expr}:
 :(ℓ ^ (-1 // 2) * g ^ (1 // 2) * τ ^ (1 // 1))
 :(θ ^ (1 // 1))

julia> eval(Π[1])
3.1304951684997055

julia> eval(Π[2])
NoDims
```

As expected, both are adimensional.

If there are more than one set of groups which are not independent, the solver will just pick one combination from the basis obtained for the null space.

Finally, one can add parameters to a given set of registered parameters and solve for the new set.

```julia
julia> v = u"m/s"
m s⁻¹

julia> @addparameters v
[ Info: Parameter(s) registered:
[ Info:  ℓ = m
[ Info:  g = 9.8 m s⁻²
[ Info:  m = g
[ Info:  τ = s
[ Info:  θ = NoDims
[ Info:  v = m s⁻¹

julia> pi_groups()
3-element Vector{Expr}:
 :(ℓ ^ (-1 // 2) * g ^ (1 // 2) * τ ^ (1 // 1))
 :(θ ^ (1 // 1))
 :(ℓ ^ (-1 // 2) * g ^ (-1 // 2) * v ^ (1 // 1))
```

## The internals

The [Buckingham-Pi Theorem](https://en.wikipedia.org/wiki/Buckingham_π_theorem) relies on the [Rank-nulity Theorem](https://en.wikipedia.org/wiki/Rank–nullity_theorem). A parameter-to-dimension matrix is composed, in which the columns correpond to the parameters and the rows to the collection of dimensions involved in the parameters. Each element in row i and column j correspond to the power of the dimension i in the parameter j.

The number of adimensional groups is the dimension of the kernel of the matrix. And the adimensional groups are obtained from a basis of the null space.

When the powers are integers or rational numbers, which is usually the case, it is desirable to keep the type of these parameters when composing the matrix and when finding the null space and the associated adimensional Π groups.

While [ian-r-rose/buckinghampy](https://github.com/ian-r-rose/buckinghampy) uses [SymPy](https://www.sympy.org/en/index.html) for symbolic manipulation of the powers of the parameters, to retain these types, we simply rely on the ability of the `LU` decomposition in the [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) standard package to retain the `Rational` eltype of the matrices.

Associated with that, we do not use the [`LinearAlgebra.nullspace`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.nullspace) since it is based on the `QR` decomposition, which changes the eltype to `Float64`. Instead, we compute the null space directly from the `U` part of the `LU` decomposition.

## License

This package is licensed under the [MIT license](https://opensource.org/licenses/MIT) (see file [LICENSE](LICENSE) in the root directory of the project).
