# UnitfulBuckinghamPi

![Main Tests Workflow Status](https://github.com/rmsrosa/UnitfulBuckinghamPi.jl/workflows/CI/badge.svg) [![codecov](https://codecov.io/gh/rmsrosa/UnitfulBuckinghamPi.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/rmsrosa/UnitfulBuckinghamPi.jl) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![GitHub repo size](https://img.shields.io/github/repo-size/rmsrosa/UnitfulCurrencies.jl) ![Lifecycle Experimental](https://img.shields.io/badge/lifecycle-experimental-orange) [![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This package is for solving for the adimensional Pi groups (or Π groups) in a given list of parameters, according to the [Buckingham-Pi Theorem](https://en.wikipedia.org/wiki/Buckingham_π_theorem).

We use the package [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) in order to facilitate the construction of the parameters and to easily handle the dimensions associated with each parameter.

This package is inspired by a similar package written in python: [ian-r-rose/buckinghampy](https://github.com/ian-r-rose/buckinghampy).

## Examples

Let us consider a couple of examples.

### Simple pendulum

We start with the period of a simple pendulum.

The parameters taken for consideration are the *length* of the rod, the *mass* of the bob, the *acceleration of gravity*, the *angle* of the rod with respect to the downards vertical direction, and the *period* of the swinging pendulum.

We define these parameters as `Unitful.FreeUnits`. Except for the acceleration of gravity, which is a constant and is given as a `Unitful.Quantity` value, and for the period, for which we do not associate any unit, only a dimension, just for fun.

To tell `UnitfulBuckinghamPi` that these are the parameters to consider, we use the macro `@setparameters`. Then, we find the adimensional Π groups with the function `pi_groups()`, which returns the groups as a vector. It can either be a vector of strings, with `pi_groups(:String)`, or of expressions, with `pi_groups(:Expr)`, which is the default.

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

julia> Π_str = pi_groups(:String)
2-element Vector{String}:
 "ℓ^(-1//2)*g^(1//2)*T^(1//1)"
 "θ^(1//1)"

julia> Π = pi_groups(:Expr)
2-element Vector{Expr}:
 :(ℓ ^ (-1 // 2) * g ^ (1 // 2) * T ^ (1 // 1))
 :(θ ^ (1 // 1)) 
```

There are two adimensional groups, `Π[1]` and `Π[2]`.

One can use [korsbo/Latexify.jl](https://github.com/korsbo/Latexify.jl) to display the groups in Latex format, but be aware that Latexify doesn't properly render Rational numbers when they appear as powers of another quantity. So, one needs to replace the double backslashes with a single backslash for a proper display, like with `latexify(replace(Π_str[1], "//" => "/"))`. Doing so, we obtain the image

![pendulum adimensional Pi group](img/pendulum_pi_group.png)

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

Adimensional groups are not unique. Even if you find a single group, any power of it is also adimensional. If there is more than one adimensional group, then you can linearly combine them to find many others. They are associated to the span of the null space of a "parameter-to-dimension" matrix. The solver here will just pick one combination from the basis obtained for the null space.

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

### Reynolds number

Another classical example of adimensional group is the Reynolds number.

What could characterize how complicate a fluid flow is? We should certainly consider the parameters characterizing the fluid, such as *density* and *viscosity*. Then, there is the *velocity* the fluid moves. Less obvious is the *lenght scale* we consider, which we can consider as a length scale associated with the injection of energy into the system, such as the width of an obstacle, the distance between the walls of channel, the distance between the bars of a grids, and so on.

With these parameters, the only possible adimensional groups is the Reynolds number (or power of it).

```julia
julia> ρ = u"g/m^3"
g m⁻³

julia> μ = u"g/m/s"
g m⁻¹ s⁻¹

julia> u = u"m/s"
m s⁻¹

julia> ℓ = u"m"
m

julia> @setparameters ρ μ u ℓ
[ Info: Parameter(s) registered:
[ Info:  ρ = g m⁻³
[ Info:  μ = g m⁻¹ s⁻¹
[ Info:  u = m s⁻¹
[ Info:  ℓ = m

julia> pi_groups()
1-element Vector{String}:
 "ρ^(1//1)*μ^(-1//1)*u^(1//1)*ℓ^(1//1)"
```

Again, we can use [korsbo/Latexify.jl](https://github.com/korsbo/Latexify.jl) to display the adimensional group in Latex format:

![Reynolds number Pi group](img/reynoldsnumber_pi_group.png)

One can recognize this as the Reynolds number

![Reynolds number](img/reynoldsnumber.png)

## The internals

The [Buckingham-Pi Theorem](https://en.wikipedia.org/wiki/Buckingham_π_theorem) relies on the [Rank-nulity Theorem](https://en.wikipedia.org/wiki/Rank–nullity_theorem). A "parameter-to-dimension" matrix is composed, in which the columns correpond to the parameters and the rows to the collection of dimensions involved in the parameters. Each element in row i and column j corresponds to the power of the dimension i in the parameter j.

The number of adimensional groups is the dimension of the kernel of the matrix. And the adimensional groups are obtained from a basis of the null space.

When the powers are integers or rational numbers, which is usually the case, it is desirable to keep the type of these parameters when composing the matrix and when finding the null space and the associated adimensional Π groups.

While [ian-r-rose/buckinghampy](https://github.com/ian-r-rose/buckinghampy) uses [SymPy](https://www.sympy.org/en/index.html) for symbolic manipulation of the powers of the parameters, to retain these types, we simply rely on the strong type system of the Julia language.

The `LinearALgebra.nullspace`, however, uses the `LinearAlgebra.svd` factorization, which does not preserve the `Rational` type. The first version of our packaged used instead the ability of the `LU` decomposition in the [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) standard package to retain the `Rational` eltype of the matrices. However `LinearAlgebra.lu` factorization implements only partial pivoting and fails with singular matrices.

For this reason, we implemented our own LU factorization algorithm [UnitfulBuckinghamPi.lu_pq()](src/UnitfulBuckinghamPi.jl#L166), with full pivoting. Then, we find the null space from the U factor of the decomposition and by properly taking into account the column permutations used in the pivoting process.

## License

This package is licensed under the [MIT license](https://opensource.org/licenses/MIT) (see file [LICENSE](LICENSE) in the root directory of the project).
