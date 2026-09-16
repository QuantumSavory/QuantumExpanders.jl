# [Getting started](@id getting-started)

This page takes a code from constructor inputs to the standard CSS interface.
The example is deliberately small so that every matrix can be inspected
directly.

## Installation

Install the package from GitHub in Julia's package mode:

```julia
pkg> add https://github.com/QuantumSavory/QuantumExpanders.jl.git
```

Then load the package together with Oscar:

```julia
using QuantumExpanders, Oscar
```

Julia 1.12 or later is required.

## Construct a lifted code

Let ``G=C_2`` and write its nonidentity element as ``g``. We use the two-element
multisets ``A=B=[1_G,g]`` and the length-two repetition code with parity-check
matrix ``H=(1\;1)``:

```julia
G = cyclic_group(2)
g = gens(G)[1]

A = [one(G), g]
B = [one(G), g]
H = [1 1]

code = QuantumTannerViaLeftRightActions(G, A, B, H, H)
```

The parity-only constructor computes the dual generator matrices internally.
The blocklength is

```math
n=|G|\,|A|\,|B|=2\cdot2\cdot2=8.
```

```julia
code_n(code), code_k(code)
# (8, 2)
```

## Inspect and validate the CSS code

`QuantumTannerViaLeftRightActions` is an `AbstractCSSCode`, so it implements the
usual QECCore interface:

```julia
hx = parity_matrix_x(code)
hz = parity_matrix_z(code)

size(hx), size(hz)
iszero(mod.(hx * hz', 2))
```

The final expression checks CSS orthogonality,
``H_XH_Z^\mathsf{T}=0``. The convenience accessor returns both matrices:

```julia
hx, hz = parity_matrix_xz(code)
```

The code can now be passed directly to functions from QECCore and
QuantumClifford.ECC.

## Accepted matrix types

The lifted constructors accept ordinary Julia integer matrices and Oscar
matrices over ``\mathbb{F}_2``. For example, this is equivalent to the previous
choice of `H`:

```julia
H = matrix(GF(2), [1 1])
code = QuantumTannerViaLeftRightActions(G, A, B, H, H)
```

The parity-only constructor also accepts the matrix returned by [`dual_code`](@ref):

```julia
H = dual_code([1 1])
code = QuantumTannerViaLeftRightActions(G, A, B, H, H)
```

Inputs are reduced modulo two, and every supplied parity-check/generator pair
is checked for duality.

## Choose the next page

- To understand the square-complex geometry, continue to
  [Quantum Tanner Codes](@ref quantum-tanner-codes).
- To understand the lifted matrices and constructor forms, continue to
  [Quantum Tanner Codes via Left-Right Actions](@ref quantum-tanner-left-right-actions).
- To rebuild the published examples, continue to
  [Reproducing the manuscript instances](@ref paper-instances).

!!! note "Distance is a separate computation"
    Constructing a code determines ``H_X``, ``H_Z``, ``n``, and ``k``. It does
    not automatically determine the minimum distance. Randomized estimators
    return upper bounds that may improve when more trials are run; exact
    algorithms can become expensive at moderate blocklength.
