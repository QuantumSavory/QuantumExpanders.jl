"""
Generate a normal Cayley subset of a finite group G. A Cayley subset S is called
normal if it is a union of [conjugacy classes](https://en.wikipedia.org/wiki/Conjugacy_class)
of G. [hirano2016ramanujan](@cite) which studies the spectral properties of Cayley graphs
built from such subsets.

We take all conjugacy classes of G and excluding the trivial conjugacy class containing
only the identity element. The resulting set S automatically satisfies several important
properties: it is symmetric (``S = S^{-1}``) because conjugacy classes are closed under
inverses, and it is normal (``g^{-1}Sg = S`` for all ``g \\in G``) by construction.

In the context of [hirano2016ramanujan](@cite), normal Cayley subsets are particularly
valuable because they allow for explicit computation of graph eigenvalues using [character
theory](https://en.wikipedia.org/wiki/Character_theory).For any normal Cayley subset S, the
eigenvalues of the Cayley graph X(S) can be expressed as:

```math
\\begin{aligned}
\\lambda_{\\chi} = \\frac{1}{\\chi(1)} \\sum_{s \\in S} \\chi(s)
\\end{aligned}
```

where ``\\chi`` ranges over the irreducible characters of G.

The work by [hirano2016ramanujan](@cite) establishes that for Frobenius groups ``G = N \\rtimes H``
with sufficiently large kernel-complement ratio (``r = (|N|-1)/|H| \\geq 4``), there exists a well-defined
bound ``l_0`` such that all normal Cayley subsets with covalency ``l \\leq l_0`` produce Ramanujan graphs.
For dihedral groups ``D_{2p}`` in particular, this bound takes the explicit form

```math
\\begin{aligned}
l_0 = 2\\lfloor \\sqrt{2p} - \\frac{1}{2} \\rfloor - 1
\\end{aligned}
```

when ``p \\geq 11``.

# Examples

```jldoctest
julia> using Oscar

julia> G = dihedral_group(10);

julia> S = normal_cayley_subset(G)
9-element Vector{PcGroupElem}:
 f1
 f1*f2^2
 f1*f2^4
 f1*f2
 f1*f2^3
 f2
 f2^4
 f2^2
 f2^3

julia> length(S)
9

julia> all(s^-1 in S for s in S)
true

julia> all(g^-1 * s * g in S for s in S, g in G)
true
```

### Arguments
- `G::Group`: A finite group
"""
function normal_cayley_subset(G::Group)
    classes = conjugacy_classes(G)
    identity_class = nothing
    for cls in classes
        if length(cls) == 1 && first(cls) == one(G)
            identity_class = cls
            break
        end
    end
    non_trivial_classes = filter(cls -> cls != identity_class, classes)
    S = union(non_trivial_classes...)
    return collect(S)
end
