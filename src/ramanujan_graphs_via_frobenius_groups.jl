"""
Generate a normal Cayley subset of a finite group G. A Cayley subset S is called
normal if it is a union of [conjugacy classes](https://en.wikipedia.org/wiki/Conjugacy_class)
of G. [hirano2016ramanujan](@cite) studies the spectral properties of Cayley graphs built from such subsets.

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

The properties of normal Cayley subsets find direct application in the construction
of Ramanujan graphs for specific group families, as demonstrated in the following results:

[mehry2023ramanujan](@cite) investigates the construction of Ramanujan Cayley graphs for
*sporadic* groups and *linear* groups, with a particular focus on the *Mathieu* groups M(9), M(10),
and M(11), as well as the *Suzuki* groups ``\\text{Sz}(q)``. Additionally, they examines the special
linear group SL(2, q). By leveraging the character tables of these groups and identifying
appropriate normal symmetric generating subsets, the [mehry2023ramanujan](@cite) derive conditions
under which the corresponding Cayley graphs satisfy the Ramanujan property.

[droll2010classification](@cite) classifies all unitary Cayley graphs ``X_n`` that are Ramanujan.
The unitary Cayley graph ``X_n`` is constructed on the additive group ``\\mathbb{Z}/n\\mathbb{Z}``,
where vertices represent integers modulo n, and edges connect two vertices if their difference is
a multiplicative unit modulo n (i.e., ``\\gcd(a-b, n) = 1``). [droll2010classification](@cite) provides
a complete characterization: ``X_n`` is Ramanujan—meaning it satisfies the optimal spectral gap
condition ``\\lambda(X_n) \\leq 2\\sqrt{k-1}`` for k-regular graphs—if and only if n belongs to
one of six explicit families based on its prime factorization, such as powers of 2, primes, or
products of two primes under specific constraints. 

[wang1998normal](cite) prove that every finite group admits a normal Cayley grap*, except for
``\\mathbb{Z}_4 \\times \\mathbb{Z}_2`` and ``Q_8 \\times \\mathbb{Z}_2^r``, and that every finite
group has a normal Cayley digraph. 

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
