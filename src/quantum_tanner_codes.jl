"""
Generate a pair of random classical codes (C_A, C_B) for quantum Tanner code construction [radebold2025explicit](@cite)

Returns a tuple (C_A, C_B) where each code is represented as a tuple (parity check matrix, generator matrix).
Both codes have block length Δ, with C_A having dimension ⌊ρΔ⌋ and C_B having dimension Δ - ⌊ρΔ⌋.

### Arguments
- `ρ::Real`: Target rate parameter, must satisfy 0 < ρ < 1/2
- `Δ::Integer`: Block length of the component codes
"""
function random_code_pair(ρ::Real, Δ::Int)
    H_A = uniformly_random_code_checkmatrix(ρ, Δ)
    G_A = dual_code(H_A)
    H_B = uniformly_random_code_checkmatrix(1-ρ, Δ)
    G_B = dual_code(H_B)
    H_A = Matrix{Int}(lift.(H_A))
    G_A = Matrix{Int}(lift.(G_A))
    H_B = Matrix{Int}(lift.(H_B))
    G_B = Matrix{Int}(lift.(G_B))
    return ((H_A, G_A), (H_B, G_B))
end

"""
The quantum Tanner code Q = (C₀, C₁) is defined by two classical Tanner codes
where Z-stabilizers: C₀ = T(Γ₀^□, (C_A ⊗ C_B)^⊥) andmX-stabilizers: C₁ = T(Γ₁^□, (C_A^⊥ ⊗ C_B^⊥)^⊥).

# Left-Right Cayley Complex

A Cayley graph `\\Gamma(V,E)`` provides a graph-theoretic representation
of a group G via a fixed generating set S that excludes the identity element.
The vertex set V corresponds to elements of G, with an edge connecting vertices
g and g' if and only if there exists ``s \\in S`` such that `g \\cdot s = g'``,
where ``\\cdot`` denotes the group operation. Edges are undirected if S is symmetric,
i.e., ``S = S^{-1}``.

A left-right Cayley complex extends this construction by incorporating both left
and right group actions. Specifically, we consider two symmetric generating sets A
and B and define a *bipartite* structure on the vertices.

Consider G be a finite group with symmetric generating sets ``A, B \\subseteq G``
such that ``\\langle A, B \\rangle = G`` and ``A = A^{-1}``, ``B = B^{-1}``. The
left-right Cayley complex ``\\Gamma(G,A,B)`` is defined as:

- Vertex set: ``V = V_0 \\cup V_1 = \\{g_i \\mid g \\in G, i \\in \\{0,1\\}\\}``
- Edge sets:
   - E_A = \\{(g_i, (ag)_j) \\mid a \\in A, g \\in G, i \\neq j\\}
   - E_B = \\{(g_i, (gb)_j) \\mid b \\in B, g \\in G, i \\neq j\\}

This construction yields a 2-dimensional complex whose faces are 4-cycles of the form:

```math
\\begin{aligned}
\\{g_i, (ag)_j, (gb)_j, (agb)_i \\mid i,j \\in \\{0,1\\}, i \\neq j\\}
\\end{aligned}
```

To ensure distinct opposite vertices in each face, we require that elements of A and B are not conjugates:

```math
\\begin{aligned}
\\forall a \\in A, b \\in B, g \\in G, \\quad ag \\neq gb
\\end{aligned}
```

Satisfying *Total Non-Conjugacy (TNC)* guarantees a proper 2D complex structure
where each vertex has degree ``\\Delta_A + \\Delta_B``, with ``\\Delta_A = |A|``
and ``\\Delta_B = |B|``. We typically take ``\\Delta_A = \\Delta_B = \\Delta`` [radebold2025explicit](@cite).

# Tensor codes

Classical linear block codes employ redundancy to encode information and
detect/correct errors. An `[n,k]`-code encodes k information bits into `n > k`
bits, described by either:
- A ``k \\times n`` generator matrix G whose rows span the code space
- An ``(n-k) \\times n`` parity check matrix H representing parity constraints

These satisfy ``GH^T = 0``.

For quantum Tanner codes, [radebold2025explicit](@cite) utilize a pair of binary
linear codes ``(C_A, C_B)`` where:
- ``C_A`` encodes ``\\rho\\Delta_A`` bits into ``\\Delta_A`` bits (``0 < \\rho < 1``)
- ``C_B`` encodes ``(1-\\rho)\\Delta_B`` bits into ``\\Delta_B`` bits

We construct tensor codes:

```math
\\begin{aligned}
C_0 = C_A \\otimes C_B, \\quad C_1 = C_A^\\perp \\otimes C_B^\\perp
\\end{aligned}
```

where ``\\dim(C_i \\otimes C_j) = \\dim(C_i)\\dim(C_j)`` and ``d(C_i \\otimes C_j) = d(C_i)d(C_j)``
for minimum distances.

# Quantum Tanner codes

To construct a quantum Tanner code, we begin with a left-right Cayley complex
``\\Gamma(G,A,B)`` built from a finite non-abelian group ``G``. Let ``Q`` denote
the complete set of faces of the complex, and for each vertex ``v \\in V``, let ``Q(v)``
be the set of faces incident to ``v`` [radebold2025explicit](@cite). We note that ``Q(v)``
is uniquely determined by pairs ``(a,b) \\in A \\times B`` for every vertex ``v``. The
physical qubits of the quantum code are placed on the faces of the complex, so the code
length is ``n = |Q|``.

We select two classical binary linear codes ``C_A`` and ``C_B``, where ``C_A``
encodes ``\\rho\\Delta_A`` logical bits into ``\\Delta_A`` bits for some
``0 < \\rho < 1``, and ``C_B`` encodes ``(1-\\rho)\\Delta_B`` logical bits into
``\\Delta_B`` bits. From these, we form the tensor codes ``C_0 = C_A \\otimes C_B``
and ``C_1 = C_A^\\perp \\otimes C_B^\\perp``. Since the number of columns of ``C_A`` is
``\\Delta_A``, we can label these columns with elements of ``A``. Having fixed this
association, we use the notation that codewords of ``C_A`` are binary vectors ``\\beta_A \\in \\mathbb{F}_2^A``.
Similarly, we associate the columns of ``C_B`` with elements of ``B``, so that
codewords of ``C_B`` are binary vectors ``\\beta_B \\subset \\mathbb{F}_2^B``. This
correspondence between the bits of the classical codes and group elements yields a natural
labeling of the columns of the tensor codes ``C_0`` and ``C_1`` by pairs ``(a,b) \\in A \\times B``
[radebold2025explicit](@cite).

To construct stabilizer generators on ``\\Gamma(G,A,B)`` using the classical codes
``C_0`` and ``C_1``, [radebold2025explicit](@cite)
define a function

```math
\\begin{aligned}
\\phi_v: A \\times B \\to Q(v)
\\end{aligned}
```

for each vertex v by ``\\phi_v(a,b) = {v, av, vb, avb}``, which maps a pair of
group generators to the unique face in ``Q(v)`` that it defines. One can verify that
``\\phi_v`` is bijective [radebold2025explicit](@cite). For each basis element
``\\beta`` of ``C_0``, we associate a set of pairs of group generators

```math
\\begin{aligned}
Z(\\beta) = \\{(a,b) \\mid \\beta_{(a,b)} = 1\\}
\\end{aligned}
```

corresponding to the nonzero entries of ``\\beta``. Each Z-stabilizer generator
of the quantum Tanner code is then specified by a choice of vertex ``v \\in V_0``
and classical codeword ``\\beta`` such that the ``Z``-stabilizer generator has
support equal to the set of faces ``\\phi_v(Z(\\beta))``. We can characterize this
stabilizer generator by a binary vector ``x \\in \\mathbb{F}2^Q``, where the |Q| qubits
are labelled by faces of the complex. For a given Z-stabilizer generator, the restriction
``x|{Q(v)}`` equals a basis element ``\\beta`` of ``C_0`` (based on a fixed ordering of
the faces), and is zero elsewhere. This yields ``\\dim(C_0)|V_0|`` Z-type stabilizer
generators, which correspond to codewords of C_0 placed locally at each vertex [radebold2025explicit](@cite).

[radebold2025explicit](@cite) repeat the same process for vertices ``v \\in V_1`` and
basis elements of ``C_1`` to produce ``\\dim(C_1)|V_1|`` X-type stabilizers at each vertex of that partition.

!!! note
    We note that there exist alternative formulations of quantum Tanner codes
    in the literature. The construction presented in [radebold2025explicit](@cite)
    utilizes the left-right Cayley complex structure where qubits are placed on the
    square faces (the 4-cycles of the form ``\\{g, ag, gb, agb\\}``) and stabilizers
    are defined via local tensor codes at vertices. In contrast, other approaches such
    as [gu2022efficient](@cite) employ a multigraph construction where qubits are
    identified with edges of the multigraphs ``\\mathcal{G}_0^\\square`` and ``\\mathcal{G}_1^\\square``
    derived from the left-right Cayley complex. These multigraphs have vertex set
    ``V_0 = G \\times \\{0\\}`` (respectively ``V_1 = G \\times \\{1\\}``). Edges correspond
    to squares ``q \\in Q`` connecting vertices ``(g,0)`` and ``(agb,0)`` via the relation
    ``g' = agb``. Stabilizers are built from Tanner codes associated with these multigraphs.

# Stabilizer Matrices

For each vertex v ∈ V₀ and basis element β ∈ C₀, we define the support set [radebold2025explicit](@cite): Z(β) = {(a,b) ∈ A×B | β_(a,b) = 1}

The corresponding Z-stabilizer generator has support φ_v(Z(β)), where φ_v: A×B → Q(v) is the bijective mapping
from generator pairs to incident faces [radebold2025explicit](@cite).

Similarly, for each vertex v ∈ V₁ and basis element β ∈ C₁, we define X-stabilizer generators
with support φ_v(Z(β)) [radebold2025explicit](@cite).

This yields dim(C₀) × |V₀| Z-type stabilizer generators and dim(C₁) × |V₁| X-type stabilizer generators

The resulting quantum code exhibits the Low-Density Parity-Check because each stabilizer generator
acts on at most Δ² qubits (where Δ = |A| = |B|) and each qubit is involved in at most 4ρ(1-ρ)Δ² stabilizer
generators. These bounds remain constant as |G| → ∞, ensuring the LDPC property [radebold2025explicit](@cite). 

# CSS Commutativity

All stabilizer generators of opposite type commute pairwise. The CSS orthogonality constraint 
C_X ⊂ C_Z^⊥ is fulfilled because when a C₀-generator (from V₀) and C₁-generator (from V₁)
have intersecting supports, their anchor vertices must be neighbors in the bipartite graph. If
connected by a B-edge, their local views share an A-set where C_A ⟂ C_A^⊥ ensures orthogonality. Note that
he A-edge case is analogous with C_B ⟂ C_B^⊥

# Quantum Tanner code parameters

For component codes C_A[Δ, ρΔ, δΔ] and C_B[Δ, (1-ρ)Δ, δΔ], the number of physical qubits is n = Δ²|G|/2, number of X-stabs is
dim(C₁) × |V₁| ≈ 2ρ(1-ρ)Δ²|G| and number of Z-stabs is dim(C₀) × |V₀| ≈ 2ρ(1-ρ)Δ²|G|. The resulting quantum code rate is
≥ (2ρ - 1)². For other properties, see [radebold2025explicit](@cite).


!!! note
    This is a newer version of the less well designed function [`gen_code`](@ref)(G, A, B, bipartite=true, use_same_local_code=false).
    It constructs the quantum Tanner code given a finite group G equipped with two *symmetric* generating sets A and B,
    alongside pairs of classical codes — comprising parity check and generator matrices — that are utilized in the
    construction of classical Tanner codes. To illustrate its application, the implementation can employ generating
    sets computed from the Morgenstern's explicit construction of Ramanujan graphs for odd prime power `q` generating sets.

Here is an example of new the `[[360, 8, 10]]` quantum Tanner code using Morgenstern generating sets

```jldoctest
julia> using QuantumExpanders; using Oscar; using QuantumClifford.ECC;

julia> l, i = 1, 2;

julia> q = 2^l;

julia> Δ = q+1;

julia> SL₂, B = morgenstern_generators(l, i);
[ Info: |SL₂(𝔽(4))| = 60

julia> A = alternative_morgenstern_generators(B, FirstOnly());

julia> H_A = [0 0 0 1;
              1 1 0 0];

julia> G_A = [1 1 0 0;
              0 0 1 0];

julia> H_B = [1 0 1;
              0 1 1];

julia> G_B = [1 1 1];

julia> classical_code_pair = ((H_A, G_A), (H_B, G_B));

julia> c = QuantumTannerCode(SL₂, A, B, classical_code_pair);

julia> import HiGHS; import JuMP;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS))
[ Info: Left-right Cayley complex Γ(G,A,B) square enumeration complete
[ Info: Group order |G| = 60, |A| = 4, |B| = 3
[ Info: Physical qubits: 360
[ Info: Left-right Cayley complex Γ(G,A,B): enumerated 360 faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)
[ Info: Left-right Cayley complex Γ(G,A,B) square enumeration complete
[ Info: Group order |G| = 60, |A| = 4, |B| = 3
[ Info: Physical qubits: 360
[ Info: Left-right Cayley complex Γ(G,A,B): enumerated 360 faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)
[ Info: Left-right Cayley complex Γ(G,A,B) square enumeration complete
[ Info: Group order |G| = 60, |A| = 4, |B| = 3
[ Info: Physical qubits: 360
[ Info: Left-right Cayley complex Γ(G,A,B): enumerated 360 faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)
[ Info: Left-right Cayley complex Γ(G,A,B) square enumeration complete
[ Info: Group order |G| = 60, |A| = 4, |B| = 3
[ Info: Physical qubits: 360
[ Info: Left-right Cayley complex Γ(G,A,B): enumerated 360 faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)
[ Info: Left-right Cayley complex Γ(G,A,B) square enumeration complete
[ Info: Group order |G| = 60, |A| = 4, |B| = 3
[ Info: Physical qubits: 360
[ Info: Left-right Cayley complex Γ(G,A,B): enumerated 360 faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)
[ Info: Left-right Cayley complex Γ(G,A,B) square enumeration complete
[ Info: Group order |G| = 60, |A| = 4, |B| = 3
[ Info: Physical qubits: 360
[ Info: Left-right Cayley complex Γ(G,A,B): enumerated 360 faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)
[ Info: Left-right Cayley complex Γ(G,A,B) square enumeration complete
[ Info: Group order |G| = 60, |A| = 4, |B| = 3
[ Info: Physical qubits: 360
[ Info: Left-right Cayley complex Γ(G,A,B): enumerated 360 faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)
[ Info: Left-right Cayley complex Γ(G,A,B) square enumeration complete
[ Info: Group order |G| = 60, |A| = 4, |B| = 3
[ Info: Physical qubits: 360
[ Info: Left-right Cayley complex Γ(G,A,B): enumerated 360 faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)
[ Info: Left-right Cayley complex Γ(G,A,B) square enumeration complete
[ Info: Group order |G| = 60, |A| = 4, |B| = 3
[ Info: Physical qubits: 360
[ Info: Left-right Cayley complex Γ(G,A,B): enumerated 360 faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)
[ Info: Squares incident to vertices: 720 at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)
(360, 8, 10)
```

### Fields
    $TYPEDFIELDS
"""
struct QuantumTannerCode <: AbstractCSSCode
    """The order of the underlying *finite* group"""
    group::Group
    """Symmetric generating set (closed under inverses) not containing the identity"""
    A::Vector{<:GroupElem}
    """Symmetric generating set (closed under inverses) not containing the identity"""
    B::Vector{<:GroupElem}
    """Tuple ((H_A, G_A), (H_B, G_B)) where (H_A, H_B) and (G_A, G_B) are parity-check and generator matrices, respectively."""
    classical_codes::Tuple{Tuple{Matrix{Int}, Matrix{Int}}, Tuple{Matrix{Int}, Matrix{Int}}}
    function QuantumTannerCode(group::Group,
                              A::Vector{<:GroupElem},
                              B::Vector{<:GroupElem}, 
                              classical_codes::Tuple{Tuple{Matrix{Int}, Matrix{Int}}, Tuple{Matrix{Int}, Matrix{Int}}})
        H_A, G_A = classical_codes[1]
        H_B, G_B = classical_codes[2]
        @assert size(H_A, 2) == length(A) "H_A parity check columns must match |A|"
        @assert size(G_A, 2) == length(A) "G_A generator columns must match |A|"
        @assert size(H_B, 2) == length(B) "H_B parity check columns must match |B|"
        @assert size(G_B, 2) == length(B) "G_B generator columns must match |B|"
        all(iszero, mod.(H_A*G_A', 2)) || @warn "C_A may not be a valid classical code: H_A*G_A^T ≠ 0"
        all(iszero, mod.(H_B*G_B', 2)) || @warn "C_B may not be a valid classical code: H_B*G_B^T ≠ 0"
        return new(group, A, B, classical_codes)
    end
end

"""
Enumerate all square incidences in the Left-Right Cayley Complex
following introduction by [dinur2022locally](@cite).

The Left-Right Cayley Complex X is an [incidence structure](https://en.wikipedia.org/wiki/Incidence_structure)
between:
- Vertices V = V₀ ∪ V₁ where V₀ = G×{0}, V₁ = G×{1}
- A-edges E_A = {(g,0), (ag,1)} for g ∈ G, a ∈ A ([double cover](https://en.wikipedia.org/wiki/Bipartite_double_cover) of left Cayley graph Cay(G,A))
- B-edges E_B = {(g,0), (gb,1)} for g ∈ G, b ∈ B ([double cover](https://en.wikipedia.org/wiki/Bipartite_double_cover) of right Cayley graph Cay(G,B))
- Squares Q = {(g,0), (ag,1), (gb,1), (agb,0)} for g ∈ G, a ∈ A, b ∈ B

Each square q ∈ Q corresponds to one physical qubit in the quantum Tanner code. Each square appears in two
natural local views [radebold2025explicit](@cite):
- From V₀ vertices: defines the graph Γ₀^□ = (V₀, Q) used for Z-stabilizers
- From V₁ vertices: defines the graph Γ₁^□ = (V₁, Q) used for X-stabilizers

We explicitly enumerates both incidences of each square to facilitate the Tanner code construction.

# Construction Framework

For each vertex v ∈ V, the set of incident faces Q(v) is uniquely determined by pairs (a,b) ∈ A×B.

The bijective mapping φ_v: A×B → Q(v) is defined as [radebold2025explicit](@cite): φ_v(a,b) = {v, av, vb, avb}

This establishes a natural labeling of qubits (*faces*) by generator pairs, allowing classical tensor codes
to be applied locally at each vertex [radebold2025explicit](@cite).

# Dihedral Ramanujan Graphs

The quantum Tanner code construction of [radebold2025explicit](@cite) utilizes a specific class
of *[Frobenius groups](https://en.wikipedia.org/wiki/Frobenius_group)*.

A finite group G is a Frobenius group if it can be expressed as a [semidirect product](https://en.wikipedia.org/wiki/Semidirect_product)

```math
\\begin{aligned}
G = N \\rtimes H
\\end{aligned}
```

where N (the *Frobenius kernel*) and H (the *Frobenius complement*) satisfy the condition that the
ratio ``r = \\frac{|N|-1}{|H|}`` is a positive integer. A canonical example, and the one primarily
utilized in ([radebold2025explicit](@cite), [hirano2016ramanujan](@cite)), is the *dihedral group*
of order 2p for an odd prime p:

```math
\\begin{aligned} 
D_{2p} = \\langle x, y \\mid x^p = y^2 = 1, yxy^{-1} = x^{-1} \\rangle.
\\end{aligned}
```

For ``D_{2p}``, the kernel is the cyclic subgroup ``N = \\langle x \\rangle \\cong \\mathbb{Z}_p`` and the
complement is ``H = \\langle y \\rangle \\cong \\mathbb{Z}_2``, yielding ``r = \\frac{p-1}{2}``.

[hirano2016ramanujan](@cite) provides a method to construct Ramanujan graphs from these groups. A Cayley
graph ``\\text{Cay}(G, S)`` is **Ramanujan** if its non-trivial eigenvalues ``\\lambda`` satisfy
``|\\lambda| \\leq 2\\sqrt{|S|-1}``. For a Frobenius group ``G = N \\rtimes H`` with ``r \\geq 4``,
the Theorem 3.3 of [hirano2016ramanujan](@cite) states that the maximum "covalency" ``\\hat{l}_{G, \\mathcal{S}_0}``
for which all corresponding normal Cayley graphs remain Ramanujan is given by the trivial bound:

```math
\\begin{aligned} 
\\hat{l}_{G, \\mathcal{S}_0} = l_0 = \\max \\{ l \\in \\mathcal{L} \\mid l \\le 2(\\sqrt{|G|}-1) \\}.
\\end{aligned}
```

For the dihedral group ``D_{2p}`` with ``p \\geq 11``, this bound specializes to Corollary 3.4 of [hirano2016ramanujan](@cite):

```math
\\begin{aligned} 
\\hat{l} = 2 \\left\\lfloor \\sqrt{2p} - \\frac{1}{2} \\right\\rfloor - 1.
\\end{aligned}
```

Thus, [radebold2025explicit](@cite) uses such dihedral groups  ``\D_{2p}`` and their symmetric generating
sets to construct the underlying expander graphs for the quantum Tanner code.

!!! note
    Through random search of clasical code pairs that are used for the construction of quantum Tanner codes,
    we found several new instances of these codes.

### Arguments
- `G`: A finite group
- `A`: Symmetric generating set (closed under inverses) not containing the identity
- `B`: Symmetric generating set (closed under inverses) not containing the identity
"""
function enumerate_squares(G::Group, A::Vector{<:GroupElem}, B::Vector{<:GroupElem})
    @assert is_symmetric_gen(A) "Definition 3.1: Set A must be symmetric generating set [dinur2022locally](@cite)"
    @assert is_symmetric_gen(B) "Definition 3.1: Set A must be symmetric generating set [dinur2022locally](@cite)"
    @assert !(one(G) in A) "Definition 3.1: Identity must not be in A [dinur2022locally](@cite)"
    @assert !(one(G) in B) "Definition 3.1: Identity must not be in A [dinur2022locally](@cite)"
    @assert is_nonconjugate(G, A, B) "Definition 3.6: ∀ a ∈ A, b ∈ B, g ∈ G, g⁻¹ag ≠ b [dinur2022locally](@cite)"
    idx_to_mat = collect(G)
    n = length(collect(G))
    mat_to_idx = Dict(mat=>i for (i,mat) in pairs(idx_to_mat))
    Q = []
    Q_redundant = ([], [])
    for g in idx_to_mat
        for a in A
            for b in B
                gᵢ = mat_to_idx[g]
                agᵢ = mat_to_idx[a*g] 
                gbᵢ = mat_to_idx[g*b]
                agbᵢ = mat_to_idx[a*g*b]
                aᵢ, bᵢ = findfirst(==(a), A), findfirst(==(b), B)
                aᵢₙᵥᵢ, bᵢₙᵥᵢ = findfirst(==(inv(a)), A), findfirst(==(inv(b)), B)
                # Typeₒ squares ({g₀, (a·g)₁, (g·b)₁, (a·g·b)₀}) represent Z-type stabilizers at V₀ vertices.
                typeₒ□ = [
                    0,
                    [g, a*g, g*b, a*g*b],
                    [gᵢ, agᵢ+n, gbᵢ+n, agbᵢ],
                    [a, b],
                    [aᵢ, bᵢ]
                ]
                typeₒ□_alt₁ = [
                    0, [a*g*b, g*b, a*g, g],
                    [agbᵢ, gbᵢ+n, agᵢ+n, gᵢ],
                    [inv(a), inv(b)], [aᵢₙᵥᵢ, bᵢₙᵥᵢ]
                ]
                typeₒ□_alt₂ = [
                    1, [g*b, a*g*b, g, a*g], 
                    [gbᵢ+n, agbᵢ, gᵢ, agᵢ+n],
                    [a, inv(b)], [aᵢ, bᵢₙᵥᵢ]
                ]
                typeₒ□_alt₃ = [
                    1, [a*g, g, a*g*b, g*b],
                    [agᵢ+n, gᵢ, agbᵢ, gbᵢ+n], 
                    [inv(a), b], [aᵢₙᵥᵢ, bᵢ]
                ]
                is_new□ = true
                for existing□ in Q
                    if (existing□[3] == typeₒ□[3] || existing□[3] == typeₒ□_alt₁[3] || existing□[3] == typeₒ□_alt₂[3] || existing□[3] == typeₒ□_alt₃[3])
                        is_new□ = false
                        break
                    end
                end
                if is_new□
                    push!(Q, typeₒ□)
                    □ᵢ = length(Q)
                    push!(typeₒ□, □ᵢ)
                    push!(typeₒ□_alt₁, □ᵢ)
                    push!(typeₒ□_alt₂, □ᵢ)
                    push!(typeₒ□_alt₃, □ᵢ)
                    push!(Q_redundant[1], typeₒ□)
                    push!(Q_redundant[1], typeₒ□_alt₁)
                    push!(Q_redundant[2], typeₒ□_alt₂)
                    push!(Q_redundant[2], typeₒ□_alt₃)
                end
                # Type 1 squares {g₁, (a·g)₀, (g·b)₀, (a·g·b)₁} represent X-stabilizers at V₁ vertices
                type□₁ = [
                    1,
                    [g, a*g, g*b, a*g*b],
                    [gᵢ + n, agᵢ, gbᵢ, agbᵢ + n],
                    [a, b],
                    [aᵢ, bᵢ]
                ]
                type□₁_alt₁ = [
                    1, [a*g*b, g*b, a*g, g],
                    [agbᵢ+n, gbᵢ, agᵢ, gᵢ+n],
                    [inv(a), inv(b)], [aᵢₙᵥᵢ, bᵢₙᵥᵢ]
                ]
                type□₁_alt₂ = [
                    0, [g*b, a*g*b, g, a*g],
                    [gbᵢ, agbᵢ+n, gᵢ+n, agᵢ],
                    [a, inv(b)], [aᵢ, bᵢₙᵥᵢ]
                ]
                type□₁_alt₃ = [
                    0, [a*g, g, a*g*b, g*b],
                    [agᵢ, gᵢ+n, agbᵢ+n, gbᵢ],
                    [inv(a), b], [aᵢₙᵥᵢ, bᵢ]
                ]
                is_new□ = true
                for existing□ in Q
                    if (existing□[3] == type□₁[3] || existing□[3] == type□₁_alt₁[3] || existing□[3] == type□₁_alt₂[3] || existing□[3] == type□₁_alt₃[3])
                        is_new□ = false
                        break
                    end
                end
                if is_new□
                    push!(Q, type□₁)
                    □ᵢ = length(Q)
                    push!(type□₁, □ᵢ)
                    push!(type□₁_alt₁, □ᵢ)
                    push!(type□₁_alt₂, □ᵢ)
                    push!(type□₁_alt₃, □ᵢ)
                    push!(Q_redundant[2], type□₁)
                    push!(Q_redundant[2], type□₁_alt₁)
                    push!(Q_redundant[1], type□₁_alt₂)
                    push!(Q_redundant[1], type□₁_alt₃)
                end
            end
        end
    end
    @info "Left-right Cayley complex Γ(G,A,B) square enumeration complete"
    @info "Group order |G| = $(length(collect(G))), |A| = $(length(A)), |B| = $(length(B))"
    @info "Physical qubits: $(length(Q))"
    @info "Left-right Cayley complex Γ(G,A,B): enumerated $(length(Q)) faces placed on 4-cycles {gᵢ, (a·g)ⱼ, (g·b)ⱼ, (a·g·b)ᵢ} where i,j ∈ {0,1}, i≠j [radebold2025explicit](@cite)"
    @info "Squares incident to vertices: $(length(Q_redundant[1])) at V₀ vertices (Z-type stabilizers) [radebold2025explicit](@cite)"
    @info "Squares incident to vertices: $(length(Q_redundant[2])) at V₁ vertices (X-type stabilizers) [radebold2025explicit](@cite)"
    return Q, Q_redundant
end

"""Convert redundant face list to matrix format for stabilizer matrix generation."""
function convert_squares_to_incidence_matrix(Q_redundant::Tuple)
    V₀□, V₁□ = Q_redundant
    total□ = length(V₀□)+length(V₁□)
    matrix□ = zeros(Int, total□, 8)
    rowᵢ = 1
    for (type□, list□) in enumerate([V₀□, V₁□])
        for square in list□
            matrix□[rowᵢ, 1] = square[1]
            matrix□[rowᵢ, 2] = square[3][1] 
            matrix□[rowᵢ, 3] = square[3][2]
            matrix□[rowᵢ, 4] = square[3][3]
            matrix□[rowᵢ, 5] = square[3][4]
            matrix□[rowᵢ, 6] = square[5][1]
            matrix□[rowᵢ, 7] = square[5][2]
            matrix□[rowᵢ, 8] = square[6]
            rowᵢ += 1
        end
    end
    return matrix□
end

"""
Construct X and Z stabilizer generators for the Quantum Tanner Code introduced in [leverrier2022quantum](@cite).

Returns the matrix for X-type stabilizer generators (dim(C₁) × |V₁| rows) and matrix for Z-type stabilizer generators (dim(C₀) × |V₀| rows).
"""
function parity_matrix_xz(c::QuantumTannerCode)
    Q, Q_red = enumerate_squares(c.group, c.A, c.B)
    squares_matrix = convert_squares_to_incidence_matrix(Q_red)
    group_order = order(c.group)
    classical_codes = c.classical_codes
    H_A, G_A = classical_codes[1]
    H_B, G_B = classical_codes[2]
    Δ_A = length(c.A)
    Δ_B = length(c.B)
    num_squares = size(squares_matrix, 1)
    # n = Δ_AΔ_B|G|/2
    n = Int(Δ_A*Δ_B*group_order/2)
    # Construct tensor codes as per [radebold2025explicit](@cite)
    # C₀ = C_A ⊗ C_B for Z-stabilizers [radebold2025explicit](@cite)
    # C₁ = C_A^⊥ ⊗ C_B^⊥ for X-stabilizers [radebold2025explicit](@cite)
    β₀ = kron(G_A, G_B) # Basis for C₀ [radebold2025explicit](@cite)
    β₁ = kron(H_A, H_B) # Basis for C₁ [radebold2025explicit](@cite)
    hz = Matrix{Int}(undef, 0, n)
    hx = Matrix{Int}(undef, 0, n)
    # Z-type stabilizers on V₀ vertices [radebold2025explicit](@cite)
    for v in 1:group_order  # V₀ vertices: 1 to |G|
        # Collect squares incident to vertex v
        squares_at_v = Matrix{Int}(undef, 0, 8)
        for square_idx in 1:num_squares 
            square = squares_matrix[square_idx, :]
            if (square[1] == 0 && square[2] == v) # Type 0 square at vertex v
                squares_at_v = [squares_at_v; square']
            end
        end
        # Generate Z-stabilizers from basis of C₀
        for basis_idx in 1:size(β₀, 1)
            β_row = β₀[basis_idx, :]
            β_matrix = transpose(reshape(β_row, (size(G_B, 2), size(G_A, 2))))
            stabs_vec = zeros(n)
            for incident_square in eachrow(squares_at_v)
                qubit_index = incident_square[8] # physical qubit index
                aᵢ = incident_square[6] # generator index from A
                bᵢ = incident_square[7] # generator index from B
                # set qubit according to classical code
                stabs_vec[qubit_index] = β_matrix[aᵢ, bᵢ]
            end
            hz = [hz; transpose(stabs_vec)]
            hz = unique(hz, dims=1)
        end
    end
    # X-Type stabilizers on V₁ vertices [radebold2025explicit](@cite)
    for v in (group_order+1):(2*group_order) # V₁ vertices: |G|+1 to 2|G|
        squares_at_v = Matrix{Int}(undef, 0, 8)
        for square_idx in 1:num_squares 
            square = squares_matrix[square_idx, :]
            if (square[1] == 1 && square[2] == v) # Type 1 square at vertex v
                squares_at_v = [squares_at_v; square']
            end
        end
        # Generate X-stabilizers from basis of C₁
        for basis_idx in 1:size(β₁, 1)
            β_row = β₁[basis_idx, :]
            β_matrix = transpose(reshape(β_row, (size(G_B, 2), size(G_A, 2))))
            stabs_vec = zeros(n)
            for incident_square ∈ eachrow(squares_at_v)
                qubit_index = incident_square[8]
                aᵢ = incident_square[6]
                bᵢ = incident_square[7]
                stabs_vec[qubit_index] = β_matrix[aᵢ, bᵢ]
            end
            hx = [hx; transpose(stabs_vec)]
            hx = unique(hx, dims=1)
        end
    end
    hx, hz = unique(hx, dims=1), unique(hz, dims=1)
    return Int.(hx), Int.(hz)
end

parity_matrix_x(c::QuantumTannerCode) = parity_matrix_xz(c)[1]

parity_matrix_z(c::QuantumTannerCode) = parity_matrix_xz(c)[2]

code_n(c::QuantumTannerCode) = Int(order(c.group))*length(c.A)*length(c.B)÷2

code_k(c::QuantumTannerCode) = code_n(c) - rank(matrix(GF(2), parity_matrix_x(c))) - rank(matrix(GF(2), parity_matrix_z(c)))
