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
    H_B = uniformly_random_code_checkmatrix(ρ, Δ)
    G_B = dual_code(H_B)
    H_A = Matrix{Int}(lift.(H_A))
    G_A = Matrix{Int}(lift.(G_A))
    H_B = Matrix{Int}(lift.(H_B))
    G_B = Matrix{Int}(lift.(G_B))
    return ((H_A, G_A), (H_B, G_B))
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

Returns incidence_matrix`: Matrix where each row represents a square incidence with:
[`bipartition_type`, `anchor_vertex`, `square_vertices`, `vertex_indices`, `generator_indices`,
`a_idx`, `b_idx`, `square_idx`] where each field means:

-  `bipartition_type::Int`: 
    - 0: Square viewed from a V₀ vertex (for Z-stabilizer assignment)
    - 1: Square viewed from a V₁ vertex (for X-stabilizer assignment)
    - Each physical square appears twice: once as type 0 and once as type 1
  
- `anchor_vertex::Int`: 
    - The reference vertex from which this square is viewed
    - For type 0: vertex index in V₀ (range: 1 to |G|)
    - For type 1: vertex index in V₁ (range: |G|+1 to 2|G|)
    - This is the vertex that "sees" this square in its local view
  
- `square_vertices::Vector{Int}`: 
    - All 4 vertices of the square in cyclic order
    - Format: [v₀₁, v₁₁, v₁₂, v₀₂] where v₀ᵢ ∈ V₀, v₁ᵢ ∈ V₁
    - For type 0: [g, ag, gb, agb] with proper bipartite indexing
    - For type 1: [ag, gb, agb, g] (rotated view)
  
- `vertex_indices::Vector{Int}`: 
    - Integer indices of the square vertices (same as square_vertices)
    - V₀ vertices: 1 to |G|, V₁ vertices: |G|+1 to 2|G|
    - Used for graph construction and neighborhood queries
  
- `generator_indices::Vector{Int}`: 
    - Pair [a_idx, b_idx] identifying which generators created this square
    - a_idx ∈ {1,...,|A|} indexes into the generating set A
    - b_idx ∈ {1,...,|B|} indexes into the generating set B
    - Together they uniquely identify the square up to the anchor vertex
  
- `a_idx::Int`: 
    - Index of the left generator a ∈ A that created this square
    - Used as the "row coordinate" in the local A×B grid at each vertex
  
- `b_idx::Int`: 
    - Index of the right generator b ∈ B that created this square  
    - Used as the "column coordinate" in the local A×B grid at each vertex
  
- `square_idx::Int`: 
    - Identifier for this square incidence
    - Note: Each physical square has two incidence records (type 0 and type 1)
    - Used as the physical qubit index in the quantum code
    - Range: 1 to 2|G||A||B| (total incidences), but physical qubits = |G||A||B|

### Arguments
- `G`: A finite group
- `A`: Symmetric generating set (closed under inverses) not containing the identity
- `B`: Symmetric generating set (closed under inverses) not containing the identity
"""
function enumerate_square_incidences(G, A, B)
    @assert is_symmetric_gen(A) "Definition 3.1: Set A must be symmetric generating set [dinur2022locally](@cite)"
    @assert is_symmetric_gen(B) "Definition 3.1: Set A must be symmetric generating set [dinur2022locally](@cite)"
    @assert !(one(G) in A) "Definition 3.1: Identity must not be in A [dinur2022locally](@cite)"
    @assert !(one(G) in B) "Definition 3.1: Identity must not be in A [dinur2022locally](@cite)"
    @assert is_nonconjugate(G, A, B) "Definition 3.6: ∀ a ∈ A, b ∈ B, g ∈ G, g⁻¹ag ≠ b [dinur2022locally](@cite)"
    group_elements = collect(G)
    element_to_idx = Dict(elem => i for (i, elem) in enumerate(group_elements))
    group_size = length(group_elements)
    square_incidences = []
    square_counter = 0
    @showprogress for g in group_elements
        g_idx = element_to_idx[g]
        for (a_idx, a) in enumerate(A)
            ag = a*g
            ag_idx = element_to_idx[ag]
            for (b_idx, b) in enumerate(B)
                square_counter += 1
                # Lets find all vertices of the square q = {(g,0), (ag,1), (gb,1), (agb,0)}
                gb = g*b
                agb = a*g*b
                gb_idx = element_to_idx[gb]
                agb_idx = element_to_idx[agb]
                # TYPE 0 INCIDENCE: Square as seen from V₀ vertex g ∈ V₀
                # This defines an edge in Γ₀^□ = (V₀, Q) connecting g and agb
                # Used for constructing Z-stabilizers via Tanner code T(Γ₀^□, C_A ⊗ C_B) see II. Section 3 of [radebold2025explicit](@cite).
                type0_incidence = Any[
                    0, # bipartition_type: V₀
                    g_idx, # anchor_vertex: v ∈ V₀ seeing this square
                    [g_idx, ag_idx + group_size, gb_idx + group_size, agb_idx],# square vertices
                    [g_idx, ag_idx + group_size, gb_idx + group_size, agb_idx], # vertex indices  
                    [a_idx, b_idx], # generator_indices identifying the square
                    a_idx,            # a_idx: coordinate in A for local view
                    b_idx,  # b_idx: coordinate in B for local view
                    square_counter  # square_idx: unique identifier for qubit
                ]
                push!(square_incidences, type0_incidence)
                # TYPE 1 INCIDENCE: Square as seen from V₁ vertex ag ∈ V₁  
                # This defines an edge in Γ₁^□ = (V₁, Q) connecting ag and gb
                # Used for constructing X-stabilizers via Tanner code T(Γ₁^□, C_A^⊥ ⊗ C_B^⊥)
                square_counter += 1
                type1_incidence = Any[
                    1, # bipartition_type: V₁
                    ag_idx + group_size, # anchor_vertex: v ∈ V₁ seeing this square
                    [ag_idx + group_size, gb_idx + group_size, agb_idx, g_idx], # rotated view
                    [ag_idx + group_size, gb_idx + group_size, agb_idx, g_idx], # vertex indices
                    [a_idx, b_idx], # same generator_indices (same physical square)
                    a_idx, # same a_idx
                    b_idx, # same b_idx  
                    square_counter  # new incidence index
                ]
                push!(square_incidences, type1_incidence)
            end
        end
    end
    num_incidences = length(square_incidences)
    incidence_matrix = Matrix{Any}(undef, num_incidences, 8)
    for (i, incidence_row) in enumerate(square_incidences)
        for j in 1:8
            incidence_matrix[i, j] = incidence_row[j]
        end
    end
    @info "Square complex construction complete" 
    @info "Group order |G| = $group_size, |A| = $(length(A)), |B| = $(length(B))"
    @info "Total square incidences: $num_incidences"
    @info "Physical squares (qubits): $(group_size * length(A) * length(B))"
    @info "Graph Γ₀^□: $(group_size) vertices, $(group_size * length(A) * length(B)) edges"
    @info "Graph Γ₁^□: $(group_size) vertices, $(group_size * length(A) * length(B)) edges"
    return incidence_matrix
end

"""Construct X and Z stabilizer generators for the Quantum Tanner Code introduced in [leverrier2022quantum](@cite).

Returns the matrix for X-type stabilizer generators (dim(C₁) × |V₁| rows) and matrix for Z-type stabilizer generators (dim(C₀) × |V₀| rows).

The quantum code Q = (C₀, C₁) is defined by two classical Tanner codes where Z-stabilizers: C₀ = T(Γ₀^□, (C_A ⊗ C_B)^⊥) and
X-stabilizers: C₁ = T(Γ₁^□, (C_A^⊥ ⊗ C_B^⊥)^⊥).

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

### Arguments
- `group_size`: The order of the underlying *finite* group
- `square_incidences`: Matrix from `enumerate_square_incidences` containing the φ_v mappings
- `classical_code_pair`: Tuple ((H_A, G_A), (H_B, G_B)) where:
  - H_A, H_B: parity-check matrices
  - G_A, G_B: generator matrices
  - C_A = ker(H_A), C_B = ker(H_B) are the classical component codes.
"""
function parity_matrix(group_size, square_incidences, classical_code_pair)
    # Extract classical code components
    parity_check_A, generator_A = classical_code_pair[1]
    parity_check_B, generator_B = classical_code_pair[2]
    num_incidences = size(square_incidences, 1)
    Δ_A = size(generator_A, 2)
    Δ_B = size(generator_B, 2)
    # Number of physical qubits = number of squares in complex
    num_qubits = Int(Δ_A*Δ_B*group_size/2)
    # Z-stabilizers use constraint code: (C_A ⊗ C_B)^⊥ = C_A^⊥ ⊗ F₂^B + F₂^A ⊗ C_B^⊥
    # But we generate them via the dual: we take generators of C_A ⊗ C_B
    z_constraint_generators = kron(generator_A, generator_B)
    # X-stabilizers use constraint code: (C_A^⊥ ⊗ C_B^⊥)^⊥ = C_A ⊗ F₂^B + F₂^A ⊗ C_B  
    # Generated via generators of C_A^⊥ ⊗ C_B^⊥
    x_constraint_generators = kron(parity_check_A, parity_check_B)
    # Z-stabilizers: From Tanner code T(Γ₀^□, (C_A ⊗ C_B)^⊥)
    # For each v ∈ V₀, constrain local view to be orthogonal to C_A ⊗ C_B
    z_stabs = Matrix{Int}(undef, 0, num_qubits)
    for v0_vertex in 1:group_size
        local_squares = [] # Squares incident to this V₀ vertex
        # Collect all Type 0 square incidences anchored at this V₀ vertex: φ_v: A×B → Q(v) for v ∈ V₀
        for incidence_idx in 1:num_incidences
            incidence_data = square_incidences[incidence_idx, :]
            if length(incidence_data) >= 2 && incidence_data[1] == 0 && incidence_data[2] == v0_vertex
                push!(local_squares, incidence_data)
            end
        end
        # Apply each generator of C_A ⊗ C_B to the local square arrangement where each generator corresponds to a basis element β ∈ C₀
        for constraint_idx in 1:size(z_constraint_generators, 1)
            constraint_vector = z_constraint_generators[constraint_idx, :]
            constraint_matrix = transpose(reshape(constraint_vector, (Δ_A, Δ_B)))
            stabilizer = zeros(num_qubits)
            for square_incidence in local_squares
                if length(square_incidence) >= 8
                    qubit_idx = square_incidence[8]
                    a_coordinate = square_incidence[6]
                    b_coordinate = square_incidence[7]
                    # Assign coefficient from constraint code based on local coordinates which implements the support set Z(β) = {(a,b) | β_(a,b) = 1} [radebold2025explicit](@cite).
                    if 1 <= qubit_idx <= num_qubits && 1 <= a_coordinate <= Δ_A && 1 <= b_coordinate <= Δ_B
                        stabilizer[qubit_idx] = constraint_matrix[a_coordinate, b_coordinate]
                    end
                end
            end
            z_stabs = [z_stabs; transpose(stabilizer)]
            z_stabs = unique(z_stabs, dims=1)
        end
    end
    # X-stabilizers: From Tanner code T(Γ₁^□, (C_A^⊥ ⊗ C_B^⊥)^⊥)  
    # For each v ∈ V₁, constrain local view to be orthogonal to C_A^⊥ ⊗ C_B^⊥
    x_stabs = Matrix{Int}(undef, 0, num_qubits)
    for v1_vertex in (group_size + 1):(2 * group_size)
        local_squares = []  # Squares incident to this V₁ vertex
        # Collect all Type 1 square incidences anchored at this V₁ vertex which implements the mapping φ_v: A×B → Q(v) for v ∈ V₁ [radebold2025explicit](@cite).
        for incidence_idx in 1:num_incidences
            incidence_data = square_incidences[incidence_idx, :]
            if length(incidence_data) >= 2 && incidence_data[1] == 1 && incidence_data[2] == v1_vertex
                push!(local_squares, incidence_data)
            end
        end
        # Apply each generator of C_A^⊥ ⊗ C_B^⊥ to the local square arrangement where generator corresponds to a basis element β ∈ C₁
        for constraint_idx in 1:size(x_constraint_generators, 1)
            constraint_vector = x_constraint_generators[constraint_idx, :]
            constraint_matrix = transpose(reshape(constraint_vector, (size(parity_check_A, 2), size(parity_check_B, 2))))
            stabilizer = zeros(num_qubits)
            for square_incidence in local_squares
                if length(square_incidence) >= 8
                    qubit_idx = square_incidence[8]
                    a_coordinate = square_incidence[6]  
                    b_coordinate = square_incidence[7]   
                    if 1 <= qubit_idx <= num_qubits && 1 <= a_coordinate <= size(constraint_matrix, 1) && 1 <= b_coordinate <= size(constraint_matrix, 2)
                        stabilizer[qubit_idx] = constraint_matrix[a_coordinate, b_coordinate]
                    end
                end
            end
            x_stabs = [x_stabs; transpose(stabilizer)]
            x_stabs = unique(x_stabs, dims=1)
        end
    end
    hx, hz = unique(x_stabs, dims=1), unique(z_stabs, dims=1)
    hx, hz = [Int.(hx), Int.(hz)]
    # expected quantum rate and LDPC parameters [radebold2025explicit](@cite).
    ρ_A = size(generator_A, 1) / Δ_A
    expected_rate = round((2ρ_A - 1)^2, digits=4)
    max_stabilizer_weight = Δ_A*Δ_B
    max_qubit_degree = 4*ρ_A*(1-ρ_A)*Δ_A*Δ_B
    @info "Physical qubits: $num_qubits"
    @info "Expected quantum rate: ≥ $expected_rate"
    @info "LDPC properties: max stabilizer weight = $max_stabilizer_weight, max qubit degree = $max_qubit_degree"
    iszero(mod.(hx*hz',2))
    return hx, hz
end
