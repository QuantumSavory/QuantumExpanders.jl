# [Reproducing Leverrier--Rozendaal--Zémor instances](@id leverrier-instances)

Leverrier, Rozendaal, and Zémor reported several small quantum Tanner codes in
*Small quantum Tanner codes from left--right Cayley complexes*
[leverrier2025small](@cite). The accompanying ancillary arXiv files provide the
resulting ``H_X`` and ``H_Z`` matrices, but not the ordered multisets ``A`` and
``B`` or the local-code column permutations used to construct them.

We provide explicit construction data for three of the reported weight-nine parameter sets. It also demonstrates that [`QuantumTannerViaLeftRightActions`](@ref) produces codes with the same length, dimension, check weight, and observed distance bounds.

!!! note "What is reproduced here?"
    The examples below reproduce the **published parameters** from explicit
    lift data. We estimated the distances using [sQetch](https://github.com/a7b/yarn): 50 million trials for each length-144 code and 100 million trials for the length-288 code.

## Recovered parameter sets

All three examples use the shortened Hamming ``[6,3,3]`` code on both sides,
so every stabilizer generator has weight at most nine.

| Published parameters | Group | Recovered lower bounds | [sQetch](https://github.com/a7b/yarn) trials |
|---|---|---|---:|
| ``[[144,8,12]]`` | ``C_2 \times C_2`` | ``(d_X,d_Z)=(\leq 12, \leq12)`` | 50,000,000 |
| ``[[144,12,11]]`` | ``C_2 \times C_2`` | ``(d_X,d_Z)=(\leq 11, \leq 11)`` | 50,000,000 |
| ``[[288,8,19]]`` | ``C_8`` | ``(d_X,d_Z)=(\leq 19, \leq 19)`` | 100,000,000 |

The distances in the paper and in this table come from randomized searches.

## Shared local code

The paper uses the following parity-check and generator matrices for the
shortened Hamming code:

```jldoctest lifted-example
julia> using QuantumExpanders, Oscar, QECCore;

julia> H633 = [1 0 0 0 1 1;
               0 1 0 1 0 1;
               0 0 1 1 1 0];

julia> G633 = [0 1 1 1 0 0;
               1 0 1 0 1 0;
               1 1 0 0 0 1];

julia> iszero(mod.(H633 * G633', 2))
true

julia> p633 = [1, 2, 6, 4, 5, 3];
```

The permutation `p633` specifies the second local code on each side. It is
passed as both `p1` and `p2` below.

## The ``[[144,8,12]]`` instance

Oscar identifies ``C_2\times C_2`` as `small_group(4, 2)`. In the permutation
representation used here, write ``x=(1\;2)`` and ``y=(3\;4)``.

```jldoctest lifted-example
julia> V4 = codomain(isomorphism(PermGroup, small_group(4, 2)));

julia> x = cperm(V4, [1,2]);

julia> y = cperm(V4, [3,4]);

julia> e = one(V4);

julia> A = [e, e, e, y, x, x*y];

julia> B = [e, e, y, y, x, x*y]; 

julia> code = QuantumTannerViaLeftRightActions(V4, A, B, H633, G633, H633, G633;p1=p633, p2=p633,);

julia> (code_n(code), code_k(code))
(144, 8)
```

## The ``[[144,12,11]]`` instance

This code uses the same group, local code, permutation, and ``A`` multiset.
Only the ordered multiset ``B`` changes.

```jldoctest lifted-example
julia> A = [e, e, e, y, x, x*y];

julia> B = [e, e, y, y, x, x];

julia> code = QuantumTannerViaLeftRightActions(V4, A, B, H633, G633, H633, G633, p1=p633, p2=p633,);

julia> (code_n(code), code_k(code))
(144, 12)
```

## The ``[[288,8,19]]`` instance

Oscar identifies ``C_8`` as `small_group(8, 1)`. The explicit permutation
representation below records the recovered multisets without depending on a
choice of abstract cyclic generator.

```jldoctest lifted-example
julia> C8 = codomain(isomorphism(PermGroup, small_group(8, 1)));

julia> a = cperm(C8, [1,7,5,3], [2,8,6,4]);

julia> ainv = cperm(C8, [1,3,5,7], [2,4,6,8]);

julia> b = cperm(C8, [1,8,7,6,5,4,3,2]);

julia> c = cperm(C8, [1,4,7,2,5,8,3,6]);

julia> d = cperm(C8, [1,6,3,8,5,2,7,4]);

julia> A = [a, a, ainv, b, c, d];

julia> B = [one(C8), a, b, b, c, d];

julia> code = QuantumTannerViaLeftRightActions(C8, A, B, H633, G633, H633, G633;p1=p633, p2=p633,);

julia> (code_n(code), code_k(code))
(288, 8)
```

## References

```@bibliography
Pages = ["leverrier_instances.md"]
Canonical = false
```
