function coefficient_matrix_and_monomials(F::Vector{<:MPolyRingElem})
    Kx = parent(first(F))
    distinct_monomials = unique(vcat(collect.(monomials.(F))...))
    coefficient_matrix = matrix([[coeff(f,mon) for mon in distinct_monomials] for f in F])
    return coefficient_matrix, distinct_monomials
end

function specialize(F::Vector{<:MPolyRingElem}, choice_of_parameters::Vector{<:Union{Int, RingElem}})
    Kax = parent(first(F))
    Ka = coefficient_ring(Kax)
    K = base_ring(Ka)
    Kx, x = polynomial_ring(K, symbols(Kax))
    phi = hom(Kax, Kx, c -> evaluate(c, choice_of_parameters), x)
    return phi.(F)
end

