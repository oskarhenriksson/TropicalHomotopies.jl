
# Define ring of parameters
A, (k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k31, c1, c2, c3, c4, c5) =
    rational_function_field(QQ, ["k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11", "k12", "k13", "k14", "k15", "k16", "k17", "k18", "k19", "k20", "k21", "k22", "k23", "k24", "k25", "k26", "k27", "k28", "k29", "k30", "k31", "c1", "c2", "c3", "c4", "c5"]);

# Define ring of parameterized polynomials
B, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19) = polynomial_ring(A, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19"]);

# Define the square(!) system of steady state equations and conservation laws
steadyStateEqs = [
    -k1 * x1 + k2 * x2, # x1
    k1 * x1 - (k2 + k26) * x2 + k27 * x3 - k3 * x2 * x4 + (k4 + k5) * x14, # x2
    k26 * x2 - k27 * x3 - k14 * x3 * x6 + (k15 + k16) * x15, # x3
    -k3 * x2 * x4 - k9 * x4 * x10 + k4 * x14 + k8 * x16 + (k10 + k11) * x18, # x4
    -k28 * x5 + k29 * x7 - k6 * x5 * x8 + k5 * x14 + k7 * x16, # x5
    -k14 * x3 * x6 - k20 * x6 * x11 + k15 * x15 + k19 * x17 + (k21 + k22) * x19, # x6
    k28 * x5 - k29 * x7 - k17 * x7 * x9 + k16 * x15 + k18 * x17, # x7
    -k6 * x5 * x8 + (k7 + k8) * x16, # x8
    -k17 * x7 * x9 + (k18 + k19) * x17, # x9
    k12 - (k13 + k30) * x10 - k9 * x4 * x10 + k31 * x11 + k10 * x18, # x10
    -k23 * x11 + k30 * x10 - k31 * x11 - k20 * x6 * x11 - k24 * x11 * x12 + k25 * x13 + k21 * x19, # x11
    -k24 * x11 * x12 + k25 * x13, # x12
    -k24 * x11 * x12 + k25 * x13, # x13
    k3 * x2 * x4 - (k4 + k5) * x14, # x14
    k14 * x3 * x6 - (k15 + k16) * x15, # x15
    -k6 * x5 * x8 + (k7 + k8) * x16, # x16
    -k17 * x7 * x9 + (k18 + k19) * x17, # x17
    k9 * x4 * x10 - (k10 + k11) * x18, # x18
    k20 * x6 * x11 - (k21 + k22) * x19]; # x19

conservationLaws = [
    (x1 + x2 + x3 + x14 + x15) - c1,
    (x4 + x5 + x6 + x7 + x14 + x15 + x16 + x17 + x18 + x19) - c2,
    (x8 + x16) - c3,
    (x9 + x17) - c4,
    (x12 + x13) - c5];

polynomialSystemParametrised = vcat([steadyStateEq for (i, steadyStateEq) in enumerate(steadyStateEqs) if i ∉ [3, 4, 8, 9, 12]],
    conservationLaws)


# Simplicifaction:
# Identify all binomial equations and use them to eliminate variables
# Note that it is important that the equations are binomial,
# as replacing x by y does not change the verticalilty of the parametrization,
# but replacing x by y+z does.
###

function construct_substitution_homomorphism(polynomialRing::MPolyRing,
    variablesToBeSubstituted::Vector{<:MPolyRingElem},
    polynomialsToSubstitute::Vector{<:MPolyRingElem})

    ###
    # Construct codomain of the substitution homomorphism
    ###
    variableIndicesRemaining = findall(x -> (x ∉ variablesToBeSubstituted), gens(polynomialRing))
    substRing, substRingVars = polynomial_ring(coefficient_ring(polynomialRing),
        symbols(polynomialRing)[variableIndicesRemaining])

    ###
    # Construct images of the variables under the substituion homomorphism
    ###
    # Step 1: construct images of the variables which are not substituted
    substImages = substRingVars
    for (i, x) in enumerate(gens(polynomialRing))
        if x in variablesToBeSubstituted
            insert!(substImages, i, polynomialRing(0))
        end
    end
    # Step 2: construct images of the variables which are substituted
    for (i, x) in enumerate(gens(polynomialRing))
        j = findfirst(isequal(x), variablesToBeSubstituted)
        if !isnothing(j)
            substImages[i] = evaluate(polynomialsToSubstitute[j], substImages)
        end
    end

    return hom(polynomialRing, substRing, substImages)

end

# Step 1: identify which variables can be substituted by which monomials
variablesToBeSubstituted = elem_type(B)[]
monomialsToSubstitute = elem_type(B)[]
for g in polynomialSystemParametrised
    if !is_binomial(g)
        continue
    end
    M = collect(monomials(g))
    C = collect(coefficients(g))

    # we are working in a polynomial ring for technical reasons
    # (lack of support for monomial orderings in Laurent polynomial rings)
    # so our simplification code only works if one of the two monomials is a variable
    @assert isone(total_degree(M[1])) || isone(total_degree(M[2]))

    if isone(total_degree(M[2])) && M[2] ∉ variablesToBeSubstituted
        push!(variablesToBeSubstituted, M[2])
        push!(monomialsToSubstitute, -C[1] / C[2] * M[1])
        continue
    elseif isone(total_degree(M[1])) && M[1] ∉ variablesToBeSubstituted
        push!(variablesToBeSubstituted, M[1])
        push!(monomialsToSubstitute, -C[2] / C[1] * M[2])
        continue
    end
end

# Step 2: make sure that substitutions are consistent,
#  i.e., that substituted monomials do not contain any variables that were substituted
for (x, s) in zip(variablesToBeSubstituted, monomialsToSubstitute)
    # eliminates x from monomialsToSubstitute
    global monomialsToSubstitute = evaluate.(monomialsToSubstitute,
        Ref([x]),
        Ref([s]))
end

# Step 3: create a new polynomial ring with the non-substituted variables
#   and construct the substitution homomorphism
phiSimplify = construct_substitution_homomorphism(B, variablesToBeSubstituted, monomialsToSubstitute)
Bsimplified = codomain(phiSimplify)
polynomialSystemSimplified = phiSimplify.(polynomialSystemParametrised)
polynomialSystemSimplified = polynomialSystemSimplified[findall(!iszero, polynomialSystemSimplified)]
