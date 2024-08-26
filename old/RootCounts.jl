
"""
    linear_and_binomial_part(polynomialSystem::Vector{<:MPolyRingElem})

    Carry out the vertical modification from [^HR23], where one introduces a new variable for each monomial appearing in the system. Outputs the linear system in the new variables, together with a binomial system.

    [^HS95]: https://browse.arxiv.org/abs/2206.07838

"""
function linear_and_binomial_part(polynomialSystem::Vector{<:MPolyRingElem}; keepLinearPolynomials::Bool=true)

    # Construct a list of distinct monomials (sorted for the sake of consistency)
    Kx = parent(first(polynomialSystem))

    if keepLinearPolynomials
        polynomialsNotToBeModified = [f for f in polynomialSystem if total_degree(f)==1]
        polynomialsToBeModified = [f for f in polynomialSystem if total_degree(f)>1]
    else
        polynomialsNotToBeModified = []
        polynomialsToBeModified = polynomialSystem
    end

    distinctMonomials = sort(unique(collect(Iterators.flatten(monomials.(polynomialsToBeModified)))))

    # Construct a new polynomial ring, with an extra variable per distinct monomial
    K = coefficient_ring(Kx)
    Kxy, x, y = polynomial_ring(K, symbols(Kx),"y"=>1:length(distinctMonomials))

    linearSystem = elem_type(Kxy)[]
    for f in polynomialsToBeModified
        # substitute monomial xAlpha by the corresponding w
        push!(linearSystem,sum([c*y[findfirst(isequal(xAlpha),distinctMonomials)] for (c,xAlpha) in zip(coefficients(f),monomials(f))]))
    end

    append!(linearSystem, evaluate.(polynomialsNotToBeModified,Ref(x)))

    distinctMonomials = evaluate.(distinctMonomials,Ref(x))
    binomialSystem = [y[i]-mon for (i,mon) in enumerate(distinctMonomials)]

    return [linearSystem, binomialSystem]

end



function tropical_stable_intersection_after_perturbation(Sigmas::Vector{<:Oscar.TropicalVarietySupertype{typeof(min), true}};
    randomDirections::Union{Vector{Vector{Int}},Nothing}=nothing, verbose=false)

    if isnothing(randomDirections)
        n = ambient_dim(first(Sigmas))
        randomDirections = vcat([zeros(Int,n)],[rand(Int8,n) for i=2:length(Sigmas)])
    end

    shiftedMaximalPolyhedra = [maximal_polyhedra(Sigma).+Ref(v) 
                                    for (Sigma,v) in zip(Sigmas,randomDirections)]

    if verbose >= 3
        println("Nunber of intersections to consider: ",prod(length.(shiftedMaximalPolyhedra)))
    end

    originalMultiplicities = [multiplicities(Sigma) for Sigma in Sigmas]

    stableIntersectionMultiplicites = Int[]
    stableIntersectionPoints = Vector{QQFieldElem}[] # Vector instead of PointVector

    for (sigmas,mults) in zip(Iterators.product(shiftedMaximalPolyhedra...),Iterators.product(originalMultiplicities...))
        sigma = intersect(sigmas...)
        @req dim(sigma)<=0 "randomDirection not generic (perturbed varieties are overlapping)"
        if dim(sigma)==0
            vertex = first(vertices(sigma))
            for s in sigmas
                @req all(is_negative, affine_inequality_matrix(facets(s))*vcat([1],vertex)) "randomDirection not generic (lower-dimensional skeleta are intersecting)"
            end
            push!(stableIntersectionPoints, Vector(vertex)) #convert from PointVector to Vector to enable componentwise operations later
            mult = prod(mults)
            for i in 2:length(sigmas)
                mult *= Oscar.tropical_intersection_multiplicity(intersect(sigmas[1:(i-1)]...),sigmas[i])
            end
            push!(stableIntersectionMultiplicites,mult)
        end
    end
            
    return  stableIntersectionPoints, stableIntersectionMultiplicites, randomDirections
end


function tropical_root_bound_with_homotopy_data(systems::Vector{Vector{QQMPolyRingElem}}; verbose=false, print_result = false)
    @req !isempty(systems) "systems must be a nonempty list"
    @req all(isequal(parent(systems[1][1])), parent.(vcat(systems...))) "all systems must have the same parent ring"
    nu = tropical_semiring_map(QQ,min)
    time_trop = @elapsed tropicalizations = [tropical_variety(ideal(system),nu)[1] for system in systems]
    if verbose >= 2
        println("Time spent tropicalizing: ", time_trop)
    end
    time_intersecting = @elapsed intersection_points, multiplicities, pertubations = tropical_stable_intersection_after_perturbation(tropicalizations, verbose=verbose)
    if verbose >= 3
        println("Time spent on stable intersection: ", time_intersecting)
        #println("Multiplicities: ", multiplicities)
    end
    if print_result
        println("Intersection points:")
        println(intersection_points)
        println()
        println("Multiplicities: ")
        println(multiplicities)
        println()
        println("Pertubations:")
        println(pertubations)
    end
    return sum(multiplicities), pertubations, intersection_points, multiplicities
end



"""
    tropical_root_bound(systems::Vector{Vector{QQMPolyRingElem}})

    Given a partition of a square system into subsystems, compute an upper bound on the number of solutions in the torus, 
    based on the stable intersection of the tropicalizations.
    
"""
tropical_root_bound(systems::Vector{Vector{QQMPolyRingElem}}) = tropical_root_bound_with_homotopy_data(systems)[1]

