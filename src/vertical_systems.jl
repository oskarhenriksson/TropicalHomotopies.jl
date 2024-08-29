

@doc raw"""
    vertical_embedding(F::Vector{<:MPolyRingElem})

Embedd a polynomial system over `QQ` into a vertically parametrized family in the canonical way.

The variables in the output will be the same as in the input system.

The parameters will be a[1],...,a[m], where m is the number of distinct monomials in the input system.

"""
function vertical_embedding(F::Vector{<:MPolyRingElem})
    coefficient_matrix, distinct_monomials = coefficient_matrix_and_monomials(F)
    m = length(distinct_monomials)
    Kx = parent(first(F))
    K = coefficient_ring(Kx)
    Ka, a = polynomial_ring(K, "a"=>1:m)
    Kax, x = polynomial_ring(Ka,symbols(Kx))
    embedding = hom(Kx,Kax,x)
    return coefficient_matrix*(a.*embedding.(distinct_monomials)), [K(1) for i=1:m]
end


@doc raw"""
    modify_vertically(F::Vector{<:MPolyRingElem})

Modify a vertically parametrized system into a linear part and a binomial part.

All the parameters will be in the binomial part.

The modification variables will be denoted by y[1],...,y[m], where m is the number of distinct monomials in the input system.

"""
function modify_vertically(F::Vector{<:MPolyRingElem})

    Kax = parent(first(F))
    Ka = coefficient_ring(Kax)
    K = base_ring(Ka)
    Kax_joint, a,x = polynomial_ring(K,symbols(Ka),symbols(Kax)) 

    make_parameters_variables = hom(Kax,Kax_joint,hom(Ka,Kax_joint,a),x)

    C, scaled_monomials = coefficient_matrix_and_monomials(make_parameters_variables.(F))

    for parameter in a
        @req sum([divides(monomial,parameter)[1] for monomial in scaled_monomials]) <= 1 "System is not vertically parametrized (parameter $parameter appears in more than one monomial)"
    end

    Kaxy, x,y = polynomial_ring(Ka,symbols(Kax),"y"=>1:length(scaled_monomials))
    embedding = hom(Kax_joint,Kaxy,vcat(Kaxy.(gens(Ka)),x))


    linear_part = C*y
    binomial_part = y-embedding.(scaled_monomials)


    return linear_part, binomial_part
end

@doc raw"""
    tropical_root_count_with_homotopy_data_vertical(F; perturbed_parameters=nothing)

Compute tropical intersection data for a vertically parametrized system through a vertical modification.

If no choice of parameters `perturbed_parameters` is given, random ones will be chosen.

The output will be a tuple consisting of:
1. The generic root count
2. The points in the tropical stable intersection
3. The initial system for each point
4. The tropical Gröbner bases for each point
5. The choice of `perturbed_parameters` used for the computation.

!!! warning
    If the choice of `perturbed_parameters` is not generic enough, an error will be thrown.
    In this case, you should try again with a different choice of `perturbed_parameters`.

"""
function tropical_root_count_with_homotopy_data_vertical(F; perturbed_parameters=nothing, verbose::Bool=false)

    linear_part, binomial_part = modify_vertically(F)

    m = length(binomial_part) # number of scaled monomials

    Kaxy = parent(first(binomial_part))
    Ka = coefficient_ring(Kaxy)
    K = base_ring(Ka)
    Kx, x = polynomial_ring(K, symbols(parent(first(F))))

    monomial_vector = -hom(Kaxy,Kx, hom(Ka,K,ones(Int,m)),vcat(gens(Kx),zeros(Kx,m))).(binomial_part)

    # Tropicalize the linear part over QQ
    linear_part_specialized = specialize(linear_part, K.(ones(Int,m)))
    nu_K = tropical_semiring_map(K)
    time_TropL = @elapsed TropL = tropical_linear_space(ideal(linear_part_specialized),nu_K)

    if verbose
        println("Time spent tropicalizing the linear part: ", time_TropL)
        println("Number of maximal polyhedra: ", length(maximal_polyhedra(TropL)))
    end

    # If no choice of perturbed parameters is given, choose random ones
    if isnothing(perturbed_parameters)
        v = rand(-100:100, m)
        Kt, t = rational_function_field(K, "t")
        perturbed_parameters = t.^v
    else
        Kt = parent(first(perturbed_parameters))
        t = gens(Kt)[1]
    end

    Ktx, x = polynomial_ring(Kt, symbols(Kx))
    Ktxy, xy = polynomial_ring(Kt, symbols(Kaxy))

    binomial_part_specialized = 
        hom(Kaxy,Ktxy,hom(Ka,Kt,perturbed_parameters),gens(Ktxy)).(binomial_part)
    
    # Tropicalize the binomial part over QQt
    nu = tropical_semiring_map(Kt,t)
    time_TropB = @elapsed TropB = Oscar.tropical_variety_binomial(ideal(binomial_part_specialized),nu)

    if verbose
        println("Time spent tropicalizing the binomial part: ", time_TropB)
    end

    # Intersect the TropL and TropB
    time_intersection = @elapsed pts, mults = tropical_stable_intersection_linear_binomial(TropL,TropB)

    if verbose
        println("Time spent on stable intersection: ", time_intersection)
    end

    # Project points to the original variable space
    projected_pts = [w[1:ngens(Kx)] for w in pts] 

    # Substitution homomorphism Kxy -> Kx with y -> target_parametrs.*monomial_vector
    Kxy, xy = polynomial_ring(K,symbols(Kaxy))
    target_parameters = (c->evaluate(c,QQ(1))).(perturbed_parameters)
    target_monomial_vector = target_parameters .* monomial_vector
    substitute_y_by_monomials = hom(Kxy,Kx,vcat(gens(Kx),target_monomial_vector)) 

    # Substitution homomorphism Kxy -> Ktx with y -> perturbed_parameters.*monomial_vector
    perturbed_monomial_vector = perturbed_parameters .* hom(Kx,Ktx,gens(Ktx)).(monomial_vector)
    substitute_y_by_perturbed_monomials = hom(Kxy,Ktx,vcat(gens(Ktx),perturbed_monomial_vector)) 

    # Compute tropical Gröbner bases and initial systems
    tropical_groebner_bases = Vector{MPolyRingElem}[]
    initial_systems = Vector{QQMPolyRingElem}[]
    for w in pts
        # Compute the tropical Gröbner basis and initial for the linear part
        G_linear = groebner_basis(ideal(linear_part_specialized),nu_K,w)
        initials_linear = initial.(G_linear, Ref(nu_K), Ref(w))
        
        # Substitute the y variables by the monomials
        initials = substitute_y_by_monomials.(initials_linear)
        @req all(isequal(2),length.(initials)) "Non-binomial initial detected"

        # Alternative computation of the initials
        G = substitute_y_by_perturbed_monomials.(G_linear)
        initials_alternative = initial.(G,Ref(nu),Ref(w[1:ngens(Kx)]))
        @req initials_alternative == initials "Initials are not correct: $(initial.(G,Ref(nu),Ref(w[1:ngens(Kx)]))) != $initials"
        
        push!(tropical_groebner_bases,G)
        push!(initial_systems,initials)
    end

    return sum(mults), projected_pts, initial_systems, tropical_groebner_bases, perturbed_parameters
end