using Oscar
import HomotopyContinuation as HC


@doc raw"""
    tropically_generic_specialization(parametrizedPolynomialSystem::Vector{<:MPolyRingElem}; choice_of_parameters::Vector{<:RingElem}=nothing, check_genericity::Bool=true)

Specialize `parametrizedPolynomialSystem` at a random choice of parameters, or `choice_of_parameters` if specified. Return error if specialization is not tropically generic, if `check_genericity==true`.

!!! Warning
Requires `parametrizedPolynomialSystem` to be either linear or binomial.

```
"""
function tropically_generic_specialization(parametrizedPolynomialSystem::Vector{<:MPolyRingElem}; choice_of_parameters::Vector{<:RingElem}=nothing, check_genericity::Bool=true)

    if all(isequal(1),total_degree.(parametrizedPolynomialSystem))
        return tropically_generic_specialization_linear(parametrizedPolynomialSystem,
                                                        choice_of_parameter=choice_of_parameter,check_genericity=check_genericity)
    elseif all(isequal(2),length.(parametrizedPolynomialSystem))
        return tropically_generic_specialization_binomial(parametrizedPolynomialSystem,
                                                         choice_of_parameter=choice_of_parameter,check_genericity=check_genericity)
    else
        error("input unsupported (neither linear nor binomial)")
    end
end

function tropically_generic_specialization_linear(parametrizedLinearSystem::Vector{<:MPolyRingElem};
                                                  genericChoiceOfParameters::Vector{<:RingElem}=nothing,
                                                  check_genericity::Bool=true)

    Kax = parent(first(parametrizedLinearSystem))
    Ka = coefficient_ring(first(parametrizedLinearSystem))
    K = base_ring(Ka)
    parametrizedMacaulayMatrix = zero_matrix(Ka,length(parametrizedLinearSystem),ngens(Kax))
    for (i,f) in enumerate(parametrizedLinearSystem)
        for (c,xAlpha) in zip(coefficients(f),monomials(f))
            j = findfirst(isequal(xAlpha),gens(Kax))
            @assert !isnothing(j)
            parametrizedMacaulayMatrix[i,j] = c
        end
    end

    if isnothing(genericChoiceOfParameters)
        genericChoiceOfParameters = rand(Int,ngens(Ka))
    end

    if check_genericity 
        macaulayMatrix = matrix(K,[[evaluate(parametrizedMacaulayMatrix[i,j],genericChoiceOfParameters) for j in 1:ncols(parametrizedMacaulayMatrix)]
                                for i in 1:nrows(parametrizedMacaulayMatrix)])

        for I in AbstractAlgebra.combinations(ncols(macaulayMatrix),nrows(macaulayMatrix))
            if det(macaulayMatrix[:,I])==0
                @req det(parametrizedMacaulayMatrix[:,I])==0 "genericChoiceOfParameters not generic"
            end
        end
    end

    Kx,x = polynomial_ring(K,symbols(Kax))
    phi = hom(Kax,Kx,c->evaluate(c,genericChoiceOfParameters),x)
    linearSystem = phi.(parametrizedLinearSystem)
    return linearSystem
end

function tropically_generic_specialization_binomial(parametrizedBinomialSystem::Vector{<:MPolyRingElem};
                                                    choice_of_parameter::Vector{<:RingElem}=nothing,
                                                    check_genericity::Bool=true)

    Kax = polynomial_ring(first(parametrizedBinomialSystem))
    Ka = coefficient_ring(first(parametrizedBinomialSystem))
    K = base_ring(Ka)
    Kxy,xy = polynomial_ring(K,symbols(Kax))

    if isnothing(choice_of_parameters)
        choice_of_parameters = rand(Int,ngens(Ka))
    end
    phi = hom(Kax,Kxy,c->evaluate(c,choice_of_parameters),xy)
    binomialSystem = phi.(parametrizedBinomialSystem)

    if check_genericity
        @req all(isequal(2),length.(binomialSystem)) "choice of parameters not generic"
    end

    return binomialSystem
end



"""
    modify_vertically(polynomialSystem::Vector{<:MPolyRingElem})

    Carry out the vertical modification from [^HR23], where one introduces a new variable for each monomial appearing in the system. Outputs the linear system in the new variables, together with a binomial system.

    [^HS95]: https://browse.arxiv.org/abs/2206.07838

"""
function modify_vertically(polynomialSystem::Vector{<:MPolyRingElem}; keepLinearPolynomials::Bool=true)

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

    return linearSystem, binomialSystem

end


"""
    tropical_vertical_root_bound_with_homotopy_data(polynomialSystem::Vector{<:MPolyRingElem})

Returns the upper bound on the root count of `tropical_vertical_root_bound`, together with auxilarly data needed to construct start systems and homotopies for solving the system. This data consists of two tuples:
- A tuple consisting of the points and multiplicities of a tropical stable intersection, as well as a pertubation used for making the intersection finite.
- A vertical modification of `polynomialSystem`, consisting of a linear system and a binomial system.

To do: Exclude certain equations (use case: converation laws)

"""
function tropical_vertical_root_bound_with_homotopy_data(polynomialSystem::Vector{<:MPolyRingElem}; 
                keepLinearPolynomials::Bool=true)
                    
    linearSystem, binomialSystem = modify_vertically(polynomialSystem; keepLinearPolynomials=keepLinearPolynomials)

    # Todo: Faster to tropicalize separately?

    nu = tropical_semiring_map(QQ,min)
    TropL = tropical_variety_linear(ideal(linearSystem),nu)
    TropV = tropical_variety(ideal(binomialSystem),nu)[1] #todo: tropical_varities

    # Todo: Use the specialized function instead?
    intersection_points, multiplicities, pertubations = tropical_stable_intersection_after_perturbation([TropL,TropV])

    #intersection_points = [intersection_points[i,:] for i=1:size(intersection_points,2)]

    #pertubations = [zeros(QQ,length(pertubation)),[pertubation...]]

    return sum(multiplicities), (linearSystem,binomialSystem), pertubations, intersection_points, multiplicities
end

"""
    tropical_vertical_root_bound(polynomialSystem::Vector{<:MPolyRingElem})

Return an upper bound on the number of roots of `polynomialSystem` in the torus of the algebraic closure of the coefficient field (counted with multiplicity). This upper bound is always at most the mixed volume of the system. To obtain auxiliary data from the intermedaite tropical intersection that is computed, use `tropical_vertical_root_bound_with_homotopy_data`. 

# Examples
```jldoctest
julia> Qx, x = polynomial_ring(QQ,"x"=>1:2);

julia> F = [5*x[1]^3*x[2] - 6*x[1]*x[2]^3 + x[1]*x[2], 5*x[1]^3*x[2] - 6*x[1]*x[2]^3 - x[1]*x[2]^2];

julia> tropical_vertical_root_bound(F)
2

```
"""
function tropical_vertical_root_bound(polynomialSystem::Vector{<:MPolyRingElem})
    return tropical_root_bound_with_intersection_data(polynomialSystem)[1]
end





"""
    tropical_homotopies_from_homotopy_data(systems,pertubations,intersection_points)

Converts tropical intersection data (for instance  the auxilarly data from the Oscar function
`tropical_vertical_root_bound_with_homotopy_data`) into HomotopyContinuation.jl homotopies.

"""
function tropical_homotopies_from_homotopy_data(systems,pertubations,intersection_points)
    homotopies = tropical_homotopies_in_oscar_from_homotopy_data(systems,pertubations,intersection_points)
    return export_homotopy_from_oscar_to_HC.(homotopies)
end

function tropical_homotopies_in_oscar_from_homotopy_data(systems,pertubations,intersection_points; verbose=false)
    T = tropical_semiring(min)
    nu = tropical_semiring_map(QQ,min)
    Qx = parent(first(first(systems)))
    Qt, t = rational_function_field(QQ,"t") #todo: change to laurent_polynomial_ring when https://github.com/oscar-system/Oscar.jl/discussions/3156 is fixed
    Qtx, x = polynomial_ring(Qt, symbols(Qx))
    homotopies = [] #todo: add type
    for w in intersection_points
        homotopy_w = []
        if verbose 
            println("\nIntersection point:")
            display(w)
        end
        for (system,u) in zip(systems,pertubations)
            r = lcm(denominator.(Vector(w-u)))
            shift = Int.(Vector(r*Vector(w-u)))
            system_GB = groebner_basis_workaround(ideal(system),nu,shift)
            system_shifted = evaluate.(system_GB,Ref(t.^(shift) .* x))
            valuation_shift = -Int.(evaluate.(tropical_polynomial.(system_GB),Ref(shift)))
            homotopy_w = vcat(homotopy_w, (t.^valuation_shift).*system_shifted )
            if verbose
                println("\nSystem:")
                display(system)
                println("\nGröbner basis:")
                display(system_GB)
                println("\nShifted system:")
                display(system_shifted)
                println("\nValuation shift:")
                display(valuation_shift)
                println("\nHomotopy:")
                display((t.^valuation_shift).*system_shifted)
                println("\n")
            end
        end
        push!(homotopies,homotopy_w)
    end
    return homotopies
end


function export_homotopy_from_oscar_to_HC(homotopy)
    zipped_homotopy = []
    for h in homotopy
        exponent_vectors = collect(exponents(h))
        coefficient_values = map(c->Rational(coeff(numerator(c),degree(numerator(c)))),coefficients(h))
        t_exponents = degree.(numerator.(coefficients(h)))
        hzipped = zip(coefficient_values,exponent_vectors,t_exponents)
        push!(zipped_homotopy, hzipped)
    end
    n = ngens(parent(first(homotopy)))
    HC.@var t z[1:n]
    return HC.Homotopy( [ sum([ c*prod(z.^e)*t^a for (c,e,a) in hzipped]) for hzipped in zipped_homotopy] , z, t)
end

start_system(H::HC.ModelKit.Homotopy) = HC.System(HC.subs(H.expressions,H.t=>0))

function solve_binomial_system(F::HC.ModelKit.System)
    system_exponents, system_coefficients = HC.support_coefficients(F)
    @req all(isequal(2),length.(system_coefficients)) "Input system must be binomial"
    # Todo: Support for monomials?
    A = hcat( (E->E[:,1]-E[:,2]).(system_exponents)... )
    b = (c->-c[2]//c[1]).(system_coefficients)
    BSS = HC.BinomialSystemSolver(A,b)
    HC.solve!(BSS)
    return [Vector(BSS.X[:,i]) for i=1:size(BSS.X,2)]
end

"""
    vertical_solve(polynomial_system::Vector{<:MPolyRingElem})

Solve a polynomial system with homotopy continuation, by carrying out a vertical modification and using the vertical root bound.
"""
#todo: Distinguish between of original variables and modification variables
function vertical_solve(polynomial_system::Vector{<:MPolyRingElem}; verbose=false, detour=false)
    n = ngens(parent(first(polynomial_system)))
    _, systems, pertubations, intersection_points, multiplicities = 
            Oscar.tropical_vertical_root_bound_with_homotopy_data(polynomial_system)
    if verbose 
        display(systems)
        display(pertubations)
        display(intersection_points)
        display(multiplicities)
    end
    homotopies = tropical_homotopies_from_homotopy_data(systems,pertubations,intersection_points)
    verbose ? display(homotopies) : nothing
    all_solutions = Vector{ComplexF64}[]
    # Todo: Should this be its own function?
    for H in homotopies
        S = start_system(H)
        verbose ? display(S) : nothing
        start_solutions = solve_binomial_system(S)
        for s in start_solutions
            verbose ? display("\nNew start solution:") : nothing
            verbose ? display(s) : nothing
            if detour
                t1 = exp(2*pi*rand()*im)
                result1 = HC.track(HC.Tracker(H),s,0.0,t1)
                verbose ? display(result1) : nothing
                s1 = result1.solution
            else 
                t1 = 0.0
                s1 = s
            end
            result = HC.track(HC.Tracker(H),s1,t1,1)
            verbose ? display(result) : nothing
            verbose ? display(HC.subs(H.expressions,H.t=>1,H.variables=>result.solution)) : nothing
            push!(all_solutions,result.solution[1:n]) 
        end
    end
    return all_solutions
end


function HC_system_from_Oscar_system(polynomial_system)
    n = ngens(parent(first(polynomial_system)))
    HC.@var x[1:n]
    zipped_system = [  zip(Rational.(coefficients(f)),exponents(f)) for f in polynomial_system ]
    return HC.System([ sum([ c*prod(x.^e) for (c,e) in fzipped]) for fzipped in zipped_system], variables=x)
end

function parametric_HC_system_from_parametric_Oscar_system(polynomial_system)
    n = ngens(parent(first(polynomial_system)))
    m = ngens(coefficient_ring(first(polynomial_system)))
    HC.@var a[1:m] x[1:n]
    zip_coefficients = c -> zip(Rational.(coefficients(c)),exponents(c))
    zipped_system = [  zip(  zip_coefficients.(coefficients(f)) ,exponents(f)) for f in polynomial_system ]
    unzip_coefficient = czipped -> sum([d*prod(a.^e) for (d,e) in czipped])
    return HC.System([ sum([ unzip_coefficient(c)*prod(x.^e) for (c,e) in fzipped]) for fzipped in zipped_system], parameters=a, variables=x)
end
