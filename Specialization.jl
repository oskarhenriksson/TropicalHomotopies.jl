

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