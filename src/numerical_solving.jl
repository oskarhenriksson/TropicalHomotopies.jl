import HomotopyContinuation as HC

@doc raw"""
    tropical_solve(F::Vector{<:PolyRingElem},target_parameters::Vector{QQFieldElem},type_of_system::Symbol;kwargs...)

Solve a system of polynomial equations with our tropical techniques. 
The input consists of a parametrized system `F`, a vector of target parameters `target_parameters`, 
and the type of system `type_of_system`. 

The `target_parameters` should be generic and rational. 
(One can use a parameter homotopy to trace the solutions to a specific choice of parameters later.) 

!!! note
    As of now, `type_of_system`` has to be `:vertical`. Implementing `:horizontal`` is work in progress.


"""
function tropical_solve(F::Vector{<:MPolyRingElem}, target_parameters::Vector{QQFieldElem}; type_of_system::Symbol, verbose::Bool=false, kwargs...)

    # Perturb parameters
    m = length(target_parameters)
    Kt, t = rational_function_field(QQ, "t")

    # Compute tropical data and start solutions
    # This part of the code is specific to each type of system
    # Repick pertubation until the data is computed successfully
    v = nothing
    perturbed_parameters = nothing
    grc = nothing
    projected_pts = nothing
    initial_systems = nothing
    tropical_groebner_bases = nothing
    perturbed_parameters = nothing
    start_solutions_for_homotopies = nothing
    for attempts = 1:10
        v = rand(-100:100, m)
        perturbed_parameters = (t .^ v) .* target_parameters

        try
            if type_of_system == :vertical
                grc, projected_pts, initial_systems, tropical_groebner_bases, perturbed_parameters =
                    tropical_root_count_with_homotopy_data_vertical(F, perturbed_parameters=perturbed_parameters, verbose=verbose)
                time_start_systems = 0
                start_solutions_for_homotopies = []
                for S in initial_systems
                    time_start_systems += @elapsed S_HC = HC_system_from_Oscar_system(S)
                    start_solutions = solve_binomial_system(S_HC)
                    push!(start_solutions_for_homotopies, start_solutions)
                end
                if verbose
                    println("Time spent computing start solutions: ", time_start_systems)
                end
            else
                error("Tropical data not implemented for $type_of_system")
            end
            break
        catch e
            println(e)
            continue
        end
    end

    # Compute homotopies from the tropical data
    homotopies = [homotopy_from_tropical_data(G, w)
                  for (w, G) in zip(projected_pts, tropical_groebner_bases)]

    all_solutions = Vector{ComplexF64}[]
    time_tracing = 0
    for (H, start_solutions) in zip(homotopies, start_solutions_for_homotopies)
        # Ensure square homotopy
        if length(H) > ngens(parent(first(H)))
            A = rand(-100:100, ngens(parent(first(H))), length(H))
            H = A * H
        elseif length(H) < ngens(parent(first(H)))
            error("Underdetermined homotopy detected")
        end
        # Export to HC format
        H_HC = export_homotopy_from_oscar_to_HC(H)
        # Trace the solutions along the homotopy
        time_tracing += @elapsed new_solutions = trace_solutions_along_homotopy(H_HC, start_solutions)
        append!(all_solutions, new_solutions)
    end
    if verbose
        println("Time spent tracing solutions: ", time_tracing)
    end
    return all_solutions
end

function trace_solutions_along_homotopy(H_HC::HC.ModelKit.Homotopy, start_solutions)
    t1 = exp(2 * pi * im * rand())
    res1 = HC.track.(HC.Tracker(H_HC), start_solutions, 0, t1)
    sol1 = [HC.solution(r) for r in res1 if r.return_code == :success]
    res = HC.track.(HC.Tracker(H_HC), sol1, t1, 1)
    new_solutions = [HC.solution(r) for r in res if r.return_code == :success]
    return new_solutions
end


"""
    solve_binomial_system(F::HC.ModelKit.System)

Solve a binomial system of equations (i.e., a system where each polynomial has two terms)
using the `BinomialSystemSolver` from the `HomotopyContinuation.jl` package.

"""
function solve_binomial_system(F::HC.ModelKit.System)
    system_exponents, system_coefficients = HC.support_coefficients(F)
    @assert all(isequal(2), length.(system_coefficients)) "Input system must be binomial"
    A = hcat((E -> E[:, 1] - E[:, 2]).(system_exponents)...)
    b = (c -> -c[2] // c[1]).(system_coefficients)
    BSS = HC.BinomialSystemSolver(A, b)
    HC.solve!(BSS)
    return [Vector(BSS.X[:, i]) for i = 1:size(BSS.X, 2)]
end


####################################
# Interface between HC and OSCAR
####################################

"""
    export_homotopy_from_oscar_to_HC(homotopy)

Convert a homotopy from OSCAR (in the form of a system over the rational function field Q(t))
to a`Homotopy` object in the `HomotopyContinuation.jl` package.

The variables in the output will always be named `x`, and the time variable will always be named `t`.

"""
function export_homotopy_from_oscar_to_HC(homotopy)
    zipped_homotopy = []
    for h in homotopy
        exponent_vectors = collect(exponents(h))
        coefficient_values = map(c -> Rational(coeff(numerator(c), degree(numerator(c)))), coefficients(h))
        t_exponents = degree.(numerator.(coefficients(h)))
        hzipped = zip(coefficient_values, exponent_vectors, t_exponents)
        push!(zipped_homotopy, hzipped)
    end
    n = ngens(parent(first(homotopy)))
    HC.@var t x[1:n]
    return HC.Homotopy([sum([c * prod(x .^ e) * t^a for (c, e, a) in hzipped]) for hzipped in zipped_homotopy], x, t)
end



"""
    HC_system_from_Oscar_system(polynomial_system)

Convert OSCAR system (in the form of a system over Q) to a `System` object in the `HomotopyContinuation.jl` package.    

The variables in the output will always be named `x`.
"""
function HC_system_from_Oscar_system(polynomial_system)
    n = ngens(parent(first(polynomial_system)))
    HC.@var x[1:n]
    zipped_system = [zip(Rational.(coefficients(f)), exponents(f)) for f in polynomial_system]
    return HC.System([sum([c * prod(x .^ e) for (c, e) in fzipped]) for fzipped in zipped_system], variables=x)
end


"""
    parametric_HC_system_from_parametric_Oscar_system(polynomial_system)

Convert a parametric system (in the form of a system over Q[a]) to a parametrized `System` object in the `HomotopyContinuation.jl` package.

The parameters in the output will always be named `a` and the variables will always be named `x`.

"""
function parametric_HC_system_from_parametric_Oscar_system(polynomial_system)
    n = ngens(parent(first(polynomial_system)))
    m = ngens(coefficient_ring(first(polynomial_system)))
    HC.@var a[1:m] x[1:n]
    zip_coefficients = c -> zip(Rational.(coefficients(c)), exponents(c))
    zipped_system = [zip(zip_coefficients.(coefficients(f)), exponents(f)) for f in polynomial_system]
    unzip_coefficient = czipped -> sum([d * prod(a .^ e) for (d, e) in czipped])
    return HC.System([sum([unzip_coefficient(c) * prod(x .^ e) for (c, e) in fzipped]) for fzipped in zipped_system], parameters=a, variables=x)
end

function parametric_Oscar_system_from_parametric_HC_system(F::HC.ModelKit.System)
    params = HC.parameters(F)
    A, a = polynomial_ring(QQ, string.(params))
    vars = HC.variables(F)
    R, x = polynomial_ring(A, string.(vars))
    polynomials = HC.expressions(F)
    polynomials_oscar = MPolyRingElem[]
    for f in polynomials
        exps, coeffs = HC.exponents_coefficients(f, vars)
        zip_coeff = c -> zip(eachcol(HC.exponents_coefficients(c, params)[1]),HC.exponents_coefficients(c, params)[2])
        fzipped = zip(eachcol(exps),zip_coeff.(coeffs))
        unzip_coeff = czipped -> sum([d * prod(a .^ e) for (e,d) in czipped])
        f_oscar = sum([unzip_coeff(c) * prod(x .^ e) for (e, c) in fzipped])
        push!(polynomials_oscar, f_oscar)
    end
    return polynomials_oscar
end