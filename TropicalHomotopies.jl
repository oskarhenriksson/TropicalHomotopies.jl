module TropicalHomotopies

using Oscar
using HomotopyContinuation

export tropical_homotopies
export tropical_start_systems_with_solutions
export tropical_solve

function groebner_basis_linear_workaround(I::MPolyIdeal,nu::TropicalSemiringMap,w::Vector{QQFieldElem})
    F = gens(I)
    R = parent(first(F))
    xand1 = vcat(gens(R),[one(R)])
    wand0 = tropical_semiring(nu).(vcat(w,[0]))
    xand1order = xand1[sortperm(wand0)]
    macaulayMatrix = matrix(coefficient_ring(R),[[coeff(f,xj) for xj in xand1order] for f in F])
    return rref(macaulayMatrix)[2]*xand1order
end

function groebner_basis_binomial_workaround(I::MPolyIdeal,nu::TropicalSemiringMap,w::Vector{QQFieldElem})
    F = gens(I)
    R = parent(first(F))
    G = groebner_basis(ideal(F),complete_reduction=true)
    # Uncommment when https://github.com/oscar-system/Oscar.jl/issues/3277 is fixed
    #inG = [initial(g,nu,w) for g in G]
    #if min(length.(inG))==1
    #    return [one(R)] 
    #end
    A = matrix(ZZ,[first(collect(expv))-last(collect(expv)) for expv in exponents.(F)])
    b = [ -QQ(last(collect(coeff))/first(collect(coeff))) for coeff in coefficients.(F)]
    S,T,_ = snf_with_transform(A)
    A_basis = T*A
    b_basis = [prod(QQ(b[i])^T[j,i] for i=1:ncols(T)) for j=1:nrows(T)]
    G = elem_type(R)[]
    for i = 1:nrows(A_basis)
        alpha = A_basis[i,:]
        if is_zero(alpha) 
            if b_basis[i]!=1
                return [one(R)]
            end
        else
            push!(G,prod(x^a for (a,x) in zip(alpha,gens(R)) if a>0)-b_basis[i]*prod(x^(-a) for (a,x) in zip(alpha,gens(R)) if a<0))
        end
    end
    return G
end

function groebner_basis_workaround(I::MPolyIdeal,nu::TropicalSemiringMap,w::Vector{QQFieldElem})
    @req is_valuation_trivial(nu) "Unsupported valuation (must be trivial)"
    F = gens(I)
    if length(F)==1
        return F
    elseif all(isequal(1),total_degree.(F)) #linear
        return groebner_basis_linear_workaround(I,nu,w)
    elseif all(isequal(2),length.(F)) #binomial case
        return groebner_basis_binomial_workaround(I,nu,w)
    else
        error("Unsupported ideal (must be principal, binomial or linear)")
    end
end




function tropical_homotopies(polynomialSystem::Vector{<:MPolyRingElem})

    _, (pts,mults,u), (linearSystem,binomialSystem) = 
                 tropical_vertical_root_bound_with_auxilary_data(polynomialSystem)

    homotopies = construct_homotopies_in_oscar(systems,u,pts)

    return export_homotopy_from_oscar_to_HC.(homotopies))

end


"""
    construct_homotopies_in_oscar(systems,pertubations,intersection_points)

TBW
"""
function construct_homotopies_in_oscar(systems,pertubations,intersection_points)
    T = tropical_semiring(min)
    nu = tropical_semiring_map(QQ,min)
    Qx = parent(first(first(systems)))
    Qt, t = rational_function_field(QQ,"t") #todo: change to laurent_polynomial_ring when https://github.com/oscar-system/Oscar.jl/discussions/3156 is fixed
    Qtx, x = PolynomialRing(Qt, symbols(Qx))
    homotopies = []
    for w in intersection_points
        display(w)
        homotopy_w = []
        for (system,u) in zip(systems,pertubations)
            display(u)
            r = lcm(denominator.(Vector(w-u)))
            shift = Int.(Vector(r*Vector(w-u)))
            system_GB = groebner_basis_workaround(ideal(system),nu,shift)
            system_shifted = evaluate.(system_GB,Ref(t.^(shift) .* x))
            valuation_shift = -Int.(evaluate.(tropical_polynomial.(system_GB),Ref(shift)))
            homotopy_w = vcat(homotopy_w, (t.^valuation_shift).*system_shifted )
        end
        push!(homotopies,homotopy_w)
    end
    return homotopies
end

function export_homotopy_from_oscar_to_HC(homotopy)
    Qtx = parent(first(homotopy))
    Qt = coefficient_ring(Qtx)
    outputRing, outputVars = polynomial_ring(QQ,vcat(symbols(Qt),symbols(Qtx)))
    outputHom = hom(Qtx, outputRing, c->evaluate(c,outputVars[1]), outputVars[2:end]  )
    outputHomotopy = outputHom.(homotopy)
    zipped_homotopy = [  zip(Rational.(coefficients(h)),exponents(h)[2:end],exponents(h)[1])  for h in outputHomotopy ]
    n = length(first(zipped_homotopy)[2])
    @var t z[1:n]
    return Homotopy( [ sum([ c*prod(z.^e)*t^a for (c,e,a) in hzipped]) for hzipped in output] , z, t)
end


"""
    start_system_with_solutions(H)

TBW
"""
function start_system_with_solutions(H)
    startSystem =  subs(H.expressions,z[1]=>0)
    startSystemExponents, startSystemCoefficients  = HomotopyContinuation.support_coefficients(System(startSystem))
    A = hcat( (E->E[:,1]-E[:,2]).(startSystemExponents)... )
    b = (c->-c[2]//c[1]).(startSystemCoefficients)
    BSS = HomotopyContinuation.BinomialSystemSolver(A,b)
    HomotopyContinuation.solve!(BSS)
    solutions = [Vector(BSS.X[:,i]) for i=1:size(BSS.X,2)]
    return startSystem, solutions
end



"""
    tropical_solve(polynomialSystem::Vector{<:MPolyRingElem})

TBW
"""
function tropical_solve(polynomialSystem::Vector{<:MPolyRingElem})
    oscar_homotopies = construct_homotopies_in_oscar(polynomialSystem)
    HC_solutions = Vector{ComplexF64}[]
    for oscar_homotopy in oscar_homotopies
        HC_homotopy = export_homotopy_from_oscar_to_HC(oscar_homotopy)
        S = start_system_with_solutions(HC_homotopy)[2]
        for s in S 
            result = HomotopyContinuation.track(Tracker(HC_homotopy),s,0.0,1.0)
            push!(HC_solutions,result.solution)
        end
    end
    #todo: project the solutions to the original coordinates
    return HC_solutions
end

end # module