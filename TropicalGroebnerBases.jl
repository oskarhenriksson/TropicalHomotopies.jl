# Todo: Replace with new groebner_basis function in Oscar 
function groebner_basis_linear_workaround(I::MPolyIdeal,nu::TropicalSemiringMap,w::Vector{QQFieldElem})
    F = gens(I)
    R = parent(first(F))
    xand1 = vcat(gens(R),[one(R)])
    wand0 = tropical_semiring(nu).(vcat(w,[0]))
    xand1order = xand1[sortperm(wand0)]
    macaulayMatrix = matrix(coefficient_ring(R),[[coeff(f,xj) for xj in xand1order] for f in F])
    rowReducedMacaulayMatrix = rref(macaulayMatrix)[2]
    return rowReducedMacaulayMatrix*xand1order
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
    A = matrix(ZZ,[first(collect(expv))-last(collect(expv)) for expv in Oscar.exponents.(F)])
    b = [ -QQ(last(collect(coeff))/first(collect(coeff))) for coeff in Oscar.coefficients.(F)]
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

groebner_basis_workaround(I::MPolyIdeal,nu::TropicalSemiringMap,w::Vector{<:Union{Int,ZZRingElem}}) = 
groebner_basis_workaround(I,nu,QQ.(w))


function groebner_basis_workaround(I::MPolyIdeal,nu::TropicalSemiringMap,w::Vector{QQFieldElem})
    #@req is_trivial_semiring(nu) "Unsupported valuation (must be trivial)"
    F = gens(I)
    if length(F)==1
        return F
    elseif all(isequal(1),total_degree.(F)) #linear
        return groebner_basis_linear_workaround(I,nu,w)
    elseif all(isequal(2),length.(F)) #binomial case
        return gens(I)
        #return groebner_basis_binomial_workaround(I,nu,w)
    else
        error("Unsupported ideal (must be principal, binomial or linear)")
    end
end

