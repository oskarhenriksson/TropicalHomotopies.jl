
"""
    vertical_solve(polynomialSystem::Vector{<:MPolyRingElem}; transversality_check=true)

Solve a polynomial system via the tropical homotopy associated to a vertical decomposition.

If `transversality_check` is true, we check if the coefficient matrix of the linear part has a
transversal column matroid. In that case, the linear part is replaced by a new linear part that
generate the same ideal, in such a way that the mixed volume of the new system is the generic root 
count of the original system. The new system is solved with the built-in tool for 
polyhedral homotopies in HomotopyContinuation.jl.
"""
function vertical_solve(polynomialSystem::Vector{<:MPolyRingElem}; transversality_check=true)
    n = ngens(parent(first(polynomialSystem)))
    decomposition = linear_and_binomial_part(polynomialSystem)
    linear_part = decomposition[1]
    binomial_part = decomposition[2]
    xy = gens(parent(first(linear_part)))
    if transversal_check
        A = linear_coefficient_matrix(linear_part)
        tp = transversal_presentation(A)
    else
        tp = false
    end
    if tp == false
        solutions = tropical_solve(decomposition)
    else
        Anew = realize_supports_from_rowspace(A,tp)
        new_linear_part = Anew*xy
        new_polynomial_system = vcat(new_linear_part,binomial_part)
        solutions = HC.solutions(HC.solve(HC_system_from_Oscar_system(new_polynomial_system)))
    end
    return [s[1:n] for s in solutions]
end


"""
    transversal_presentation(A::QQMatrix)

Checks if the column matroid of a rational matrix `A` is transversal.
If it is, a transversal presentation of the matroid is returned; if it is not, `false` is returned.

"""
function transversal_presentation(A::QQMatrix)
    M = matroid_from_matrix_columns(A;check=false).pm_matroid;
    transversality_witness = Polymake.matroid.check_transversality(M)
    if transversality_witness == false
        return false
    else
        return [w.+1 for w in transversality_witness]
    end
end


"""
    linear_coefficient_matrix(linearSystem::Vector{<:MPolyRingElem})

Return the matrix of coefficients of a linear polynomial system.
"""
function linear_coefficient_matrix(linearSystem::Vector{<:MPolyRingElem})
    @assert all(isequal(1),total_degree.(lp)) "Input system must be linear"
    R = parent(first(linearSystem))
    K = coefficient_ring(R)
    variables = gens(R)
    rows = [transpose(coeff.(f,variables)) for f in linearSystem]
    return matrix(K,vcat(rows...))
end


"""
    realize_supports_from_rowspace(A,supports)

Given a matrix `A` and a list of supports, return a matrix `B` with the same row space as `A` such that
the supports of the rows of `B` are given by `supports`. Returns an error if impossible.

"""
function realize_supports_from_rowspace(A,supports)
    d,n = size(A)
    @assert length(supports)==rank(A) "Number of supports must equal rank(A)"
    complement_of_row_space = transpose(A)*inv(A*transpose(A))*A-identity_matrix(QQ,n)
    B = zero_matrix(QQ,length(supports),n)
    for (i,supp) in enumerate(supports)
        submatrix = complement_of_row_space[:,supp]
        nullspace_basis = nullspace(submatrix)[2]
        vector_in_nullspace = nullspace_basis*rand(-100:100,size(nullspace_basis,2))
        B[i,supp] = vector_in_nullspace
    end
    @req rref(A)[2]==rref(B)[2] "Failed to obtain same row space"
    return B
end
