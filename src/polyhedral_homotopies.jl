@doc raw"""
    polyhedral_homotopies_and_start_systems(F::HomotopyContinuation.ModelKit.System)

Compute polyhedral homotopies and start systems for a system of polynomial equations in the format
of the `HomotopyContinuation.jl` package. This is done by computing a mixed subdivision for a randomly
chosen lifting. We obtain one polyhedral homotopy and one start system for each mixed cell.

"""
function polyhedral_homotopies_and_start_systems(F::HC.ModelKit.System)
    HC.@var t
    vars = HC.variables(F)
    n = length(vars)
    @assert n == length(F.expressions) "System needs to be square"

    # Compute the supports for each polynomial in the system
    supports = [transpose(HC.exponents_coefficients(f, vars)[1]) for f in F.expressions]

    # Convert the supports from LinearAlgebra.Transpose{Int, Matrix{Int}} to Matrix{Int}
    supports = [hcat(eachrow(v)...) for v in supports]

    # Make a choice of liftings for all terms in in the system
    liftings = [rand(Int8, size(supports[i])[2]) for i = 1:n]

    # Compute the mixed cells
    mc = HC.MixedSubdivisions.mixed_cells(supports, liftings)

    # Given the choice of liftings, compute perturbed polynomials 
    F_perturbed = []
    for j = 1:n
        f = F.expressions[j]
        supps, coeffs = HC.exponents_coefficients(f, vars)
        append!(F_perturbed, sum([t^liftings[j][i] * coeffs[i] * prod((F.variables) .^ supps[:, i]) for i = 1:size(supps)[2]]))
    end

    # For each mixed cell, compute the polyhedral homotopies and start systems
    homotopies = HC.ModelKit.Homotopy[]
    start_systems = HC.ModelKit.System[]
    for cell in mc

        # Tropical intersection point
        w = rationalize.(cell.normal, tol=1e-6)
        
        # Clear denominators
        N = lcm(map(x -> x.den, w))
        w = Int.(N * w)
        F_perturbed_reparametrized = HC.subs.(F_perturbed, Ref(t => t^N))
        normalizations = Int.(-N * rationalize.(cell.β, tol=1e-6))

        # Compute polyhedral homotopy
        Hw = [t^(normalizations[j]) * HC.subs(F_perturbed_reparametrized[j], F.variables => (t .^ (w)) .* (F.variables)) for j = 1:length(F_perturbed)]
        Hw = HC.expand.(Hw)
        append!(homotopies, [HC.Homotopy(Hw, vars, t)])

        # Obtain start system by setting t=0
        Sw = HC.subs(Hw, t => 0)
        append!(start_systems, [HC.System(Sw)])
    end
    return homotopies, start_systems
end