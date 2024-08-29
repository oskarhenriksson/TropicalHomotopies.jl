# todo: add HC for commands from the HC package
# how does this fit with the rest of the code?


@doc raw"""
    polyhedral_homotopies_and_start_systems(F::HomotopyContinuation.ModelKit.System)

Compute polyhedral homotopies and start systems for a system of polynomial equations.
One for each mixed cell.

"""
function polyhedral_homotopies_and_start_systems(F::HomotopyContinuation.ModelKit.System)
    HC.@var t
    vars = HC.variables(F)
    n = length(vars)
    @assert n==length(F.expressions) "System needs to be square"
    supports = [transpose(HC.exponents_coefficients(f,vars)[1]) for f in F.expressions]
    supports = [hcat(eachrow(v)...) for v in supports]
    liftings = [rand(Int8,size(supports[i])[2]) for i=1:n]
    F_perturbed = []
    for j=1:n
        f = F.expressions[j]
        supps, coeffs = HC.exponents_coefficients(f,vars)
        append!(F_perturbed,sum([t^liftings[j][i]*coeffs[i]*prod((F.variables).^supps[:,i]) for i=1:size(supps)[2]]))
    end
    mc = HC.MixedSubdivisions.mixed_cells(supports, liftings)
    homotopies = []
    start_systems = []
    for cell in mc
        w = Rational.(cell.normal)
        N = lcm(map(x->x.den,w)) #for clearing denominators
        w = Int.(N*w)
        F_perturbed_reparametrized = HC.subs.(F_perturbed,Ref(t=>t^N))
        normalizations = Int.(-N*cell.β)
        Hw = [t^(normalizations[j])*HC.subs(F_perturbed_reparametrized[j],F.variables=>(t.^(w)).*(F.variables)) for j=1:length(F_perturbed)]
        Hw = HC.expand.(Hw)
        Sw = HC.subs(Hw,t=>0)
        append!(homotopies,[HC.Homotopy(Hw,vars,t)])
        append!(start_systems,[HC.System(Sw)])
    end
    return homotopies, start_systems
end