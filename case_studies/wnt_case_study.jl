include("../src/main.jl")

# Set up the parmetric WNT pathway system

include("construct_wnt_system.jl")
system = polynomialSystemSimplified
number_of_parameters = ngens(coefficient_ring(parent(first(system))))

# Make a choice of parameters
target_parameters = collect(1:number_of_parameters)
target_system = specialize(polynomialSystemSimplified, target_parameters)

# Embedding in vertical family
F, target_parameters = vertical_embedding(target_system)

# TODO: Explicitly comptue TropL and TropB


# Computation of tropical data
# Redraw the perturbed parameters until we get a transverse intersection
m = length(target_parameters)
Kt, t = rational_function_field(QQ, "t")

v = rand(-100:100, m)
perturbed_parameters = (t.^v) .* target_parameters
grc, projected_pts, initial_systems, tropical_groebner_bases, perturbed_parameters = 
            tropical_root_count_with_homotopy_data_vertical(F, perturbed_parameters=perturbed_parameters, verbose=true);
            
# TODO: Include this 'try' trick in the general code
grc = nothing
projected_pts = nothing
initial_systems = nothing
tropical_groebner_bases = nothing
perturbed_parameters = nothing
for attempts = 1:5
    v = rand(-100:100, m)
    perturbed_parameters = (t.^v) .* target_parameters
    try 
        grc, projected_pts, initial_systems, tropical_groebner_bases, perturbed_parameters = 
            tropical_root_count_with_homotopy_data_vertical(F, perturbed_parameters=perturbed_parameters, verbose=true);
        break
    catch e
        println(e)
        continue
    end
end


# Computation of homotopies 
homotopies = [  homotopy_from_tropical_data(G,w) 
        for (w, G) in zip(projected_pts, tropical_groebner_bases)  ]  

# Print out one of the homotopies for tracing
H = homotopies[1]
S = initial_systems[1]
println(H)
println(S)


target_parameters = QQ.(rand(-10:10,m))

# BUG: Initials not correct for 
# target_parameters = [-81//25, 613//100, -47//5, 46//25, 86//25, 61//25, -7//50, 467//100, 743//100, -169//100, -42//5, -39//25, -69//100, -387//50, 493//50, -707//100, 893//100, -209//100, 329//100]
tropical_solve(F, target_parameters, type_of_system=:vertical)
