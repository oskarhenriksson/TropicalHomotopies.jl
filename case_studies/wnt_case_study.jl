include("../src/main.jl")

# Set up the parmetric WNT pathway system

include("construct_wnt_system.jl")
#system = polynomialSystemSimplified
system = polynomialSystemParametrised
number_of_parameters = ngens(coefficient_ring(parent(first(system))))

# Make a choice of parameters
target_parameters = collect(1:number_of_parameters)
target_system = specialize(system, target_parameters)

# Embedding in vertical family
F, target_parameters = vertical_embedding(target_system)

# Computation of tropical data
m = length(target_parameters)
Kt, t = rational_function_field(QQ, "t")

v = rand(-100:100, m)
perturbed_parameters = (t.^v) .* target_parameters
grc, projected_pts, initial_systems, tropical_groebner_bases, perturbed_parameters = 
            tropical_root_count_with_homotopy_data_vertical(F, perturbed_parameters=perturbed_parameters, verbose=true);
            

# Computation of homotopies 
homotopies = [  homotopy_from_tropical_data(G,w) 
        for (w, G) in zip(projected_pts, tropical_groebner_bases)  ]  

# Print out one of the homotopies
i = 1
H = homotopies[i]
S = initial_systems[i]
println(H)
println(S)


for i=1:5
        sols = tropical_solve(F, target_parameters; type_of_system=:vertical, verbose=true)
        println(length(sols))
end