include("../src/main.jl")

# Set up the parmetric WNT pathway system

# Define ring of parameters
A, (k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k31, c1, c2, c3, c4, c5) =
    polynomial_ring(QQ, ["k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11", "k12", "k13", "k14", "k15", "k16", "k17", "k18", "k19", "k20", "k21", "k22", "k23", "k24", "k25", "k26", "k27", "k28", "k29", "k30", "k31", "c1", "c2", "c3", "c4", "c5"]);

# Define ring of parameterized polynomials
B, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19) = polynomial_ring(A, ["x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "x15", "x16", "x17", "x18", "x19"]);

# Define the square(!) system of steady state equations and conservation laws
steadyStateEqs = [
    -k1 * x1 + k2 * x2, # x1
    k1 * x1 - (k2 + k26) * x2 + k27 * x3 - k3 * x2 * x4 + (k4 + k5) * x14, # x2
    k26 * x2 - k27 * x3 - k14 * x3 * x6 + (k15 + k16) * x15, # x3
    -k3 * x2 * x4 - k9 * x4 * x10 + k4 * x14 + k8 * x16 + (k10 + k11) * x18, # x4
    -k28 * x5 + k29 * x7 - k6 * x5 * x8 + k5 * x14 + k7 * x16, # x5
    -k14 * x3 * x6 - k20 * x6 * x11 + k15 * x15 + k19 * x17 + (k21 + k22) * x19, # x6
    k28 * x5 - k29 * x7 - k17 * x7 * x9 + k16 * x15 + k18 * x17, # x7
    -k6 * x5 * x8 + (k7 + k8) * x16, # x8
    -k17 * x7 * x9 + (k18 + k19) * x17, # x9
    k12 - (k13 + k30) * x10 - k9 * x4 * x10 + k31 * x11 + k10 * x18, # x10
    -k23 * x11 + k30 * x10 - k31 * x11 - k20 * x6 * x11 - k24 * x11 * x12 + k25 * x13 + k21 * x19, # x11
    -k24 * x11 * x12 + k25 * x13, # x12
    -k24 * x11 * x12 + k25 * x13, # x13
    k3 * x2 * x4 - (k4 + k5) * x14, # x14
    k14 * x3 * x6 - (k15 + k16) * x15, # x15
    -k6 * x5 * x8 + (k7 + k8) * x16, # x16
    -k17 * x7 * x9 + (k18 + k19) * x17, # x17
    k9 * x4 * x10 - (k10 + k11) * x18, # x18
    k20 * x6 * x11 - (k21 + k22) * x19]; # x19

conservationLaws = [
    (x1 + x2 + x3 + x14 + x15) - c1,
    (x4 + x5 + x6 + x7 + x14 + x15 + x16 + x17 + x18 + x19) - c2,
    (x8 + x16) - c3,
    (x9 + x17) - c4,
    (x12 + x13) - c5];

system = vcat([steadyStateEq for (i, steadyStateEq) in enumerate(steadyStateEqs) if i ∉ [3, 4, 8, 9, 12]],
    conservationLaws)

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

# Solve the system multiple times and check what the root counts are
for i=1:5
        sols = tropical_solve(F, target_parameters; type_of_system=:vertical, verbose=true)
        println(length(sols))
end


F_target_HC = HC_system_from_Oscar_system(target_system)
HC.certify(F_target_HC,sols)