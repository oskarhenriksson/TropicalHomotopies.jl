using TropicalHomotopies, Oscar

using Random
Random.seed!(1234)
## Construction of the system
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

# Make a choice of parameters
number_of_parameters = ngens(coefficient_ring(parent(first(system))))
#target_parameters = collect(1:number_of_parameters)
target_parameters = rand(-100:100, number_of_parameters)
target_system = specialize(system, target_parameters)
## Embedding in vertical family
F, target_parameters = vertical_embedding(target_system)
## Pertubation of parameters
m = length(target_parameters)
Kt, t = rational_function_field(QQ, "t")

v = [-2, 5, -51, -98, -81, 73, 88, 94, 60, 30, 31, -95, 26, 91, -13, -31, -83, 42, 70, -19, 5, 54, 23, -92, 31, -26, -19]
#v = rand(-100:100, m)

perturbed_parameters = (t .^ v) .* target_parameters;

linear_part, binomial_part = modify_vertically(F)
display(linear_part)
display(binomial_part)
## Tropicalizations
Kaxy = parent(first(binomial_part));
Ka = coefficient_ring(Kaxy);
K = base_ring(Ka);
Kx, x = polynomial_ring(K, symbols(parent(first(F))));
# Tropicalize the linear part over QQ
linear_part_specialized = specialize(linear_part, K.(ones(Int, ngens(Ka))))
nu_K = tropical_semiring_map(K)
@time TropL = tropical_linear_space(ideal(linear_part_specialized), nu_K)
length(maximal_polyhedra(TropL))
lineality_space(TropL)
Ktx, x = polynomial_ring(Kt, symbols(Kx))
Ktxy, xy = polynomial_ring(Kt, symbols(Kaxy))

binomial_part_specialized =
    hom(Kaxy, Ktxy, hom(Ka, Kt, perturbed_parameters), gens(Ktxy)).(binomial_part)

# Tropicalize the binomial part over QQt
nu = tropical_semiring_map(Kt, t)
@time TropB = Oscar.tropical_variety_binomial(ideal(binomial_part_specialized), nu)
## Computation of tropical intersection points
@time pts, mults = tropical_stable_intersection_linear_binomial(TropL, TropB)
projected_pts = [w[1:ngens(Kx)] for w in pts]
grc = sum(mults)

display(projected_pts)
display(mults)
display(grc)
## Computation of initials and homotopies for a given point
# Substitution homomorphism Kxy -> Kx with y -> target_parametrs.*monomial_vector
Kxy, xy = polynomial_ring(K, symbols(Kaxy))
target_parameters = (c -> evaluate(c, QQ(1))).(perturbed_parameters)
monomial_vector = -hom(Kaxy, Kx, hom(Ka, K, ones(Int, m)), vcat(gens(Kx), zeros(Kx, m))).(binomial_part)
target_monomial_vector = target_parameters .* monomial_vector
substitute_y_by_monomials = hom(Kxy, Kx, vcat(gens(Kx), target_monomial_vector))

# Substitution homomorphism Kxy -> Ktx with y -> perturbed_parameters.*monomial_vector
perturbed_monomial_vector = perturbed_parameters .* hom(Kx, Ktx, gens(Ktx)).(monomial_vector)
substitute_y_by_perturbed_monomials = hom(Kxy, Ktx, vcat(gens(Ktx), perturbed_monomial_vector));

# Choice of tropical point
#w = pts[1] # this works okay
w = pts[2] # this gives nonbinomial initials
#w = pts[6] # this works okay

# Compute the tropical Gröbner basis and initial for the linear part
G_linear = groebner_basis(ideal(linear_part_specialized), nu_K, w)
initials_linear = initial.(G_linear, Ref(nu_K), Ref(w))

# Substitute the y variables by the monomials
initials = substitute_y_by_monomials.(initials_linear)

for f in initials
   @req terms(f) |> collect |> length == 2 "Nonbinomial initial obtained: $(f)"
end

# Tropical Gröbner basis
G = substitute_y_by_perturbed_monomials.(G_linear)
# Homotopy in OSCAR format
projected_point = w[1:ngens(Kx)]
H = homotopy_from_tropical_data(G,projected_point)
## Solve start system and trace along homotopy
# Solve start system
S_HC = HC_system_from_Oscar_system(initials)
start_solutions = solve_binomial_system(S_HC)
# Homotopy in HC format
H_HC = export_homotopy_from_oscar_to_HC(H)
trace_solutions_along_homotopy(H_HC, start_solutions)