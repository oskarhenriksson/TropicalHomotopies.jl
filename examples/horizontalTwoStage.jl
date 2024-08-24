using Revise
using Oscar

###
# Choosing parameters (and their valuations)
###
v = [4,2,0,0,11,7,0,1]
a = [1,1,1,1,2,3,5,7] .* t.^v


###
# Constructing the tropicalizations
###
K,t = rational_function_field(QQ,"t")
nu = tropical_semiring_map(K,t)

R,(x1,x2,y1) = K["x1","x2","y1"];
y2 = x1
z1 = y1^3
z2 = y1^2
z3 = y2
z4 = 1

f1 = a[1]*z1+a[2]*z2+a[3]*z3+a[4]*z4
f2 = a[5]*z1+a[6]*z2+a[7]*z3+a[8]*z4
h = y1 - (1+x1+x2)

TropH1 = tropical_hypersurface(f1,nu)
TropH2 = tropical_hypersurface(f2,nu)
TropH3 = tropical_hypersurface(h,nu)


###
# verifying that the Newton subdivisions are maximal.
# this means that the initial forms w.r.t. weight vectors
# in the relative interior of the maximal polyhedra
# will be binomial
###
@assert n_maximal_polyhedra.([TropH1,TropH2,TropH3])==[5,5,6] "Newton subdivisions not maximal"

intersectionDimensions = [dim(reduce(intersect,[sigma,tau,rho])) for sigma in maximal_polyhedra(TropH1), tau in maximal_polyhedra(TropH2), rho in maximal_polyhedra(TropH3)]


###
# Computing the intersection point
###
TropI = reduce(stable_intersection,[TropH1,TropH2,TropH3])


###
# Computing the initial forms of the intersection points
###
for w in vertices(TropI)
    display(initial.(gens(I),Ref(nu),Ref(w)))
end



###
# Extra: constructing the homotopies
###
R,(x1,x2,y1,y2,z1,z2,z3,z4) = K["x1","x2","y1","y2","z1","z2","z3","z4"];

v = [4,2,0,0,11,7,0,1]
a = [1,1,1,1,2,3,5,7] .* t.^v

f1 = a[1]*z1+a[2]*z2+a[3]*z3+a[4]*z4
f2 = a[5]*z1+a[6]*z2+a[7]*z3+a[8]*z4

g1 = z1 - y1^3
g2 = z2 - y1^2
g3 = z3 - y2
g4 = z4 - 1

h1 = y1 - (1+x1+x2)
h2 = y2 - x1

F = [f1,f2,g1,g2,g3,g4,h1,h2]

w = [1,-1,-1,1,-3,-2,1,0]

Hw = evaluate.(F,Ref(gens(R) .* t.^w))

evaluate.(Hw,Ref([x1,x2,(t+t^2*x1+x2),x1,(t+t^2*x1+x2)^3,(t+t^2*x1+x2)^2,x1,1]))
