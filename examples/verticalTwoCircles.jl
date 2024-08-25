using Revise
using Oscar


###
# Choosing parameters (and their valuations)
###
K,t = rational_function_field(QQ,"t")
nu = tropical_semiring_map(K,t)

v = [1,0,0,1,0]
a = t .^v


###
# Constructing the tropicalizations
###
R,(x1,x2,y1,y2) = K["x1","x2","y1","y2"];
y3 = x1
y4 = x2
y5 = 1


Ilin = ideal(R,[a[1]*y1+a[2]*y2+a[3]*y3+a[4]*y4+a[5]*y5,
                3*a[1]*y1+3*a[2]*y2+5*a[3]*y3+7*a[4]*y4+11*a[5]*y5])
TropL = tropical_variety(Ilin,nu)

g1 = y1-x1^2
TropH1 = tropical_hypersurface(g1,nu)
g2 = y2-x2^2
TropH2 = tropical_hypersurface(g2,nu)


###
# verifying that the intersection is zero-dimensional and
# that only one pair of maximal polyhedra is intersecting
# the intersection necessarily has to be in their relative interiors
###
intersectionDimensions = [dim(reduce(intersect,[sigma,tau,rho])) for sigma in maximal_polyhedra(TropL), tau in maximal_polyhedra(TropH1), rho in maximal_polyhedra(TropH2)]
@assert isempty(findall(d->(d>0),intersectionDimensions)) "positive dimensional intersection"
@assert length(findall(d->(d>=0),intersectionDimensions))==1 "intersection potentially non-transversal"


###
# Computing the intersection point
###
TropI = reduce(stable_intersection,[TropH1,TropH2,TropL])
println(vertices(TropI))
println(multiplicities(TropI))


###
# Computing the initial form of the intersection point
###
S,(X1,X2) = K["X1","X2"]
I = ideal(S,[a[1]*X1^2+a[2]*X2^2+a[3]*X1+a[4]*X2+a[5],
             3*a[1]*X1^2+3*a[2]*X2^2+5*a[3]*X1+7*a[4]*X2+11*a[5]])
w = first(vertices(TropI))
display(initial.(gens(I),Ref(nu),Ref(w[1:2])))



###
# Extra: computing tropical Groebner basis of the homogenized ideal
###
S,(X0,X1,X2) = K["X0","X1","X2"]
Ih = ideal(S,[a[1]*X1^2+a[2]*X2^2+a[3]*X1*X0+a[4]*X2*X0+a[5]*X0^2,
              3*a[1]*X1^2+3*a[2]*X2^2+5*a[3]*X1*X0+7*a[4]*X2*X0+11*a[5]*X0^2])
Gh = groebner_basis(Ih,nu,vcat(0,w[1:2]))
G = evaluate.(Gh,Ref([1,X1,X2]))
