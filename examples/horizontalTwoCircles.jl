using Revise
using Oscar


###
# Choosing parameters (and their valuations)
###
v = [3,0,2,3,2,2,3,2] # first choice
# v = [0,1,2,0,0,1,2,0] # second choice
b = t .^v


###
# Constructing the tropicalizations
###
K,t = rational_function_field(QQ,"t")
nu = tropical_semiring_map(K,t)

R,(x1,x2,y1) = K["x1","x2","y1"];
y2 = x1
y3 = x2
y4 = 1


Ilin = ideal(R,[b[1]*y1+b[2]*y2+b[3]*y3+b[4]*y4, 3*b[5]*y1+5*b[6]*y2+7*b[7]*y3+11*b[8]*y4])
TropL = tropical_variety(Ilin,nu)

f = x1^2+x2^2+y1
TropH = tropical_hypersurface(f,nu)


###
# verifying that the intersection is zero-dimensional and
# that only one pair of maximal polyhedra is intersecting
# the intersection necessarily has to be in their relative interiors
###
intersectionDimensions = [dim(intersect(sigma,tau)) for sigma in maximal_polyhedra(TropL), tau in maximal_polyhedra(TropH)]
@assert isempty(findall(d->(d>0),isequal(1),intersectionDimensions)) "positive dimensional intersection"
@assert length(findall(d->(d>=0),intersectionDimensions))==1 "intersection potentially non-transversal"


###
# Computing the intersection point
###
TropI = stable_intersection(TropL,TropH)
println(vertices(TropI))
println(multiplicities(TropI))


###
# Computing the initial form of the intersection point
###
S,(X1,X2) = K["X1","X2"]
I = ideal(S,[b1*X1^2+b1*X2^2+b2*X1+b3*X2+b4,
             3*b5*X1^2+3*b5*X2^2+5*b6*X1+7*b7*X2+11*b8])
w = first(vertices(TropI))
display(initial.(gens(I),Ref(nu),Ref(w[1:2])))



###
# Extra: computing tropical Groebner basis of the homogenized ideal
###
S,(X0,X1,X2) = K["X0","X1","X2"]
I = ideal(S,[b1*X1^2+b1*X2^2+b2*X1*X0+b3*X2*X0+b4*X0^2, 3*b5*X1^2+3*b5*X2^2+5*b6*X1*X0+7*b7*X2*X0+11*b8*X0^2])
initial(I,nu,vcat(0,w[1:2]))
