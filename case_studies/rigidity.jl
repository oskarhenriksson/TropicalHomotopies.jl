using Oscar

################################################################################
#
#  Computing Functions
#
################################################################################

function atomic_macaulay_matrices(G::Graph, d::Int)
    ###
    # Constructing atomic Macaulay Matrix of W:
    #   ( 1    1        ... 1 )
    #     |    |            |
    #     w_ij u_{ij,1}     u_{ij,d}
    # To assemble: take a product, one per edge (ij)
    ###
    macW = matrix(QQ,ones(Int,1,d+1))

    ###
    # Constructing atomic Macaulay Matrix of L:
    #   ( I_{|E(G)|}         | I_G )
    #      | ... |              |
    #     columns for u_ij,k   columns for x_{i,k}
    # where
    # - I_{|E(G)|} identity matrix with a row/column for every edge in G
    # - I_G 0/±1 matrix encoding the edges of G, one column per vertex of G, one row per edge of G
    # To assemble: take a product, one per dimension
    ###
    I1 = identity_matrix(QQ,ne(G))
    I2 = matrix(QQ,zeros(QQ,ne(G),nv(G)))
    for (i,edge) in enumerate(edges(G))
        I2[i,dst(edge)] = 1
        I2[i,src(edge)] = -1
    end
    macL = hcat(I1,I2)

    return macW, macL
end


function matroid_from_macaulay_matrix(linearEquationsMatrix::QQMatrix)
    ###
    # Compute the projection matrix onto the (orthogonal) complement of the linear space,
    # or in other words the images of the unit vectors projected onto the complement
    ###
    basisOfComplementTransposed = nullspace(linearEquationsMatrix)[2]
    basisOfComplement = transpose(basisOfComplementTransposed)
    projectionMatrix = basisOfComplementTransposed*inv(basisOfComplement*basisOfComplementTransposed)*basisOfComplement
    return matroid_from_matrix_columns(projectionMatrix)
end



################################################################################
#
#  Graph Constructors
#
################################################################################

function triangles_graph(n::Int)
    @req n>=3 "unsupported number of vertices"
    G = Graph{Undirected}(n);
    # subdivided square
    add_edge!(G,1,2);
    add_edge!(G,2,3);
    add_edge!(G,1,3);
    for i in 4:n
        add_edge!(G,i,i-1)
        add_edge!(G,i,i-2)
    end
    return G
end

function laman_graph(n::Int)
    if n==6
        G = Graph{Undirected}(6);
        # top triangle
        add_edge!(G,1,2);
        add_edge!(G,2,3);
        add_edge!(G,3,1);
        # bottom triangle
        add_edge!(G,4,5);
        add_edge!(G,5,6);
        add_edge!(G,6,4);
        # other edges
        add_edge!(G,1,4);
        add_edge!(G,2,5);
        add_edge!(G,3,6);
        return G
    elseif n==7
        G = Graph{Undirected}(7);
        # top left rectangle
        add_edge!(G,1,2);
        add_edge!(G,2,3);
        add_edge!(G,3,4);
        add_edge!(G,4,1);
        # top right rectangle
        add_edge!(G,3,4);
        add_edge!(G,4,5);
        add_edge!(G,5,6);
        add_edge!(G,6,3);
        # bottom left triangle
        add_edge!(G,1,2);
        add_edge!(G,2,7);
        add_edge!(G,7,1);
        # bottom right triangle
        add_edge!(G,5,6);
        add_edge!(G,6,7);
        add_edge!(G,7,5);
        return G
    end
    @req false "unsupported number of vertices"
end

function geiringer_graph(n::Int)
    if n==6
        G = Graph{Undirected}(6);
        # outer edges of hexagon
        add_edge!(G,1,2);
        add_edge!(G,2,3);
        add_edge!(G,3,4);
        add_edge!(G,4,5);
        add_edge!(G,5,6);
        add_edge!(G,6,1);
        # odd triangle
        add_edge!(G,1,3);
        add_edge!(G,3,5);
        add_edge!(G,5,1);
        # even triangle
        add_edge!(G,2,4);
        add_edge!(G,4,6);
        add_edge!(G,6,2);
        return G
    elseif n==7
        G = Graph{Undirected}(7);
        # edges of pentagon
        add_edge!(G,1,2);
        add_edge!(G,2,3);
        add_edge!(G,3,4);
        add_edge!(G,4,5);
        add_edge!(G,5,1);
        # connecting pentagon to vertex 6 on top
        add_edge!(G,1,6);
        add_edge!(G,2,6);
        add_edge!(G,3,6);
        add_edge!(G,4,5);
        add_edge!(G,5,6);
        # connecting pentagon to vertex 7 on bottom
        add_edge!(G,1,7);
        add_edge!(G,2,7);
        add_edge!(G,3,7);
        add_edge!(G,4,7);
        add_edge!(G,5,7);
        return G
    end
    @req false "unsupported number of vertices"
end

import Oscar.cone
function cone(G::Graph)
    edgesG = Vector{Int}.(edges(G))
    n = n_vertices(G)
    edgesConeG = vcat(edgesG,[[i,n+1] for i in 1:n])
    return graph_from_edges(edgesConeG)
end

import Oscar.tropical_linear_space
function tropical_linear_space(G::Graph, nu::Union{Nothing,TropicalSemiringMap}=nothing)
    M = zero_matrix(QQ,nv(G),ne(G))
    for (i,edge) in enumerate(edges(G))
        M[src(edge),i] = 1
        M[dst(edge),i] = -1
    end
    return tropical_linear_space(M,nu)
end


import Base.*
function multiply_by_scalar(u::PointVector{QQFieldElem}, c::QQFieldElem)
    return u .* c
end
function multiply_by_scalar(u::RayVector{QQFieldElem}, ::QQFieldElem)
    return u
end
function *(c::QQFieldElem, Sigma::PolyhedralComplex)
    SigmaVertsAndRays = vertices_and_rays(Sigma)
    SigmaRayIndices = findall(vr -> vr isa RayVector, SigmaVertsAndRays)
    SigmaLineality = lineality_space(Sigma)
    SigmaIncidence = maximal_polyhedra(IncidenceMatrix,Sigma)
    return polyhedral_complex(SigmaIncidence, multiply_by_scalar.(SigmaVertsAndRays,c), SigmaRayIndices, SigmaLineality)
end
*(c::ZZRingElem, Sigma::PolyhedralComplex) = QQ(c)*Sigma
*(c::Rational, Sigma::PolyhedralComplex) = QQ(c)*Sigma
*(c::Int, Sigma::PolyhedralComplex) = QQ(c)*Sigma

*(Sigma::PolyhedralComplex,c::QQFieldElem) = c*Sigma
*(Sigma::PolyhedralComplex,c::ZZRingElem) = QQ(c)*Sigma
*(Sigma::PolyhedralComplex,c::Rational) = QQ(c)*Sigma
*(Sigma::PolyhedralComplex,c::Int) = QQ(c)*Sigma


function *(c::ZZRingElem, TropV::Oscar.TropicalVarietySupertype)
    Sigma = polyhedral_complex(TropV)
    mults = multiplicities(TropV)
    minOrMax = convention(TropV)
    return tropical_variety(c*Sigma,abs(c)*mults,minOrMax)
end
*(c::Int, TropV::Oscar.TropicalVarietySupertype) = ZZ(c)*TropV

*(TropV::Oscar.TropicalVarietySupertype,c::ZZRingElem) = c*TropV
*(TropV::Oscar.TropicalVarietySupertype,c::Int) = ZZ(c)*TropV


import Base.-
function -(Sigma::PolyhedralComplex)
    SigmaVertsAndRays = vertices_and_rays(Sigma)
    SigmaRayIndices = findall(vr -> vr isa RayVector, SigmaVertsAndRays)
    SigmaLineality = lineality_space(Sigma)
    SigmaIncidence = maximal_polyhedra(IncidenceMatrix,Sigma)
    return polyhedral_complex(SigmaIncidence, -SigmaVertsAndRays, SigmaRayIndices, SigmaLineality)
end
function -(TropV::Oscar.TropicalVarietySupertype)
    Sigma = polyhedral_complex(TropV)
    mults = multiplicities(TropV)
    minOrMax = convention(TropV)
    return tropical_variety(-Sigma,mults,minOrMax)
end


import Base.+
function translate_by_vector(u::PointVector{QQFieldElem}, v::Vector{QQFieldElem})
    return u .+ v
end
function translate_by_vector(u::RayVector{QQFieldElem}, ::Vector{QQFieldElem})
    return u
end
function +(v::Vector{QQFieldElem}, Sigma::PolyhedralComplex)
    @req length(v)==ambient_dim(Sigma) "ambient dimension mismatch"
    SigmaVertsAndRays = vertices_and_rays(Sigma)
    SigmaRayIndices = findall(vr -> vr isa RayVector, SigmaVertsAndRays)
    SigmaLineality = lineality_space(Sigma)
    SigmaIncidence = maximal_polyhedra(IncidenceMatrix,Sigma)
    return polyhedral_complex(SigmaIncidence, translate_by_vector.(SigmaVertsAndRays,Ref(v)), SigmaRayIndices, SigmaLineality)
end
+(v::Vector{ZZRingElem}, Sigma::PolyhedralComplex) = QQ.(v)+Sigma
+(v::Vector{Rational}, Sigma::PolyhedralComplex) = QQ.(v)+Sigma
+(v::Vector{Int}, Sigma::PolyhedralComplex) = QQ.(v)+Sigma

+(Sigma::PolyhedralComplex, v::Vector{QQFieldElem}) = v+Sigma
+(Sigma::PolyhedralComplex, v::Vector{ZZRingElem}) = QQ.(v)+Sigma
+(Sigma::PolyhedralComplex, v::Vector{Rational}) = QQ.(v)+Sigma
+(Sigma::PolyhedralComplex, v::Vector{Int}) = QQ.(v)+Sigma

function +(v::Vector{QQFieldElem}, TropV::Oscar.TropicalVarietySupertype)
    Sigma = polyhedral_complex(TropV)
    mults = multiplicities(TropV)
    minOrMax = convention(TropV)
    return tropical_variety(v+Sigma,mults,minOrMax)
end
+(v::Vector{ZZRingElem}, TropV::Oscar.TropicalVarietySupertype) = QQ.(v)+TropV
+(v::Vector{Rational}, TropV::Oscar.TropicalVarietySupertype) = QQ.(v)+TropV
+(v::Vector{Int}, TropV::Oscar.TropicalVarietySupertype) = QQ.(v)+TropV

+(TropV::Oscar.TropicalVarietySupertype, v::Vector{QQFieldElem}) = v+TropV
+(TropV::Oscar.TropicalVarietySupertype, v::Vector{ZZRingElem}) = QQ.(v)+TropV
+(TropV::Oscar.TropicalVarietySupertype, v::Vector{Rational}) = QQ.(v)+TropV
+(TropV::Oscar.TropicalVarietySupertype, v::Vector{Int}) = QQ.(v)+TropV


import Oscar.direct_sum
function direct_sum(Sigma1::PolyhedralComplex, Sigma2::PolyhedralComplex)

    n1 = ambient_dim(Sigma1)
    n2 = ambient_dim(Sigma2)

    ###
    # Construct vertices and rays of the direct sum
    ###
    Sigma1VertsAndRays = vertices_and_rays(Sigma1)
    Sigma2VertsAndRays = vertices_and_rays(Sigma2)
    Sigma1VertsAndRaysPadded = vcat.(Sigma1VertsAndRays, Ref(zeros(QQ,n2)))
    Sigma2VertsAndRaysPadded = vcat.(Ref(zeros(QQ,n1)), Sigma2VertsAndRays)
    Sigma12VertsAndRays = vcat(Sigma1VertsAndRaysPadded, Sigma2VertsAndRaysPadded)

    ###
    # Construct the ray indices of the direct sum
    ###
    Sigma1RayIndices = findall(vr -> vr isa RayVector, Sigma1VertsAndRays)
    Sigma2RayIndices = findall(vr -> vr isa RayVector, Sigma2VertsAndRays)
    Sigma2RayIndicesShifted = Sigma2RayIndices .+ length(Sigma1VertsAndRays)
    Sigma12RayIndices = vcat(Sigma1RayIndices, Sigma2RayIndicesShifted)

    ###
    # Construct the lineality space of the direct sum
    ###
    Sigma1Lineality = lineality_space(Sigma1)
    Sigma2Lineality = lineality_space(Sigma2)
    Sigma1LinealityPadded = vcat.(Sigma1Lineality, Ref(zeros(QQ,n2)))
    Sigma2LinealityPadded = vcat.(Ref(zeros(QQ,n1)), Sigma2Lineality)
    Sigma12Lineality = vcat(Sigma1LinealityPadded, Sigma2LinealityPadded)

    ###
    # Construct the incidence matrix of the direct sum
    ###
    Sigma1Incidence = maximal_polyhedra(IncidenceMatrix,Sigma1)
    Sigma2Incidence = maximal_polyhedra(IncidenceMatrix,Sigma2)
    Sigma1Intcidence = [findall(Sigma1Incidence[i,:]) for i in 1:nrows(Sigma1Incidence)]
    Sigma2IntcidenceShifted = [findall(Sigma1Incidence[i,:]) .+ ncols(Sigma1Incidence) for i in 1:nrows(Sigma2Incidence)]
    Sigma12Intcidence = Vector{Int}[]
    for i1 in Sigma1Intcidence
        for i2 in Sigma2IntcidenceShifted
            push!(Sigma12Intcidence, vcat(i1,i2))
        end
    end
    Sigma12Incidence = IncidenceMatrix(Sigma12Intcidence)

    Sigma12 = polyhedral_complex(Sigma12Incidence, Sigma12VertsAndRays, Sigma12RayIndices, Sigma12Lineality)
    return Sigma12

end

function direct_sum(TropL1::TropicalLinearSpace{minOrMax,true}, TropL2::TropicalLinearSpace{minOrMax,true}) where {minOrMax <: Union{typeof(min),typeof(max)}}
    Sigma12 = direct_sum(polyhedral_complex(TropL1), polyhedral_complex(TropL2))
    mults = ZZ.(ones(n_maximal_polyhedra(Sigma12)))
    return TropicalLinearSpace{minOrMax,true}(Sigma12,mults)
end


import Oscar.⊕
⊕(Sigma1::PolyhedralComplex, Sigma2::PolyhedralComplex) = direct_sum(Sigma1,Sigma2)
⊕(TropL1::TropicalLinearSpace, TropL2::TropicalLinearSpace) = direct_sum(TropL1,TropL2)

contains_in_interior(v::AbstractVector, P::Polyhedron) =
  Polymake.polytope.contains_in_interior(Oscar.pm_object(P), coefficient_field(P).([1; v]))::Bool

function is_intersection_transverse(sigma1::Polyhedron, sigma2::Polyhedron, sigma12::Polyhedron)
    # check that sigma12 intersects the interior of sigma1 and sigma2
    p = relative_interior_point(sigma12)
    if !contains_in_interior(p,sigma1) || !contains_in_interior(p,sigma2)
        return false
    end
    # check that sigma1 and sigma2 span the entire space
    L1 = kernel(affine_equation_matrix(affine_hull(sigma1))[:,2:end]; side=:right)
    L2 = kernel(affine_equation_matrix(affine_hull(sigma2))[:,2:end]; side=:right)
    if rank(hcat(L1,L2))<ambient_dim(sigma1)
        return false
    end
    return true
end

function stable_intersection_transverse(TropV1::Oscar.TropicalVarietySupertype, TropV2::Oscar.TropicalVarietySupertype)
    Sigma12 = Polyhedron{QQFieldElem}[]
    mults12 = ZZRingElem[]

    for (j1,(sigma1,m1)) in enumerate(maximal_polyhedra_and_multiplicities(TropV1))
        for (j2,(sigma2,m2)) in enumerate(maximal_polyhedra_and_multiplicities(TropV2))
            println("intersecting (", j1, ",", j2, ") of (", n_maximal_polyhedra(TropV1), ",", n_maximal_polyhedra(TropV2),"), current mults ",mults12)
            sigma12 = intersect(sigma1, sigma2)
            if dim(sigma12)<0
                continue
            end

            @req is_intersection_transverse(sigma1, sigma2, sigma12) "intersection is not transverse"

            i = findfirst(isequal(sigma12), Sigma12)
            if isnothing(i)
                push!(Sigma12, sigma12)
                push!(mults12, m1*m2*Oscar.tropical_intersection_multiplicity(sigma1,sigma2))
            else
                mults12[i] += m1*m2*Oscar.tropical_intersection_multiplicity(sigma1,sigma2)
            end
        end
    end

    if isempty(Sigma12)
        # return empty tropical variety
        Sigma = polyhedral_complex(IncidenceMatrix(),zero_matrix(QQ,0,ambient_dim(TropV1)))
        mults = ZZRingElem[]
        return tropical_variety(Sigma,mults,convention(TropV1))
    end

    return tropical_variety(Sigma12,mults12,convention(TropV1))
end


function chain_of_flats_differences(v::PointVector{QQFieldElem})
    E = sortperm(v, rev=true)
    chainOfFlatsDifferences = Vector{Int}[]
    currentFlatsDifference = Int[E[1]]
    for i in 2:length(E)
        if v[E[i]]==v[E[i-1]]
            push!(currentFlatsDifference,E[i])
        else
            push!(chainOfFlatsDifferences,currentFlatsDifference)
            currentFlatsDifference = Int[E[i]]
        end
    end
    push!(chainOfFlatsDifferences,currentFlatsDifference)
    return chainOfFlatsDifferences
end


function permute_entries(vv::Vector{Vector{Int}}, perm::Vector{Int})
    return [[perm[i] for i in v] for v in vv]
end

# ################################################################################
# #
# #  Example code
# #
# ################################################################################

# G = geiringer_graph(7)
# macW, macL = atomic_macaulay_matrices(G,3)
# matL = matroid_from_macaulay_matrix(macL)
# transL = Oscar.Polymake.matroid.check_transversality(matL.pm_matroid)
