# Extension of multiplication, addition and subtraction
# Will eventually be implemented in Oscar

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