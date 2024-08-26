# We need this to compute initials for rational points
function Base.:(*)(x::TropicalSemiringElem, y::QQFieldElem)
    iszero(x) && return x # if x is zero, return it
    return parent(x)(data(x) + y) # otherwise, return their sum
end