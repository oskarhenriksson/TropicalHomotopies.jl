function homotopy_from_tropical_data(G,w)
    T = tropical_semiring(min)
    Ktx = parent(first(G))
    x = gens(Ktx)
    Kt = coefficient_ring(Ktx)
    t = gens(Kt)[1]
    nu = tropical_semiring_map(Kt,t) 
    @req length(G) == ngens(parent(first(G))) "Input system needs to be square"
    r = Int(lcm(denominator.(w))) # for clearing denominators in the t exponents
    G_reparametrized = hom(Ktx,Ktx,c->c(t^r),gens(Ktx)).(G) 
    Htilde = hom(Ktx,Ktx,t.^(Int.(r*w)).*x).(G_reparametrized)
    valuation_shift = -Int.(evaluate.(tropical_polynomial.(G_reparametrized,Ref(nu)),Ref(T.(r*w))))
    H = (t.^valuation_shift).*Htilde
    return H
end