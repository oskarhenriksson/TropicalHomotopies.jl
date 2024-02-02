using TropicalHomotopies

@testset "groebner_basis_workaround(I::MPolyIdeal,nu::TropicalSemiringMap,w::Vector{QQFieldElem})" begin

    R, x = polynomial_ring(QQ,"x"=>1:3)
    G = groebner_basis_workaround(ideal([x[1]+x[2]+x[3]+1,x[1]-2*x[2]+x[3]+1]),tropical_semiring_map(QQ),QQ.([1,0,1]))
    @test issetequal(G,[x[2], x[1] + x[3] + 1])

    R, x = polynomial_ring(QQ,"x"=>1:3)
    I = ideal([x[1]^2-x[2]*x[3],
    x[3]^2-2*x[2]^2,
    x[1]^2*x[3]^2-2*x[2]^3*x[3]])
    G = groebner_basis_workaround(I,tropical_semiring_map(QQ),QQ.([1,0,1]))
    @test issetequal(G,[x[1]^2*x[2] - 1//2*x[3]^3, -2*x[1]^4 + x[3]^4])
end
