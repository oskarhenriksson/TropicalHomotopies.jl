module TropicalHomotopies

    using Oscar
    import HomotopyContinuation as HC

    include("helper_functions.jl")
    include("parametric_systems.jl")
    include("tropical_geometry.jl")
    include("tropical_homotopies.jl")
    include("vertical_systems.jl")
    include("numerical_solving.jl")
    include("polyhedral_homotopies.jl")

end

