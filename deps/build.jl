if VERSION < v"0.7-"
    Pkg.clone("https://github.com/vincentcp/WaveletsCopy.jl.git")
    Pkg.build("WaveletsCopy")
    Pkg.clone("https://github.com/daanhb/BasisFunctions.jl.git")
    Pkg.checkout("BasisFunctions","extract")
    Pkg.build("BasisFunctions")
    Pkg.clone("https://github.com/daanhb/Domains.jl.git")
end
