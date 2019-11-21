using FrameFun: ProductPlatform
export NdCDFPlatform
"""
    NdCDFPlatform(moments::NTuple{N,Int), types=ntuple(k=>CDFPlatform, Val(N))) where N

Return a ProductPlatform that has B-spline platforms of vanishing moments `moments` as its elements.
"""
NdCDFPlatform(moments::NTuple{N,NTuple{2,Int}}) where N =
    ProductPlatform(map(m->CDFPlatform(m...), moments)...)
NdCDFPlatform(d::Int, P::Int, Q::Int) =
    NdCDFPlatform(ntuple(k->(P,Q),Val(d)))

export NdDaubechiesPlatform
"""
    NdDaubechiesPlatform(moments::NTuple{N,Int), types=ntuple(k=>CDFPlatform, Val(N))) where N

Return a ProductPlatform that has Daubechies platforms of vanishing moments `moments` as its elements.
"""
NdDaubechiesPlatform(moments::NTuple{N,Int}) where N =
    ProductPlatform(map(DaubechiesPlatform, moments)...)
NdDaubechiesPlatform(d::Int, P::Int) =
    NdDaubechiesPlatform(ntuple(k->P,Val(d)))
