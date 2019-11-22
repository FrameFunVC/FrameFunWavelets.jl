module DyadicPeriodicEquispacedGrids

using GridArrays, DomainSets

using DomainSets: endpoints

import Base: promote_eltype, ==, size
import GridArrays: name, support, isperiodic, length, similargrid, hasextension,
    extend, resize, PeriodicEquispacedGrid, _extension_size, rescale

export isdyadic
isdyadic(n::Int) = n≈(2.0^round(log2(n)))
export DyadicPeriodicEquispacedGrid
"""
    struct DyadicPeriodicEquispacedGrid{T} <: AbstractEquispacedGrid{T}

A dyadic periodic equispaced grid is an equispaced grid that omits the right
endpoint and has length `n = 2^l`.
It has step `(b-a)/n`.
"""
struct DyadicPeriodicEquispacedGrid{T} <: AbstractEquispacedGrid{T}
    range   :: LinRange{T}
    l   ::  Int
    a   ::  T
    b   ::  T

    DyadicPeriodicEquispacedGrid{T}(n::Int, a, b) where {T} =
        (@assert isdyadic(n);new(LinRange(T(a),T(b),n+1)[1:end-1],Int(log2(n)),T(a),T(b)))
end


name(g::DyadicPeriodicEquispacedGrid) = "Dyadic periodic equispaced grid"
support(grid::DyadicPeriodicEquispacedGrid) = Interval(grid.a, grid.b)
isperiodic(::DyadicPeriodicEquispacedGrid) = true

export dyadic_length
"""
    dyadic_length(grid::DyadicPeriodicEquispacedGrid)

The dyadic length of the grid
"""
dyadic_length(g::DyadicPeriodicEquispacedGrid) = g.l
length(g::DyadicPeriodicEquispacedGrid) = length(g.range)

DyadicPeriodicEquispacedGrid(n::Int, d::AbstractInterval) =
    DyadicPeriodicEquispacedGrid(n, endpoints(d)...)
similargrid(grid::DyadicPeriodicEquispacedGrid, ::Type{T}, n::Int) where {T} =
    DyadicPeriodicEquispacedGrid{T}(n, map(T, endpoints(support(grid)))...)
rescale(grid::DyadicPeriodicEquispacedGrid, a, b) =
    DyadicPeriodicEquispacedGrid{promote_type(typeof(a/2),typeof(b/2),eltype(grid))}(length(grid), a, b)
DyadicPeriodicEquispacedGrid(n::Int, a, b) =
    DyadicPeriodicEquispacedGrid{promote_type(typeof(a/2),typeof(b/2))}(n, a, b)
mapped_grid(grid::DyadicPeriodicEquispacedGrid, map::AffineMap) =
    DyadicPeriodicEquispacedGrid(length(grid), endpoints(map*support(grid))...)

"""
    isdyadic(grid::AbstractIntervalGrid)

Does the grid have a dyadic length?
"""
isdyadic(::DyadicPeriodicEquispacedGrid) = true
isdyadic(grid::PeriodicEquispacedGrid) = isdyadic(length(grid))
isdyadic(::AbstractIntervalGrid) = false

PeriodicEquispacedGrid(g::DyadicPeriodicEquispacedGrid{T}) where {T} = PeriodicEquispacedGrid{T}(length(g), g.a, g.b)
DyadicPeriodicEquispacedGrid(g::PeriodicEquispacedGrid{T}) where {T} =
    (@assert isdyadic(g);DyadicPeriodicEquispacedGrid(length(g),g.a,g.b))

_extension_size(::DyadicPeriodicEquispacedGrid, n::Int, factor::Int) = (@assert isdyadic(factor);factor*n)
hasextension(::DyadicPeriodicEquispacedGrid) = true
extend(grid::DyadicPeriodicEquispacedGrid, factor::Int) =
    resize(grid, _extension_size(grid, length(grid), factor))

similargrid(g::DyadicPeriodicEquispacedGrid, a, b, T) = DyadicPeriodicEquispacedGrid{T}(g.l, a, b)
resize(g::DyadicPeriodicEquispacedGrid, n::Int) = DyadicPeriodicEquispacedGrid(n, g.a, g.b)


# We need this basic definition, otherwise equality does not seem to hold when T is BigFloat...
==(g1::DyadicPeriodicEquispacedGrid, g2::DyadicPeriodicEquispacedGrid) =
    (g1.l == g2.l) && (g1.a == g2.a) && (g1.b==g2.b)
==(g1::DyadicPeriodicEquispacedGrid, g2::PeriodicEquispacedGrid) =
    (length(g1) == length(g2)) && (g1.a == g2.a) && (g1.b==g2.b)
end
