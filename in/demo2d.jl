# in/demo2d.jl
#
# Run:
#   julia> include("in/demo2d.jl")
#
# Output:
#   out/out2d.vts   (open in ParaView)

module DemoCurvilinearGrids2D

using CurvilinearGrids
using WriteVTK

"""
Build a 2D O-grid around a cylinder in the xy-plane.

ni: points around circumference
nj: points radially outward

Returns: x, y as (ni, nj) arrays.
"""
function cylinder_ogrid(ni::Int, nj::Int;
                        R::Float64 = 1.0,
                        Rfar::Float64 = 6.0,
                        radial_stretch::Float64 = 2.8)
    x = zeros(Float64, ni, nj)
    y = zeros(Float64, ni, nj)

    # Radial clustering near the wall (tanh)
    function radial_map(η::Float64)
        a = radial_stretch
        if a <= 1.0
            return η
        end
        return tanh(a * η) / tanh(a)
    end

    @inbounds for j in 1:nj
        η = (j - 1) / (nj - 1)
        s = radial_map(η)
        r = R + (Rfar - R) * s

        for i in 1:ni
            ξ = (i - 1) / (ni - 1)
            θ = 2π * ξ
            x[i, j] = r * cos(θ)
            y[i, j] = r * sin(θ)
        end
    end

    return x, y
end

function main()
    # Big enough for MEG6 metrics
    ni, nj = (129, 65)

    scheme_primary  = :meg6
    scheme_fallback = :meg6_symmetric

    x, y = cylinder_ogrid(ni, nj; R=1.0, Rfar=6.0, radial_stretch=2.8)

    # --- sanity prints from raw coordinates (not mesh)
    println("Raw grid sanity:")
    println("  x[1,1], y[1,1] = ", (x[1,1], y[1,1]), "  (should be ~ (R, 0))")
    println("  x[1,end], y[1,end] = ", (x[1,end], y[1,end]), "  (should be ~ (Rfar, 0))")
    println("  x[mid,1], y[mid,1] = ", (x[cld(ni,2),1], y[cld(ni,2),1]))

    local mesh
    local used_scheme = scheme_primary

    try
        mesh = CurvilinearGrid2D(x, y, scheme_primary)
    catch err
        println("WARNING: 2D mesh build failed with scheme = $(scheme_primary).")
        println("         Error: ", err)
        println("         Retrying with scheme = $(scheme_fallback) ...")
        used_scheme = scheme_fallback
        mesh = CurvilinearGrid2D(x, y, scheme_fallback)
    end

    println("OK: built 2D mesh = $(typeof(mesh)) using scheme = $(used_scheme)")

    # NOTE: coord(mesh, ...) indexing conventions vary; don't rely on it yet.
    # We'll trust x/y for visualization first.

    # --- write outputs under repo/out/
    outdir = normpath(joinpath(@__DIR__, "..", "out"))
    mkpath(outdir)

    vtkbase = joinpath(outdir, "out2d")

    # IMPORTANT: WriteVTK structured 2D grid takes (x, y) matrices (NOT x,y,z)
    vtk = vtk_grid(vtkbase, x, y)
    vtk_save(vtk)

    println("OK: wrote VTK structured grid: $(vtkbase).vts")
    println("    ParaView: Representation -> 'Surface With Edges' or 'Wireframe'")
    println("    If you want only lines: Filters -> Extract Edges")

    return nothing
end

end # module

DemoCurvilinearGrids2D.main()
nothing
