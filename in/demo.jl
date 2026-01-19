# in/demo.jl
#
# Run:
#   julia> include("in/demo.jl")
#
# Output:
#   out/sphere_grid.vts   (open in ParaView)

module DemoCurvilinearGrids

using CurvilinearGrids
using WriteVTK

# ------------------------------
# Grid definition (spherical shell patch)
# ------------------------------
function sphere_grid(nr::Int, ntheta::Int, nphi::Int;
                     r0::Float64 = 1.0, r1::Float64 = 3.0,
                     θ0_deg::Float64 = 35.0, θ1_deg::Float64 = 180.0 - 35.0,
                     ϕ0_deg::Float64 = 45.0, ϕ1_deg::Float64 = 360.0 - 45.0)

    θ0, θ1 = deg2rad(θ0_deg), deg2rad(θ1_deg)
    ϕ0, ϕ1 = deg2rad(ϕ0_deg), deg2rad(ϕ1_deg)

    # Linear spacing
    r(ξ) = r0 + (r1 - r0) * ((ξ - 1) / (nr - 1))
    θ(η) = θ0 + (θ1 - θ0) * ((η - 1) / (ntheta - 1))
    ϕ(ζ) = ϕ0 + (ϕ1 - ϕ0) * ((ζ - 1) / (nphi - 1))

    x = zeros(Float64, nr, ntheta, nphi)
    y = zeros(Float64, nr, ntheta, nphi)
    z = zeros(Float64, nr, ntheta, nphi)

    @inbounds for idx in CartesianIndices(x)
        i, j, k = idx.I
        rr = r(i)
        th = θ(j)
        ph = ϕ(k)
        s = sin(th)
        x[idx] = rr * s * cos(ph)
        y[idx] = rr * s * sin(ph)
        z[idx] = rr * cos(th)
    end

    return x, y, z
end

# ------------------------------
# Main runner (avoids REPL soft-scope issues)
# ------------------------------
function main()
    # Make these reasonably large for 6th-order metrics (MEG6)
    ni, nj, nk = (21, 33, 33)

    scheme_primary  = :meg6
    scheme_fallback = :meg6_symmetric

    x, y, z = sphere_grid(ni, nj, nk)

    local mesh
    local used_scheme = scheme_primary

    try
        mesh = CurvilinearGrid3D(x, y, z, scheme_primary)
    catch err
        println("WARNING: mesh build failed with scheme = $(scheme_primary).")
        println("         Error: ", err)
        println("         Retrying with scheme = $(scheme_fallback) ...")
        used_scheme = scheme_fallback
        mesh = CurvilinearGrid3D(x, y, z, scheme_fallback)
    end

    println("OK: built mesh = $(typeof(mesh)) using scheme = $(used_scheme)")
    println("OK: example coord(mesh, (1,1,1)) = ", coord(mesh, (1, 1, 1)))

    # --- write outputs under repo/out/
    outdir = normpath(joinpath(@__DIR__, "..", "out"))
    mkpath(outdir)

    vtkbase = joinpath(outdir, "sphere_grid")
    vtk = vtk_grid(vtkbase, x, y, z)
    vtk_save(vtk)

    println("OK: wrote VTK structured grid: $(vtkbase).vts")
    println("    ParaView tip: Representation -> 'Surface With Edges' or 'Wireframe'")

    return nothing
end

end # module

DemoCurvilinearGrids.main()
nothing
