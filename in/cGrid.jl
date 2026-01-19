# in/cGrid.jl
#
# Run:
#   julia> include("in/cGrid.jl")
#
# Output:
#   out/cGrid.vts   (open in ParaView)

module DemoCurvilinearCGrid

using CurvilinearGrids
using WriteVTK

"""
2D C-grid around a cylinder in the xy-plane.

- i direction: wraps from theta_start to theta_end (NOT periodic), leaving a cut line
- j direction: radial from R to Rfar

Returns x,y as (ni,nj) arrays.
"""
function cylinder_cgrid(ni::Int, nj::Int;
                        R::Float64 = 1.0,
                        Rfar::Float64 = 6.0,
                        radial_stretch::Float64 = 2.8,
                        # C-grid cut: leave a wedge around theta_cut_center (usually downstream)
                        theta_cut_center_deg::Float64 = 180.0,
                        theta_cut_width_deg::Float64 = 30.0)

    θc = deg2rad(theta_cut_center_deg)
    Δ  = deg2rad(theta_cut_width_deg)

    # theta range excludes the cut wedge: [θc+Δ/2, θc+2π-Δ/2]
    θ0 = θc + 0.5 * Δ
    θ1 = θc + 2π - 0.5 * Δ

    # Radial clustering near the wall (tanh)
    function radial_map(η::Float64)
        a = radial_stretch
        if a <= 1.0
            return η
        end
        return tanh(a * η) / tanh(a)
    end

    x = zeros(Float64, ni, nj)
    y = zeros(Float64, ni, nj)

    @inbounds for j in 1:nj
        η = (j - 1) / (nj - 1)
        s = radial_map(η)
        r = R + (Rfar - R) * s

        for i in 1:ni
            ξ = (i - 1) / (ni - 1)
            θ = θ0 + (θ1 - θ0) * ξ
            x[i, j] = r * cos(θ)
            y[i, j] = r * sin(θ)
        end
    end

    return x, y
end

function main()
    # Similar density as your O-grid
    ni, nj = (129, 65)

    scheme_primary  = :meg6
    scheme_fallback = :meg6_symmetric

    x, y = cylinder_cgrid(ni, nj;
                          R=1.0,
                          Rfar=6.0,
                          radial_stretch=2.8,
                          theta_cut_center_deg=180.0,   # downstream
                          theta_cut_width_deg=35.0)     # size of wake cut

    println("Raw grid sanity:")
    println("  x[1,1], y[1,1]     = ", (x[1,1], y[1,1]), "  (near one cut edge at wall)")
    println("  x[end,1], y[end,1] = ", (x[end,1], y[end,1]), "  (near other cut edge at wall)")
    println("  x[1,end], y[1,end] = ", (x[1,end], y[1,end]), "  (farfield near cut edge)")

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

    println("OK: built C-grid mesh = $(typeof(mesh)) using scheme = $(used_scheme)")

    outdir = normpath(joinpath(@__DIR__, "..", "out"))
    mkpath(outdir)

    vtkbase = joinpath(outdir, "cGrid")
    vtk = vtk_grid(vtkbase, x, y)   # 2D structured grid export
    vtk_save(vtk)

    println("OK: wrote VTK structured grid: $(vtkbase).vts")
    println("    ParaView: Surface With Edges / Wireframe")
    println("    Tip: the C-grid 'cut' is the non-periodic i=1 and i=end boundaries.")

    return nothing
end

end # module

DemoCurvilinearCGrid.main()
nothing
