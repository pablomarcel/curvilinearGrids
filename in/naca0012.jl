# in/naca0012.jl
#
# Run:
#   julia> include("in/naca0012.jl")
#
# Output:
#   out/naca0012.vts

module DemoNACA0012

using CurvilinearGrids
using WriteVTK

# ------------------------------
# NACA0012 geometry (closed trailing edge)
# ------------------------------
# Thickness distribution (closed TE uses -0.1036)
@inline function naca_00xx_thickness(x::Float64, t::Float64)
    # Guard for sqrt near x=0
    xs = max(x, 0.0)
    return 5.0 * t * (0.2969 * sqrt(xs) - 0.1260 * xs - 0.3516 * xs^2 +
                      0.2843 * xs^3 - 0.1036 * xs^4)
end

"""
Create a closed NACA0012 surface as a single periodic loop (TE duplicated at end).

Returns:
  xs, ys :: Vector{Float64} length ni
where:
  ni = 2*n_side - 1
and point 1 and point end are both TE (x=1,y=0).
"""
function naca0012_loop(n_side::Int; chord::Float64 = 1.0, t::Float64 = 0.12)
    @assert n_side >= 5

    # Cosine spacing x in [0,1] (clusters LE and TE)
    s = range(0.0, 1.0; length = n_side)
    x = 0.5 .* (1 .- cos.(π .* s))  # 0..1

    # Upper: TE->LE inclusive (reverse x)
    xU = reverse(x)
    yU = [naca_00xx_thickness(xi, t) for xi in xU]

    # Lower: LE->TE excluding LE and TE endpoints (avoid duplicates)
    xL = x[2:end-1]
    yL = [-naca_00xx_thickness(xi, t) for xi in xL]

    # Assemble loop and duplicate TE at end
    xs = vcat(chord .* xU, chord .* xL, chord * 1.0)
    ys = vcat(chord .* yU, chord .* yL, 0.0)

    return xs, ys
end

# ------------------------------
# O-grid generator around the airfoil
# ------------------------------
"""
Body-fitted O-grid around NACA0012:
- i: around airfoil (periodic with duplicate TE endpoint)
- j: outward to circular farfield

Returns x,y as (ni, nj) matrices.
"""
function naca0012_ogrid(n_side::Int, nj::Int;
                        chord::Float64 = 1.0,
                        t::Float64 = 0.12,
                        Rfar::Float64 = 20.0,
                        radial_stretch::Float64 = 3.0,
                        smooth_iters::Int = 150,
                        smooth_omega::Float64 = 1.2)

    xs, ys = naca0012_loop(n_side; chord=chord, t=t)
    ni = length(xs)

    # Center the farfield circle at mid-chord for a clean O-grid
    xc, yc = 0.5 * chord, 0.0

    # Angle of each surface point around (xc,yc)
    θ = [atan(ys[i] - yc, xs[i] - xc) for i in 1:ni]

    # Outer boundary (circle)
    x_out = [xc + Rfar * cos(θi) for θi in θ]
    y_out = [yc + Rfar * sin(θi) for θi in θ]

    # Radial stretching (clusters near wall if >1)
    radial_map(η::Float64) = radial_stretch <= 1.0 ? η : tanh(radial_stretch * η) / tanh(radial_stretch)

    x = zeros(Float64, ni, nj)
    y = zeros(Float64, ni, nj)

    # Boundary j=1 (airfoil) and j=nj (farfield)
    @inbounds for i in 1:ni
        x[i, 1]  = xs[i]
        y[i, 1]  = ys[i]
        x[i, nj] = x_out[i]
        y[i, nj] = y_out[i]
    end

    # Interior via algebraic interpolation
    @inbounds for j in 2:nj-1
        η = (j - 1) / (nj - 1)
        sη = radial_map(η)
        for i in 1:ni
            x[i, j] = (1 - sη) * xs[i]    + sη * x_out[i]
            y[i, j] = (1 - sη) * ys[i]    + sη * y_out[i]
        end
    end

    # Optional Laplace smoothing (keeps boundaries fixed)
    # periodic in i via duplicated TE point (i=1 and i=ni)
    if smooth_iters > 0
        @inbounds for it in 1:smooth_iters
            for j in 2:nj-1
                for i in 2:ni-1
                    xn = 0.25 * (x[i-1,j] + x[i+1,j] + x[i,j-1] + x[i,j+1])
                    yn = 0.25 * (y[i-1,j] + y[i+1,j] + y[i,j-1] + y[i,j+1])
                    x[i,j] = (1 - smooth_omega) * x[i,j] + smooth_omega * xn
                    y[i,j] = (1 - smooth_omega) * y[i,j] + smooth_omega * yn
                end
            end
            # enforce seam duplicate explicitly
            for j in 1:nj
                x[ni,j] = x[1,j]
                y[ni,j] = y[1,j]
            end
        end
    end

    return x, y
end

function main()
    # --- resolution
    n_side = 129           # TE->LE points (upper). Total ni = 2*n_side - 1
    nj     = 129           # radial points

    # --- grid params
    chord = 1.0
    Rfar  = 25.0
    radial_stretch = 3.2

    # smoothing: increase iters for prettier grid
    smooth_iters = 200
    smooth_omega = 1.15

    x, y = naca0012_ogrid(n_side, nj;
                          chord=chord,
                          t=0.12,
                          Rfar=Rfar,
                          radial_stretch=radial_stretch,
                          smooth_iters=smooth_iters,
                          smooth_omega=smooth_omega)

    # quick sanity
    ni = size(x,1)
    println("NACA0012 O-grid:")
    println("  size(x) = ", size(x))
    println("  TE seam check: (x[1,1],y[1,1]) = ", (x[1,1],y[1,1]), "  |  (x[end,1],y[end,1]) = ", (x[ni,1],y[ni,1]))
    println("  farfield sample: (x[1,end],y[1,end]) = ", (x[1,end],y[1,end]))

    # build metrics object (optional but good to validate)
    scheme_primary  = :meg6
    scheme_fallback = :meg6_symmetric
    try
        _mesh = CurvilinearGrid2D(x, y, scheme_primary)
        println("  OK: CurvilinearGrid2D metrics built (scheme=$(scheme_primary))")
    catch err
        println("  WARNING: metrics failed with scheme=$(scheme_primary): ", err)
        println("           retrying scheme=$(scheme_fallback)")
        _mesh = CurvilinearGrid2D(x, y, scheme_fallback)
        println("  OK: CurvilinearGrid2D metrics built (scheme=$(scheme_fallback))")
    end

    # write VTK
    outdir = normpath(joinpath(@__DIR__, "..", "out"))
    mkpath(outdir)

    vtkbase = joinpath(outdir, "naca0012")
    vtk = vtk_grid(vtkbase, x, y)   # 2D structured grid export
    vtk_save(vtk)

    println("  OK: wrote ", vtkbase, ".vts")
    println("  ParaView: Surface With Edges / Wireframe; view from +Z for 2D look")

    return nothing
end

end # module

DemoNACA0012.main()
nothing
