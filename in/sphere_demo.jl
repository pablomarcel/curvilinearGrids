using CurvilinearGrids

function sphere_grid(nr, ntheta, nphi)
    r0, r1 = (1, 3)
    (θ0, θ1) = deg2rad.((35, 180 - 35))
    (ϕ0, ϕ1) = deg2rad.((45, 360 - 45))

    r(ξ) = r0 + (r1 - r0) * ((ξ - 1) / (nr - 1))
    θ(η) = θ0 + (θ1 - θ0) * ((η - 1) / (ntheta - 1))
    ϕ(ζ) = ϕ0 + (ϕ1 - ϕ0) * ((ζ - 1) / (nphi - 1))

    x = zeros(nr, ntheta, nphi)
    y = zeros(nr, ntheta, nphi)
    z = zeros(nr, ntheta, nphi)

    for idx in CartesianIndices(x)
        i, j, k = idx.I
        x[idx] = r(i) * sin(θ(j)) * cos(ϕ(k))
        y[idx] = r(i) * sin(θ(j)) * sin(ϕ(k))
        z[idx] = r(i) * cos(θ(j))
    end
    return x, y, z
end

ni, nj, nk = (5, 9, 11)
scheme = :meg6_symmetric

x, y, z = sphere_grid(ni, nj, nk)

# Constructor signature can vary by version; check methods if unsure:
# methods(CurvilinearGrid3D)

mesh = CurvilinearGrid3D(x, y, z, scheme)

println(typeof(mesh))
println(coord(mesh, (1, 1, 1)))
