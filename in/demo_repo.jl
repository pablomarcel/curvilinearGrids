using CurvilinearGrids

"""Create a spherical grid as a function of (r,θ,ϕ)"""
function sphere_grid(nr, ntheta, nphi)
  r0, r1 = (1, 3) # min/max radius
  (θ0, θ1) = deg2rad.((35, 180 - 35)) # min/max polar angle
  (ϕ0, ϕ1) = deg2rad.((45, 360 - 45)) # min/max azimuthal angle

  # Linear spacing in each dimension
  # Sometimes (ξ, η, ζ) is used instead of (i, j, k), depending on preference
  r(ξ) = r0 + (r1 - r0) * ((ξ - 1) / (nr - 1))
  θ(η) = θ0 + (θ1 - θ0) * ((η - 1) / (ntheta - 1))
  ϕ(ζ) = ϕ0 + (ϕ1 - ϕ0) * ((ζ - 1) / (nphi - 1))

  x = zeros(nr, ntheta, nphi)
  y = zeros(nr, ntheta, nphi)
  z = zeros(nr, ntheta, nphi)
  # simple spherical to cartesian mapping
  for idx in CartesianIndices(x)
    i,j,k = idx.I
    x[idx] = r(i) * sin(θ(j)) * cos(ϕ(k))
    y[idx] = r(i) * sin(θ(j)) * sin(ϕ(k))
    z[idx] = r(i) * cos(θ(j))
  end

  return (x, y, z)
end

ni, nj, nk = (5, 9, 11) # number of nodes/vertices in each dimension
nhalo = 4 # halo cells needed for stencils (can be set to 0)

# Obtain the x, y, and z coordinate functions
x, y, z = sphere_grid(ni, nj, nk)

# Create the mesh
scheme = :meg6_symmetric # the symmetric scheme is more robust but more costly than :meg6
mesh = CurvilinearGrid3D(x, y, z, scheme)