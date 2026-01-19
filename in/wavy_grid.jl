# Define a simple 2D wavy grid
ni, nj = (21, 21)
x = zeros(ni, nj)
y = zeros(ni, nj)
for j in 1:nj, i in 1:ni
    x[i,j] = (i-1) + 0.2 * sin(pi * (j-1)/(nj-1))
    y[i,j] = (j-1) + 0.2 * sin(pi * (i-1)/(ni-1))
end

# Create the CurvilinearGrid object (used for metrics/physics)
mesh = CurvilinearGrid2D(x, y, :meg6)
