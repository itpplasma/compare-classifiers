using Plots

include("standard_map.jl")

# Number of iterations and initial conditions
tmax = 1000
norb = 100

# Arrays to hold the results
q = Array{Float64, 2}(undef, norb, tmax+1)
p = Array{Float64, 2}(undef, norb, tmax+1)

# Initialize q at to cut all orbits and p on a regular grid
q[:,1] .= pi
p[:,1] = range(0, 2*pi, length=norb)

for k in 1:tmax
    for i in 1:norb
        q[i,k+1], p[i,k+1] = standard_map(q[i,k], p[i,k])
    end
end

# Plot the results
plot(q, p, seriestype=:scatter, markersize=0.1, legend=false)
