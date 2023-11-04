include("standard_map.jl")

tmax = 128

q = Array{Float64}(undef, tmax+1)
p = Array{Float64}(undef, tmax+1)

q[1] = pi
p[1] = 3.6

for k in 1:tmax
    q[k+1], p[k+1] = standard_map(q[k], p[k])
end

# Now take points in some interval and sort them
kpoi_interval = findall(q -> pi < q < 1.2*pi, q)
q_interval = q[kpoi_interval]
p_interval = p[kpoi_interval]

# Get the sorting order according to q
sort_order = sortperm(q_interval)

# Sort q_interval and p_interval using the sorting order
q_interval = q_interval[sort_order]
p_interval = p_interval[sort_order]

z2 = standard_map.(q_interval, p_interval)
z3 = standard_map.(first.(z2), last.(z2))
z4 = standard_map.(first.(z3), last.(z3))
z5 = standard_map.(first.(z4), last.(z4))
z6 = standard_map.(first.(z5), last.(z5))
z7 = standard_map.(first.(z6), last.(z6))
z8 = standard_map.(first.(z7), last.(z7))


# Plot the results
plot(q, p, seriestype=:scatter, markersize=0.5, legend=false)
plot!(q_interval, p_interval, linewidth=5, linecolor=:red, linealpha=0.5)
plot!(z2, linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z3, linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z4, linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z5, linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z6, linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z7[1:4], linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z7[5:end], linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z8, linewidth=5, linecolor=:blue, linealpha=0.5)

xlims!(0, 2*pi)
ylims!(0, 2*pi)
