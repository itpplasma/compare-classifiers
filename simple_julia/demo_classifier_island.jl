#%%
using Plots

include("standard_map.jl")

tmax = 128
K = 0.6

q = Array{Float64}(undef, tmax+1)
p = Array{Float64}(undef, tmax+1)

q[1] = 0.5*pi
p[1] = 0.4

for k in 1:tmax
    q[k+1], p[k+1] = standard_map(q[k], p[k], K)
end

# Now take points in some interval and sort them
kpoi_interval = findall(qi -> q[1] <= qi <= q[2], q)
q_interval = q[kpoi_interval]
p_interval = p[kpoi_interval]

z = [q';p']

krec = find_recurrences(q)
println("Classified as ideal: ", classify(krec))

#%%

# Get the sorting order according to q
sort_order = sortperm(q_interval)

# Sort q_interval and p_interval using the sorting order
q_interval = q_interval[sort_order]
p_interval = p_interval[sort_order]

z2 = standard_map.(q_interval, p_interval, K)
z3 = standard_map.(first.(z2), last.(z2), K)
z4 = standard_map.(first.(z3), last.(z3), K)
z5 = standard_map.(first.(z4), last.(z4), K)
z6 = standard_map.(first.(z5), last.(z5), K)
z7 = standard_map.(first.(z6), last.(z6), K)
z8 = standard_map.(first.(z7), last.(z7), K)

# Plot the results
plot(q, p, seriestype=:scatter, markersize=0.5, legend=false)
xlims!(0, 2*pi)
ylims!(0, 2*pi)
savefig("island01.png")
plot!(q_interval[p_interval.<pi], p_interval[p_interval.<pi], linewidth=5, linecolor=:red, linealpha=0.5)
plot!(q_interval[p_interval.>=pi], p_interval[p_interval.>=pi], linewidth=5, linecolor=:red, linealpha=0.5)
savefig("island02.png")
plot!(z2[last.(z2).<pi], linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z2[last.(z2).>=pi], linewidth=5, linecolor=:green, linealpha=0.5)
savefig("island03.png")
plot!(z3[last.(z3).<pi .&& first.(z3).<pi], linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z3[last.(z3).<pi .&& first.(z3).>=pi], linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z3[last.(z3).>=pi .&& first.(z3).<pi], linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z3[last.(z3).>=pi .&& first.(z3).>=pi], linewidth=5, linecolor=:green, linealpha=0.5)
savefig("island04.png")
plot!(z4[last.(z4).<pi .&& first.(z4).<pi], linewidth=5, linecolor=:blue, linealpha=0.5)
plot!(z4[last.(z4).<pi .&& first.(z4).>=pi], linewidth=5, linecolor=:blue, linealpha=0.5)
plot!(z4[last.(z4).>=pi .&& first.(z4).<pi], linewidth=5, linecolor=:blue, linealpha=0.5)
plot!(z4[last.(z4).>=pi .&& first.(z4).>=pi], linewidth=5, linecolor=:blue, linealpha=0.5)
savefig("island05.png")
plot!(z5[last.(z5).<pi .&& first.(z5).<4], linewidth=5, linecolor=:yellow, linealpha=0.5)
plot!(z5[last.(z5).<pi .&& first.(z5).>=4], linewidth=5, linecolor=:yellow, linealpha=0.5)
plot!(z5[last.(z5).>=pi .&& first.(z5).<4], linewidth=5, linecolor=:yellow, linealpha=0.5)
plot!(z5[last.(z5).>=pi .&& first.(z5).>=4], linewidth=5, linecolor=:yellow, linealpha=0.5)
savefig("island06.png")
