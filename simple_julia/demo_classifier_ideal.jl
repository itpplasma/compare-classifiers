using Plots

include("standard_map.jl")
include("classifier.jl")

tmax = 128
K = 0.6

q = Array{Float64}(undef, tmax+1)
p = Array{Float64}(undef, tmax+1)

q[1] = 0.5*pi
p[1] = 0.9

for k in 1:tmax
    q[k+1], p[k+1] = standard_map(q[k], p[k], K)
end

z = [q';p']

N1,too_short_time = count_recurrences(z)
println("Recurrence numbers N1:")
println(N1[.!too_short_time])

N2,too_short_time = count_recurrences(z, 2)
println("Recurrence numbers N2:")
println(N2[.!too_short_time])

# Now take points in some interval and sort them
kpoi_interval = findall(qi -> q[1] <= qi < q[2], q)
q_interval = q[kpoi_interval]
p_interval = p[kpoi_interval]

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
z8full = standard_map.(first.(z7), last.(z7), K)
z8 = standard_map.(first.(z7[1:15]), last.(z7[1:15]), K)

# Plot the results
plot(q, p, seriestype=:scatter, markersize=0.5, legend=false)
xlims!(0, 2*pi)
ylims!(0, 2*pi)
xlabel!("q")
ylabel!("p")
savefig("1.png")
plot!(q_interval, p_interval, linewidth=5, linecolor=:red, linealpha=0.5)
savefig("2.png")
plot!(z2, linewidth=5, linecolor=:green, linealpha=0.5)
savefig("3.png")
plot!(z3, linewidth=5, linecolor=:green, linealpha=0.5)
savefig("4.png")
plot!(z4, linewidth=5, linecolor=:green, linealpha=0.5)
savefig("5.png")
plot!(z5[1:8], linewidth=5, linecolor=:green, linealpha=0.5)
plot!(z5[9:end], linewidth=5, linecolor=:green, linealpha=0.5)
savefig("6.png")
plot!(z6, linewidth=5, linecolor=:green, linealpha=0.5)
savefig("7.png")
plot!(z7, linewidth=5, linecolor=:yellow, linealpha=0.5)
savefig("8.png")
plot!(z8, linewidth=5, linecolor=:orange, linealpha=0.5)
savefig("9.png")
#plot!(z7[invperm(sort_order)][N1.==5], seriestype=:scatter, color=:yellow)
savefig("earlypoints.png")
#plot!(z8full[invperm(sort_order)][N1.==6], seriestype=:scatter, color=:orange)
savefig("latepoints.png")
