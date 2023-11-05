#%%
using Plots

include("standard_map.jl")
include("classifier.jl")

tmax = 128
K = 1.0

q = Array{Float64}(undef, tmax+1)
p = Array{Float64}(undef, tmax+1)

q[1] = 0.5*pi
p[1] = 1.2

for k in 1:tmax
    q[k+1], p[k+1] = standard_map(q[k], p[k], K)
end

krec = find_recurrences(q)
println("Classified as ideal: ", classify(krec))

#%%
# Plot the results
plot(q, p, seriestype=:scatter, markersize=0.5, legend=false)
xlims!(0, 2*pi)
ylims!(0, 2*pi)
xlabel!("q")
ylabel!("p")
savefig("chaos.png")

tmax2 = 1024

q2 = Array{Float64}(undef, tmax2+1)
p2 = Array{Float64}(undef, tmax2+1)

q2[1] = q[1]
p2[1] = p[2]

for k in 1:tmax2
    q2[k+1], p2[k+1] = standard_map(q2[k], p2[k], K)
end

krec2 = find_recurrences(q2)
println("Classified as ideal: ", classify(krec))

# Plot the results
plot(q2, p2, seriestype=:scatter, markersize=0.5, legend=false)
xlims!(0, 2*pi)
ylims!(0, 2*pi)
xlabel!("q")
ylabel!("p")
savefig("chaos2.png")
