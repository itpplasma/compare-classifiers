using Plots

include("standard_map.jl")
include("classifier.jl")

K = 0.4

# Number of iterations and initial conditions
tmax = 1024
norb = 128

# Arrays to hold the results
q = Array{Float64, 2}(undef, norb, tmax+1)
p = Array{Float64, 2}(undef, norb, tmax+1)

# Initialize q at to cut all orbits and p on a regular grid
q[:,1] .= pi
p[:,1] = range(0, 2*pi, length=norb)

for k in 1:tmax
    for i in 1:norb
        q[i,k+1], p[i,k+1] = standard_map(q[i,k], p[i,k], K)
    end
end


k_ideal_old = Int[]
# Plot the results
for kclass = (16, 32, 64, 128, 256, 512)
    k_ideal = Int[]
    k_non_ideal = Int[]

    for i in 1:norb
        ideal = classify(find_recurrences(q[i,1:kclass]))
        ideal_reverse = classify(find_recurrences(q[i,kclass:-1:1]))
        if ideal || ideal_reverse
            push!(k_ideal, i)
        else
            push!(k_non_ideal, i)
        end
    end

    println(length(k_ideal_old) - length(k_ideal))
    k_ideal_old = k_ideal

    plot(q[k_ideal,:], p[k_ideal,:], seriestype=:scatter, color=:blue, markersize=0.5, markerstrokewidth=0, legend=false)

    plot!(q[k_non_ideal,:], p[k_non_ideal,:], seriestype=:scatter, color=:red, markersize=0.5, markerstrokewidth=0, legend=false)

    xlabel!("q")
    ylabel!("p")

    title!("$kclass footprints")
    savefig("footprints_$kclass.png")
end
