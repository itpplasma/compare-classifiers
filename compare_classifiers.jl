include("simple_julia/classifier.jl")
include("simple_julia/util.jl")
include("simple_julia/standard_map.jl")

using Plots
using StatsBase
using Statistics


# we will estimate the proportion of trajectories that leave p<5


K = 0.9

# Number of iterations and initial conditions
tmax = 1024
norb = 512

### Uniform Sampling

# Arrays to hold the results
q = Array{Float64, 2}(undef, norb, tmax+1)
p = Array{Float64, 2}(undef, norb, tmax+1)

# Initialize q at to cut all orbits and p on a regular grid
q[:,1] = rand(norb)*2*pi
p[:,1] = rand(norb)*2*pi

leaves_region = fill(false, norb)
for k in 1:tmax
    for i in 1:norb
        q[i,k+1], p[i,k+1] = standard_map(q[i,k], p[i,k], K)
        if i > 1
            leaves_region[i] = leaves_region[i] || p[i, k+1] > 5
        end
    end
end





cols = [:red, :blue]
plot(q[leaves_region,:], p[leaves_region,:], seriestype=:scatter, color =:red, markersize=0.5, markerstrokewidth=0, legend=false)

plot!(q[.!leaves_region,:], p[.!leaves_region,:], seriestype=:scatter, color=:blue, markersize=0.5, markerstrokewidth=0, legend=false)

xlabel!("q")
ylabel!("p")

title!("Confined particles (blue)")
savefig("confinement_map.png")

plot(1:norb, cumsum(leaves_region) ./ (1:norb), series_type=:line)
xlabel!("n_particles")
ylabel!("Proportion leaving")
title!("LLN for proportion leaving (Unweighted)")
savefig("classifer_orbits/unweighted_confinement.png")




### Importance Sampling by past conclusions

# define weight function
function sampling_weight(prob, mean)
    #return abs(prob-mean)
    if prob > mean
        return 1 - (prob - mean)/(1-mean)
    end
    return prob / mean
end

function sample_square(idx, l, n)
    y_idx = idx % n
    x_idx = idx ÷ n

    x_val = rand()*l + x_idx*l
    y_val = rand()*l + y_idx*l
    return x_val, y_val
end


function build_cond_ev(p, q, leaves_region, max_idx, l, n)
    apriori_mean = sum(leaves_region[1:max_idx])/max_idx
    interp_matrix = fill(0.0, n, n)
    
    counts = fill(0, n, n)
    for i in 1:max_idx
        for k in 1:tmax
            x_idx = min(Int(ceil(q[i, k] / l)), n)
            y_idx = min(Int(ceil(p[i, k] / l)), n) # fixes an odd bug for now.
            
            if leaves_region[i]
                interp_matrix[y_idx, x_idx] += 1
            end
            counts[y_idx, x_idx] += 1
        end
    end 
    
    # naive probability if no better information
    interp_matrix = interp_matrix ./ counts
    interp_matrix[isnan.(interp_matrix)] .= apriori_mean
    weight_mat = map(x -> sampling_weight(x, apriori_mean), interp_matrix) 
    return interp_matrix , weight_mat
end


# Arrays to hold the results
q = Array{Float64, 2}(undef, norb, tmax+1)
p = Array{Float64, 2}(undef, norb, tmax+1)
warm_start = 50
sample_size = 100

# Initialize q at to cut all orbits and p on a regular grid
q[1:warm_start,1] = rand(warm_start)*2*pi
p[1:warm_start,1] = rand(warm_start)*2*pi

leaves_region = fill(false, norb) 
outcome_weight = fill(0.0, norb) # observations reweighted from importance sampling

norb_done = 0
# create initial data
for i in 1:warm_start
    for k in 1:tmax
        q[i,k+1], p[i,k+1] = standard_map(q[i,k], p[i,k], K)
        if i > 1
            leaves_region[i] = leaves_region[i] || p[i, k+1] > 5
        end
    end
end
outcome_weight[1:warm_start] .= (1/(4*pi^2)) # both are equal weighted (changed to 1 later)
global norb_done += warm_start

plot(q[1:norb_done,:], p[1:norb_done, :], seriestype=:scatter, zcolor=leaves_region[1:warm_start], color=:roma, markersize=0.5, markerstrokewidth=0, legend=false)
xlabel!("q")
ylabel!("p")

title!("$norb_done orbits complete")
savefig("classifer_orbits/confinement_map_orbits$norb_done.png") 


while norb_done < norb
    # calculate sampling probabilities

    # record known outcomes
    n = Int(floor(sqrt(norb)))
    l = 2*pi / n
    
    interp_matrix, weight_mat = build_cond_ev(p, q, leaves_region, norb_done, l, n) 
    
    #println("calc weights")
    #weight_mat = map(x -> sampling_weight(x, apriori_mean), interp_matrix)   
    
    # record current weights as a heatmap
    heatmap(weight_mat)
    title!("Particle outcome confidence")
    savefig("classifer_orbits/sampling_weights$norb_done.png") 

 


    # proceed with potentially smaller batches
    start_idx = norb_done+1 
    end_idx = min(norb_done+sample_size, norb)
    
    # uniform sampling
    #q[start_idx:end_idx,1] = rand(end_idx - start_idx + 1)*2*pi
    #p[start_idx:end_idx,1] = rand(end_idx - start_idx + 1)*2*pi

    # importance sampling
    # use sampling weights to get indices
    weight_vec = vec(weight_mat)
    total_weight = sum(weight_vec)
    next_indices = StatsBase.wsample(1:n^2, vec(weight_mat), end_idx - start_idx + 1)
    observed_weights = weight_vec[next_indices] ./ (total_weight * l^2)
    outcome_weight[start_idx:end_idx] = observed_weights

    next_start_pts = map(idx -> sample_square(idx, l, n), next_indices)
    q[start_idx:end_idx,1]= getfield.(next_start_pts, 1)
    p[start_idx:end_idx,1] = getfield.(next_start_pts, 2)


    for i in start_idx:end_idx
        for k in 1:tmax
            q[i,k+1], p[i,k+1] = standard_map(q[i,k], p[i,k], K)
            if i > 1
                leaves_region[i] = leaves_region[i] || p[i, k+1] > 5
            end
        end
    end

    global norb_done = end_idx
    plot(q[1:norb_done,:], p[1:norb_done, :], seriestype=:scatter, zcolor=leaves_region[1:norb_done], color=:roma, markersize=0.5, markerstrokewidth=0, legend=false)
    xlabel!("q")
    ylabel!("p")
    
    title!("$norb_done orbits complete")
    savefig("classifer_orbits/confinement_map_orbits$norb_done.png") 


end

outcome_weight .= (1/(4*pi^2)) ./ outcome_weight
println(outcome_weight'leaves_region)

plot(1:norb, cumsum(outcome_weight .* leaves_region) ./ (1:norb), series_type=:line)
xlabel!("n_particles")
ylabel!("Proportion leaving")
title!("LLN for proportion leaving (Weighted)")
savefig("classifer_orbits/weighted_confinement.png")

n = Int(floor(sqrt(norb)))
l = 2*pi / n

cond_ev = fill(0.0, norb)
for i=1:norb
    interp_matrix, _ = build_cond_ev(p, q, leaves_region, i, l, n)
    evs = vec(interp_matrix)
    cond_ev[i] = sum(evs) / length(evs)
end
plot(1:norb, cond_ev, series_type=:line)
xlabel!("n_particles")
ylabel!("Proportion leaving")
title!("LLN for proportion leaving (Conditional EVs)")
savefig("classifer_orbits/cond_ev_confinement.png")

println("done")