function unwrap(x)
    y = Array{Float64}(undef, size(x))
    y[1] = x[1]
    for i in 1:size(x,1)-1
        if x[i+1] < x[i]
            y[i+1] = x[i+1] + 2*pi
        else
            y[i+1] = y[i] + x[i+1] - x[i]
        end
    end
    return y
end

function count_recurrences(z, order=1)
    nstep = size(z, 2)

    q = unwrap(z[1,:])
    p = z[2,:]

    kpoi_interval = findall(qi -> q[1] <= mod(qi,2*pi) < q[2], q)
    npoi_interval = length(kpoi_interval)

    N = Array{Int32}(undef, npoi_interval)
    too_short_time = Array{Bool}(undef, npoi_interval)
    for i in 1:npoi_interval
        k = kpoi_interval[i]
        too_short_time[i] = 1
        for k2 in k:nstep-1
            counter = 0
            if mod(q[k2],2*pi) < q[1] && mod(q[k2+1],2*pi) >= q[1]
                counter+=1
                if counter == order
                    N[i] = k2-k
                    too_short_time[i] = 0
                    break
                end
            end
        end
    end

    return N, too_short_time
end
