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

# https://math.stackexchange.com/questions/2275439/check-if-point-on-circle-is-between-two-other-points-on-circle
function is_x3_between_x1_and_x2(x1, x2, x3)
    x1_x2 = mod(x2 - x1 + 2π, 2π)
    x1_x3 = mod(x3 - x1 + 2π, 2π)

    return (x1_x2 <= π) ⊻ (x1_x3 > x1_x2)
end


function count_recurrences(z, order=1)
    nstep = size(z, 2)

    q = unwrap(z[1,:])
    p = z[2,:]

    kpoi_interval = findall(qi -> is_x3_between_x1_and_x2(q[1], q[2], qi), q)
    npoi_interval = length(kpoi_interval)

    N = Array{Int32}(undef, npoi_interval)
    too_short_time = Array{Bool}(undef, npoi_interval)
    for i in 1:npoi_interval
        k = kpoi_interval[i]
        too_short_time[i] = 1
        for k2 in k:nstep-1
            counter = 0
            if is_x3_between_x1_and_x2(q[1], q[2], q[k2+1])
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
