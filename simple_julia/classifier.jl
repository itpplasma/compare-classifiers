function count_single_recurrences(z)
    nstep = size(z, 2)
    q0 = z[1,1]
    q1 = z[1,2]

    kpoi_interval = findall(qi -> q0 <= qi < q1, z[1,:])
    npoi_interval = length(kpoi_interval)

    N1 = Array{Int32}(undef, npoi_interval)
    too_short_time = Array{Bool}(undef, npoi_interval)
    for i in 1:npoi_interval
        k = kpoi_interval[i]
        too_short_time[i] = 1
        for k2 in k:nstep-1
            if z[1,k2] < q0 && z[1,k2+1] >= q0
                N1[i] = k2-k
                too_short_time[i] = 0
                break
            end
        end
    end

    return N1, too_short_time
end
