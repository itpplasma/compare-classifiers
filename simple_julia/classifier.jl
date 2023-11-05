include("util.jl")

function classify(krec, maxorder=size(krec,1))
    for order in 1:maxorder
        recnum = get_recurrence_numbers(krec, order)
        if (all(recnum .<= 0))
            break
        end
        difference = recnum[recnum.>0] .- recnum[1]
        if !all(x -> x in (0,1), difference)
            return false
        end
    end
    return true
end


function find_recurrences(q)
    nstep = length(q)

    kpoi_interval = findall(qi -> is_x3_between_x1_and_x2(q[1], q[2], qi), q)
    npoi_interval = length(kpoi_interval)

    if npoi_interval < 2
        return Array{Int32}(undef, 0, 0)
    end
    krec = Array{Int32}(undef, length(q), npoi_interval-1)
    krec .= 0
    for i in 1:npoi_interval-1
        kfirst = kpoi_interval[i+1]
        krec[1, i] = kfirst
        counter = 1
        for k2 in kfirst+1:nstep
            if is_x3_between_x1_and_x2(q[1], q[2], q[k2])
                counter+=1
                krec[counter, i] = k2
            end
        end
    end

    return krec
end


function get_recurrence_numbers(krec, order=1)
    recnum = Array{Int32}(undef, size(krec, 1), size(krec, 2))
    recnum .= 0
    for i in 1:size(krec, 1) - order
        recnum[i, :] = krec[i+order, :] - krec[i, :]
        if (all(recnum[i, :] .<= 0))
            break
        end
    end
    return recnum
end
