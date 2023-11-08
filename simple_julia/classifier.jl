include("util.jl")

"""
    classify(krec, maxorder=size(krec,1))

Classify a set of points as ideal or non-ideal.

This function takes the output of find_recurrences and checks if it takes
the same number of maps for points return to the interval between the first two points of the sequence for the first, second, third, etc. time.

# Arguments
- `krec::Matrix{Int32}`: A matrix of size `number of maps x (npoi_interval-1)`,
  where `npoi_interval` is the number of points in the interval between
  q[0] and q[1] and `length(q)` is the number of iterations of the map.
- `maxorder::Int = size(krec,1)`: The maximum order of recurrence numbers to
  check inside `get_recurrence_numbers`.

# Returns
- `is_ideal::Bool`: True if the set of points is classified as ideal, false
  otherwise.
"""
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

"""
    find_recurrences(q)

Counts the number of iterative applications of a map until which a set of
points returns to its original interval.

This function first finds all points between the first two points of the
sequence `q` and then counts the number of iterations until which each of
these points returns to the interval between the first two points.

# Arguments
- `q::Vector{Float64}`: A sequence of points for subsequent application
  of a map. Initial conditions: q[0], after one iteration: q[1], etc.

# Returns
- `krec::Matrix{Int32}`: A matrix of size `length(q) x (npoi_interval-1)`,
  where `npoi_interval` is the number of points in the interval between
  q[0] and q[1]. The value of `krec` is the index
  of the point in the sequence `q` for which the point returns the `n`th time.
  The first dimension corresponds to the `n`-th time
  for which the point returns to the interval. The second dimension is
  the index of the point in the interval.

"""
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


"""
    get_recurrence_numbers(krec, order=1)

Calculate the number of mappings that each of some values.

Based on the output `krec` of find_recurrences, this function computes the
number of iterations a set of points takes to return into an interval for
the first, second, third, etc. time. In addition, an order parameter can be
specified to count `order` iterations as a single mapping.

# Arguments
- `krec::Matrix{Int32}`: A matrix of size `number of maps x (npoi_interval-1)`,
  where `npoi_interval` is the number of points in the interval between
  q[0] and q[1] and `length(q)` is the number of iterations of the map.
- `order::Int = 1`:

# Returns
- `recnum::Matrix{Int32}`: A matrix of the same size as `krec`, containing the
  calculated recurrence numbers, i.e the number of maps taken between the
  `n`-th and `n+1`-th visit of the initial interval. In case not enough points
  are available, the remaining entries are set to zero or negative values.

# Example
```julia


```
"""
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
