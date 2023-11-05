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
