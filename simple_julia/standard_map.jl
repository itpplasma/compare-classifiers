
function standard_map(q, p, K=0.6)
    pnew = mod(p + K*sin(q), 2*pi)
    qnew = mod(q + pnew, 2*pi)
    return qnew, pnew
end
