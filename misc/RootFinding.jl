using Symbolics

function composite(order::Int)
    if order==1
        return z^2+c
    else
        substitute(composite(order-1),z=>z^2+c)
    end
end

function gleason(order::Int)
    return expand(substitute(composite(order),z=>0))
end

