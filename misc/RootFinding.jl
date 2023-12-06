using Symbolics
using Primes

#This seems to tap out at order 8, where it returns a polynomial with teeny tiny coefficients

@variables z,c

function composite(order::Int)
    if order==1
        return z^2+c
    else
        substitute(composite(order-1),z=>z^2+c)
    end
end

function divisors(n)
    d = Int64[1]
    for (p, e) in factor(n)
        r = 1
        l = length(d)
        for i in 1:e
            r *= p
            for j in 1:l
                push!(d, d[j]*r)
            end
        end
    end
    return sort(d)
end

function gleason(order::Int)
    if order==1
        return c
    end
    f = substitute(composite(order),z=>0)
    D = divisors(order)
    pop!(D)
    for d in D
        f = simplify(f/gleason(d))
    end
    return f

end

