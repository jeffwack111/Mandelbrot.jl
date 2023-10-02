using IterTools
using Primes

function farey_sequence(depth)
    if depth == 1
        return [0//1,1//1]
    else
        parents = farey_sequence(depth - 1)
        for parent in partition(parents,2,1)
            child_denominator = denominator(parent[1]) + denominator(parent[2])
            if child_denominator <= depth
                child_numerator = numerator(parent[1]) + numerator(parent[2])
                ratio = child_numerator//child_denominator
                push!(parents,ratio)
            end
        end
        return sort!(parents)
    end
end


function how_many(n)
    return length(findall(x->x.joint==1&&length(x.orbit)==n,F))
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

function p(n)
    if n == 1
        return 2
    end
    div = divisors(n)
    pop!(div)
    composite = sum(p,div)
    return 2^n - composite
end