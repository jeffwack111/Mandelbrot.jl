function internal_address(K::PeriodicSequence)
    S = [1]
    v = PeriodicSequence([1])
    for ii in eachindex(K)
        if K[ii] != v[ii]
            push!(S,ii)
            v = PeriodicSequence(K[1:ii])
        end
    end
    return S
end

function kneading_sequence(A::Vector)
    if A == [1]
        return PeriodicSequence([1])
    else
        s = pop!(A)
        K = kneading_sequence(A)
        R = K[1:s-1]
        if K[s] == 0
            push!(R,1)
        else
            push!(R,0)
        end
        return PeriodicSequence(R)
    end
end
