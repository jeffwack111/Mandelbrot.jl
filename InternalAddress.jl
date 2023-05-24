function internal_address(K::Vector)
    S = [1]
    v = ["A"]
    for ii in eachindex(K)
        if K[ii] != v[mod1(ii,end)]
            v = K[1:ii]
            push!(S,ii)
        end
    end
    return S
end

