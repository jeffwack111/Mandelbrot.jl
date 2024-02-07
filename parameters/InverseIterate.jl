function inverseiterate(c::Complex,steps::Int)
    preimages = [3]
    for ii in 1:steps
        newpreimages = [] 
        for point in preimages
            push!(newpreimages,sqrt(point-c))
            push!(newpreimages,-sqrt(point-c))
        end
        preimages = newpreimages
    end
    return preimages
end