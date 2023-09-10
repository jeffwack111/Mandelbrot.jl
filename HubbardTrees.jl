using NetworkLayout

struct Sequence
    items::Vector
    pointsto::Vector{Int}
end

function shift(seq::Sequence)

    if seq.pointsto[end] != 1
        newpointsto = Int[]
        #we conclude the topology is preperiodic
        for index in Iterators.drop(seq.pointsto,1)
            push!(newpointsto,index-1)
        end

        newitems = collect(Iterators.drop(seq.items,1))

        return Sequence(newitems,newpointsto)
    else
        newitems = empty(seq.items)
        for index in seq.pointsto
            push!(newitems,seq.items[index])
        end
        
        return Sequence(newitems,seq.pointsto)
    end  

end




#We will represent the hubbard tree as an adjacency matrix of marked points