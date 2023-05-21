using IterTools

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
        
