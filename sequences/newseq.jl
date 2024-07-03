abstract type Seq{T} end

#can only make things which isbits() into type parameters.

struct Seque{T}
    
    items::Vector

    function Seq(items::Vector,alp::Set) 
        for item in items
            if !(item in alp)
                error("$item is not in $alp")
            end
        end
        new{alp}(items)
    end

end

#This is promising!
struct Ary{N}
    digits::Vector
    function Ary{N}(digits::Vector) where N
        for x in digits
            if x<0 || x>=N
                error("$x is out of range $N")
            end
        end
        new{N}(digits)
    end
end