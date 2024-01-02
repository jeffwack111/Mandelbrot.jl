abstract type RationalAngle end

struct PeriodicRationalAngle <: RationalAngle
    angle::Rational   
end

struct PrePeriodicRationalAngle <: RationalAngle
    angle::Rational   
end

function RationalAngle(theta::Rational)
    num = theta.num%theta.den
    if theta.den%2 == 1
        return PeriodicRationalAngle(num//theta.den)
    else 
        return PrePeriodicRationalAngle(num//theta.den)
    end
end

#this seems wrong because the below struct in a kneading sequence, not an internal address
#=
struct IntAdd
    K::Sequence
end

function Base.iterate(I::IntAdd)
    return (1,1)
end

function Base.iterate(I::IntAdd,x)
    z = rho(I.K,x)    
    if isnothing(z)
        return z
    else
        return (z,z)
    end
end
=#

#=
function Base.in(I::IntAdd,n::Int)
    next = iterate(I)
    while next !== nothing && next[1] < n
        next = iterate(I, next[1])
    end

    if next == n
        return true
    else
        return false
    end

end
=#