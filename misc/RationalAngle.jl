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