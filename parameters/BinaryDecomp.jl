using Images
using ColorSchemes
include("JuliaSet.jl")

#From page 40 of "The Beauty of Fractals"
epsilon = 0.001
#We will consider a point arrived at zero if its modulus is less than epsilon, and arrived at infinity if its modulus is greater than 1/epsilon

#first, set up a grid of test points.
J = julia_patch(0.0+0.0im,2.0+0.0im)

pic = zeros(RGB{Float64},size(J))

maxiter = 100

function binarycolor(z::Complex,c::Complex,eps::Real,maxiter::Int)
    iter = 0

    while iter<maxiter && abs2(z) > epsilon && abs2(z) < 1/epsilon 
        z = z^2 + c
        iter += 1
    end
    
    if angle(z)/(2*pi) <= 0
        return RGB{Float64}(1,1,1)
    else
        return RGB{Float64}(0,0,0)
    end
end

pic = binarycolor.(J,-2.0+0.0im,epsilon,maxiter)