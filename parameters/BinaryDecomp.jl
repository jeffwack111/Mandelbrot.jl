include("RenderFractal.jl")

#From page 40 of "The Beauty of Fractals"
epsilon = 0.001
#We will consider a point arrived at zero if its modulus is less than epsilon, and arrived at infinity if its modulus is greater than 1/epsilon

#first, set up a grid of test points.
J = julia_patch(0.0+0.0im,2.0+0.0im)

pic = zeros(RGB{Float64},size(J))

maxiter = 100
f(z) = z*z - 1

PA = problem_array(J,f,escape(1/epsilon),100)

binarycolor.(PA, blackwhite(0))

#=
#This needs to be less of a mess
for ii in 1:200
    save("binarydecomp$ii.png", binarycolor.(J, -0.12256+0.74486im,escape(1/epsilon),blackwhite(ii/200),maxiter))
end
=#