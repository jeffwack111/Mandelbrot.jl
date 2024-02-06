include("RenderFractal.jl")

function assignbinary((iter,z))
    if isnan(iter)
        return NaN
    else
        turns = angle(z)/(2*pi) + 0.5
        if  turns >= 0.5 
            return 1
        else
            return 0
        end
    end
end

#From page 40 of "The Beauty of Fractals"
epsilon = 0.01
#We will consider a point arrived at zero if its modulus is less than epsilon, and arrived at infinity if its modulus is greater than 1/epsilon

#first, set up a grid of test points.
J = julia_patch(0.0+0.0im,3.0+0.0im)

maxiter = 100
f(z) = z*z + 0.0+0.0im
PA = jproblem_array(J,f,escape(1/epsilon),maxiter)

pic = assignbinary.(escapetime.(PA))

fig = Figure()
scene = Scene(camera=campixel!, resolution=(1000,1000))
ax = Axis(scene,aspect = 1)
hidedecorations!(ax)
hidespines!(ax)

heatmap!(scene,pic,colormap = :grayC, nan_color = RGBf(13/255, 99/255, 0/255))
scene
