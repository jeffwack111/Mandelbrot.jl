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

function juliamovie(frames::Int)
    angles = LinRange(0,1,frames+1)[1:frames]
    params = [0.35*(1-cos(2*pi*angle))*exp(2*pi*1.0im*angle) for angle in angles]
    for (ii,param) in enumerate(params)
        save("frame$ii.png",binaryframe(param))
    end
end

function juliaframe(param)
    scene = Scene(camera=campixel!, resolution=(1000,1000))
    J = julia_patch(0.0+0.0im,2.0+0.0im)
    f(z) = z*z + param
    PA = jproblem_array(J,f,escapeorconverge(100),100)
    pic = [x[1] for x in escapetime.(PA)]
    heatmap!(scene,pic,nan_color = RGBAf(0,0,0,1),colormap = :dense)
    return scene
end

function binaryframe(param)
    epsilon = 0.0001
    #We will consider a point arrived at zero if its modulus is less than epsilon, and arrived at infinity if its modulus is greater than 1/epsilon

    #first, set up a grid of test points.
    J = julia_patch(0.0+0.0im,4.0+0.0im)

    maxiter = 100
    f(z) = z*z + param
    PA = jproblem_array(J,f,escape(1/epsilon),maxiter)

    pic = assignbinary.(escapetime.(PA))

    scene = Scene(camera=campixel!, resolution=(1000,1000))

    heatmap!(scene,pic,colormap = :grayC, nan_color = RGBf(172/255, 88/255, 214/255))
    return scene
end

#=
z_rabbit = -0.12256117 + 0.74486177im

J = julia_patch(0.0+0.0im,2.0+0.0im)
f(z) = z*z + 0.0+0.0im
PA = jproblem_array(J,f,escapeorconverge(100),100)
#C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
#colors = fill(C,size(PA))
pic = [x[1] for x in escapetime.(PA)]

scene = Scene(camera=campixel!, resolution=(1000,1000))
heatmap!(scene,pic,nan_color = RGBAf(0,0,0,1),colormap = :dense)
scene
=#