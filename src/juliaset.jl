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

function juliaframe!(scene,param,colors)
    J = julia_patch(0.0+0.0im,2.0+0.0im)
    f(z) = z*z + param
    PA = jproblem_array(J,f,escapeorconverge(100),100)
    pic = [x[1] for x in escapetime.(PA)]
    n_colors = length(colors[])
    heatmap!(scene,pic,nan_color = RGBAf(0,0,0,1),colormap = colors,interpolate = false,colorrange = (1,n_colors),colorscale = x->mod1(x,n_colors))
    return scene
end

function juliaframe(param)
    scene = Scene(size= (500,500))#,camera = campixel!)
    return juliaframe!(scene,param)
end

function binaryframe!(scene,param)
    epsilon = 0.0001
    #We will consider a point arrived at zero if its modulus is less than epsilon, and arrived at infinity if its modulus is greater than 1/epsilon

    #first, set up a grid of test points.
    J = julia_patch(0.0+0.0im,2.0+0.0im)

    maxiter = 100
    f(z) = z*z + param
    PA = jproblem_array(J,f,escape(1/epsilon),maxiter)

    pic = assignbinary.(escapetime.(PA))

    heatmap!(scene,pic,colormap = :grayC, nan_color = RGBf(172/255, 88/255, 214/255))
    return scene
end

function binaryframe(param)
    scene = Scene(size= (500,500),camera = campixel!)
    return binaryframe!(scene,param)
end