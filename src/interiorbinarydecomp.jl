function assignbinary((iter,z),destination,alpha)
    if isnan(iter)
        return NaN
    else
        turns = angle(z-destination)/(2*pi) + 0.5
        if alpha < 0.5
            if  turns >= alpha && turns <alpha+0.5
                return 1
            else
                return 0
            end
        else
            if  turns >= alpha || turns < alpha-0.5
                return 1
            else
                return 0
            end
        end
    end
end

function arrive(destination::Number,radius::Real)
    R = radius*radius
    function arrived(z::Number)
        if abs2(z-destination) < R
            return true
        else
            return false
        end
    end
    return arrived
end




function alphamovie(frames::Int)

    #From page 40 of "The Beauty of Fractals"
    epsilon = 0.0001
    #We will consider a point arrived at zero if its modulus is less than epsilon, and arrived at infinity if its modulus is greater than 1/epsilon

    #first, set up a grid of test points.
    J = julia_patch(0.0+0.0im,2.0+0.0im)

    c = -0.3+0.0im

    #find the stable fixed point
    star_1 = (1+sqrt(1-4*c))/2
    star_2 = (1-sqrt(1-4*c))/2

    if abs2(star_1*2) < 1
        print("$star_1 is the attracting fixed point")
        star = star_1
    elseif abs2(star_2*2) < 1
        print("$star_2 is the attracting fixed point")
        star = star_2
    else
        print("no attracting fixed point!")
    end


    maxiter = 100
    f(z) = z*z + c
    PA = jproblem_array(J,f,arrive(star,epsilon),maxiter)
    ET = escapetime.(PA)

    fig = Figure()
    scene = Scene(camera=campixel!, resolution=(1000,1000))
    ax = Axis(scene,aspect = 1)
    hidedecorations!(ax)
    hidespines!(ax)

    record(scene,"test.gif",1:frames; framerate = 10) do ii
        empty!(ax)
        pic = assignbinary.(ET,star,ii/frames)
        heatmap!(scene,pic,colormap = :grayC, nan_color = RGBf(152/255, 109/255, 227/255))
    end 
end