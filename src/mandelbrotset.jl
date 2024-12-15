function mandelbrotpatch(A::Complex, B::Complex, scale::Real)
    #we will first compute the patch at the correct scale and orientation, centered at the origin
    #we will then translate the patch to to correct location and return it

    center = 0.5*(A + B)  
    #The center of the patch

    to_side = scale*(A - B) 
    #The complex number which points from the origin to the right side of the origin-centered patch

    to_top = 1.0im*to_side
    #points from the origin to the top of the frame

    horizontal_axis = LinRange(-to_side,to_side,1000)
    vertical_axis = LinRange(-to_top,to_top,1000)

    origin_patch = transpose(vertical_axis) .+ horizontal_axis
    #we want a matrix whose i,jth element is H[i] + V[j]

    return origin_patch .+ center
end

function showmandelbrot((A, B), scale::Real)
    M = mandelbrot_patch(A,B,scale)
    PA = mproblem_array(M,escape(100),100)
    pic = [x[1] for x in escapetime.(PA)]
    return heatmap(pic,nan_color = RGBAf(0,0,0,1),colormap = :PRGn_9)
end

function showm(angle1::Rational,angle2::Rational)
    @time c1 = spideriterate(angle1,500)
    println(c1)
    @time c2 = spideriterate(angle2,500)
    println(c2)

    M = mandelbrotpatch(c1,c2,10)
    PA = mproblem_array(M,escape(100),500)

    fig = Figure()
    scene = Scene(camera=campixel!, resolution=(1000,1000))
    ax = Axis(scene,aspect = 1)
    hidedecorations!(ax)
    hidespines!(ax)

    @time heatmap!(scene,mod.([x[1] for x in escapetime.(PA)],50),nan_color = RGBAf(0,0,0,1),colormap = :Set3_12)

    return scene
end


function showm(aia::AngledInternalAddress)
    @time theta1 = first(anglesof(aia))
    println(theta1)
    @time theta2 = first(anglesof(bifurcate(aia)))
    println(theta2)
    return showm(theta1,theta2)
end
    
function showm(theta::Rational)
    aia = AngledInternalAddress(theta)
    @time theta2 = angleof(first(criticalanglesof(bifurcate(aia))))
    println(theta2)
    return showm(theta,theta2)
end

function movie(list)
    for (ii,item) in enumerate(list)
        save("frameb$ii.png",showm(item))
    end
end

function overlay(maxiter::Int)
    M = mandelbrotpatch(0+0im,-2+0im,0.2)

    MA = mproblem_array(M,escape(100),maxiter)
    LA = lproblem_array(M,escape(10000),maxiter)

    EMA = escapetime.(MA)
    ELA = escapetime.(LA)

    CEMA = zeros(size(EMA))

    for ii in eachindex(EMA)
        if isnan(EMA[ii][1])
            CEMA[ii] = 1
        end
    end

    CELA = zeros(size(ELA))

    for ii in eachindex(ELA)
        if isnan(ELA[ii][1])
            CELA[ii] = 1
        end
    end

    fig = Figure()
    scene = Scene(camera=campixel!, resolution=(1000,1000))
    ax = Axis(scene,aspect = 1)

    #heatmap!(scene,CEMA, colormap = [RGBAf(c.r,c.g,c.b,0.5) for c in to_colormap(:reds)] )
    heatmap!(scene,CELA,colormap = [RGBAf(c.r,c.g,c.b,0.5) for c in to_colormap(:blues)])

    return scene

end