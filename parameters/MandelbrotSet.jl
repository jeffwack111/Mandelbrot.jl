include("RenderFractal.jl")
include("../spiders/Spiders.jl")
include("../trees/HubbardTrees.jl")
using ColorSchemes 

function mandelbrot_patch(A::Complex, B::Complex, scale::Real)
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
    pic = escapetime.(PA)
    return heatmap(pic,nan_color = RGBAf(0,0,0,1),colormap = :PRGn_9)
end

function showm(intadd::Vector{Int},numerators)
    head = push!(copy(intadd),2*intadd[end])

    @time theta1 = angleof(binary(hubbardtree(intadd),numerators))
    println(theta1)
    @time theta2 = angleof(binary(hubbardtree(head),numerators))
    println(theta2)

    @time c1 = spideriterate(theta1,500)
    println(c1)
    @time c2 = spideriterate(theta2,500)
    println(c2)

    M = mandelbrot_patch(c1,c2,50)
    PA = mproblem_array(M,escape(100),500)

    fig = Figure()
    scene = Scene(camera=campixel!, resolution=(1000,1000))
    ax = Axis(scene,aspect = 1)
    hidedecorations!(ax)
    hidespines!(ax)

    @time heatmap!(scene,mod.(escapetime.(PA),50),nan_color = RGBAf(0,0,0,1),colormap = :Set3_12)

    return scene
end

function showm(intadd::Vector{Int})
    n = length(characteristicpoints(hubbardtree(intadd)))
    nums = fill(1,n)
    return showm(intadd,nums)
end

function movie(intaddlist::Vector{Vector{Int}})
    for (ii,intadd) in enumerate(intaddlist)
        save("frame$ii.png",showm(intadd))
    end
end

function anglepairlist(n::Int)
    list = []
    for ii in 1:n
        intadd = [1,2,2*(ii+2)]
        head = [1,2,2*(ii+2),4*(ii+2)]
        num = [1]
        theta1 = angleof(binary(hubbardtree(intadd),num))
        theta2 = angleof(binary(hubbardtree(head),num))
        push!(list,(theta1,theta2))
    end
    return list
end

function anglepair(intadd::Vector{Int},numerators)
    head = push!(copy(intadd),2*intadd[end])
    theta1 = angleof(binary(hubbardtree(intadd),numerators))
    theta2 = angleof(binary(hubbardtree(head),numerators))
    return (theta1,theta2)
end

function parameterpair(angles)
    c1 = spideriterate(angles[1],500)
    c2 = spideriterate(angles[2],500)
    return (c1,c2)
end


function oldmovie(frames::Int)

    fig = Figure()
    scene = Scene(camera=campixel!, resolution=(1000,1000))
    ax = Axis(scene,aspect = 1)
    hidedecorations!(ax)
    hidespines!(ax)

    record(scene,"test.gif",1:frames; framerate = 3) do ii
        empty!(ax)
        intadd = [1,2]
        push!(intadd,(ii+1)*2+1)
        theta1 = angleof(binary(hubbardtree(intadd),[1]))
        head = copy(intadd)
        push!(head,intadd[end]*2)
        theta2 = angleof(binary(hubbardtree(head),[1]))

        c1 = spideriterate(theta1,100+50*ii)
        c2 = spideriterate(theta2,100+50*ii)

        M = mandelbrot_patch(c1,c2,50)
        PA = mproblem_array(M,escape(10000),500)
        heatmap!(scene,mod.(escapetime.(PA),50),nan_color = RGBAf(0,0,0,1),colormap = :Set3_12)

    end 
end