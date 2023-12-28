include("RenderFractal.jl")
include("../spiders/Spiders.jl")
include("../trees/HubbardTrees.jl")

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

function showmandelbrot(A::Complex, B::Complex, scale::Real)
    M = mandelbrot_patch(A,B,scale)
    PA = mproblem_array(M,escape(4),100)
    pic = escapetime(PA)
    return heatmap(pic)
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

    M = mandelbrot_patch(c1,c2,2)
    PA = mproblem_array(M,escape(4),500)
    return @time heatmap(escapetime.(PA))

end