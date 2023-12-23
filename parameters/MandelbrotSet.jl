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

    origin_patch = transpose(horizontal_axis) .+ vertical_axis
    #we want a matrix whose i,jth element is H[i] + V[j]

    return origin_patch .+ center
end

function showmandelbrot(A::Complex, B::Complex, scale::Real)
    M = mandelbrot_patch(A,B,scale)
    PA = mproblem_array(M,escape(4),100)
    C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
    colors = fill(C,size(PA))
    pic = escape_time.(PA,colors)
    return pic
end

function showm(intadd::Vector{Int})
    head = push!(copy(intadd),2*intadd[end])
    print(head)
    H1 = hubbardtree(intadd)
    H2 = hubbardtree(head)
    n = length(characteristicpoints(H1))
    numerators = fill(1,n)

    theta1 = angleof(binary(H1,numerators))
    theta2 = angleof(binary(H2,numerators))
    println(theta1)
    println(theta2)

    c1 = spideriterate(theta1,100)
    c2 = spideriterate(theta2,100)
    println(c1)
    println(c2)

    M = mandelbrot_patch(c1,c2,4)
    PA = mproblem_array(M,escape(2),100)
    C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
    colors = fill(C,size(PA))
    return escape_time.(PA,colors)

end