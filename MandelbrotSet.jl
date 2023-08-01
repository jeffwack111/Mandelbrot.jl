include("RenderFractal.jl")

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

M = mandelbrot_patch(0.38693644424000584 + 0.09839851183492665im,0.38693644424000806 + 0.09839851183492719im,6)
PA = problem_array(M,4,1000)
ET = escape_time.(PA)
C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
pic = apply_color(ET,C)