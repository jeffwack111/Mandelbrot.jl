include("RenderFractal.jl")

function julia_patch(center::Complex, right_center::Complex)
    #Overall strategy: 
    #we will first compute the patch at the correct scale and orientation, centered at the origin
    #we will then translate the patch to to correct location and return it

    top_center = 1.0*im*right_center
    #points from the origin to the top of the frame

    horizontal_axis = LinRange(-right_center,right_center,1000)
    vertical_axis = LinRange(-top_center,top_center,1000)

    origin_patch = transpose(horizontal_axis) .+ vertical_axis
    #we want a matrix whose i,jth element is H[i] + V[j]

    return origin_patch .+ center
end

J = julia_patch(0.0+0.0im,2.0+0.0im)
f(z) = z*z - 1
PA = problem_array(J,f,4,100)
ET = escape_time.(PA)
C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
pic = apply_color(ET,C)