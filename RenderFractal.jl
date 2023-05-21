using Images
using ColorSchemes
using ImageView

function mandelbrot_patch(A::Complex, B::Complex, scale::Real)
    #Overall strategy: 
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

function escape_time(z0::Number,f::Function,radius_squared::Real,max_iter::Int)
    #here f must be a function of a single variable
    z = z0
    for iter in 1:max_iter
        if abs2(z) < radius_squared
            z = f(z)
        else
            return (iter - 1) % UInt8
        end
    end
    return max_iter % UInt8
end

function julia_pic(patch::Matrix,f::Function,radius_squared::Real,max_iter::Int)
    pic = zeros(RGB{Float64},size(patch))
        for i in eachindex(patch)
            pic[i] = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256][escape_time(patch[i],f,radius_squared,max_iter)+1]
        end
    return pic
end

function mandelbrot_pic(A::Complex, B::Complex, scale::Real, max_iter::Int)
    patch = mandelbrot_patch(A::Complex, B::Complex, scale::Real)
    pic = zeros(RGB{Float64},size(patch))
    for i in eachindex(patch)
        f(z) = z*z + patch[i]
        pic[i] = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256][escape_time(0,f,4,max_iter)+1]
    end
    return pic
end


