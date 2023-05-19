using Images
using ColorSchemes

function patch(A::Complex, B::Complex, scale::Real)
    #Overall strategy: 
    #we will first compute the patch at the correct scale and orientation, centered at the origin
    #we will then translate the patch to to correct location and return it

    center = 0.5*(A + B) 
    #The center of the patch

    to_side = scale*(A - B) 
    #The complex number which points from the origin to the right side of the origin-centered patch

    to_top = -1.0im*to_side
    #points from the origin to the top of the frame

    horizontal_axis = LinRange(-to_side,to_side,1000)
    vertical_axis = LinRange(-to_top,to_top,1000)

    origin_patch = transpose(x) .+ y
    #we want a matrix whose i,jth element is H[i] + V[j]

    return origin_patch .+ center
end

function escape_time(z::Number,f::Function,radius_squared::Real,max_iter::Int)
    #here f must be a function of a single variable
    iter_count::UInt8 = 0
    while iter_count < max_iter && abs2(z) < radius_squared
        z = f(z)
        try
            iter_count += 1
        catch
            iter_count = 0
        end
    end
    return iter_count
end



