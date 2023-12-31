using CairoMakie

struct EscapeTimeProblem
    z0::Number
    f::Function
    stop::Function
    maxiter::Int
end 

function escapetime(prob::EscapeTimeProblem)
    z = prob.z0
    for iter in 1:prob.maxiter
        if prob.stop(z) == true
            if imag(z) > 0
                return iter
            else
                return -iter
            end
        else
            z = prob.f(z)
        end
    end
    return NaN
end


function jproblem_array(patch::Matrix,f::Function,s::Function,maxiter::Int)
    PA = Array{EscapeTimeProblem}(undef,size(patch)...)
    for i in eachindex(patch)
        PA[i] = EscapeTimeProblem(patch[i],f,s,maxiter)
    end
    return PA
end

function mproblem_array(patch::Matrix,s::Function,maxiter::Int)

    PA = Array{EscapeTimeProblem}(undef,size(patch)...)
    for i in eachindex(patch)
        PA[i] = EscapeTimeProblem(patch[i],z->z*z+patch[i],s,maxiter)
    end
    return PA
end

function normalize_escape_time!(patch::Matrix)
    patch = patch .- patch[1,1]
end

function blackwhite(alpha::Real)

    function color(z::Complex)
        turns = angle(z)/(2*pi) + 0.5
        if  turns >= alpha/2 && turns <alpha/2+0.5
            return 1
        else
            return 0
        end
    end

    return color

end

function escape(Radius::Real)

    function escaped(z::Number)
        if abs2(z) > Radius
            return true
        else
            return false
        end
    end

    return escaped
end


function julia_patch(center::Complex, right_center::Complex)
    #Overall strategy: 
    #we will first compute the patch at the correct scale and orientation, centered at the origin
    #we will then translate the patch to to correct location and return it

    top_center = 1.0*im*right_center
    #points from the origin to the top of the frame

    horizontal_axis = LinRange(-right_center,right_center,1000)
    vertical_axis = LinRange(-top_center,top_center,1000)

    origin_patch = transpose(vertical_axis) .+ horizontal_axis

    return origin_patch .+ center
end



