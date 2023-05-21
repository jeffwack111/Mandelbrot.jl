using Images
using ColorSchemes
using ImageView
include("MakePatch.jl")

struct EscapeTimeProblem
    z0::Number
    f::Function
    R2::Real
    max_iter::Int
end

function problem_array(patch::Matrix,f::Function,R2::Real,max_iter::Int)
    PA = Array{EscapeTimeProblem}(undef,size(patch)...)
    for i in eachindex(patch)
        PA[i] = EscapeTimeProblem(patch[i],f,R2,max_iter)
    end
    return PA
end

function problem_array(patch::Matrix,R2::Real,max_iter::Int)
    PA = Array{EscapeTimeProblem}(undef,size(patch)...)
    for i in eachindex(patch)
        f(z) = z*z + patch[i]
        PA[i] = EscapeTimeProblem(0.0+0.0im,f,R2,max_iter)
    end
    return PA
end


function escape_time(prob::EscapeTimeProblem)
    z = prob.z0
    for iter in 1:prob.max_iter
        if abs2(z) < prob.R2
            z = prob.f(z)
        else
            return (iter - 1) % UInt8
        end
    end
    return prob.max_iter % UInt8
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

function normalize_escape_time!(patch::Matrix)
    patch = patch .- patch[1,1]
end

function apply_color(patch::Matrix,colors::Vector{RGB{Float64}})
    pic = zeros(RGB{Float64},size(patch))
    for i in eachindex(patch)
        pic[i] = colors[patch[i]+1]
    end
    return pic
end

J = julia_patch(0.0+0.0im,2.0+0.0im)
f(z) = z*z - 1
PA = problem_array(J,f,4,100)
ET = escape_time.(PA)
C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
pic = apply_color(ET,C)