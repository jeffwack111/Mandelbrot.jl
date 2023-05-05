using Plots

include("Spiders.jl")

K = KneadingSequence(1//240)
frames = 20

L = [standardlegs(K.orbit)]


for i in 1:frames-1
    push!(L,push!(spider_map(L[end],K)...))
end


@userplot SpiderPlot
@recipe function f(sp::SpiderPlot)
    L, i = sp.args
    xlim --> (-4,2)
    ylim --> (-3,3)
    aspect_ratio --> 1
    real(L[i]), imag(L[i])
end

anim = @animate for i ∈ 1:frames
    spiderplot(L, i)
end

gif(anim, "anim_fps1.gif", fps = 1)

