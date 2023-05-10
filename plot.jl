using Plots

include("Spiders.jl")

K = KneadingSequence(1//6)
frames = 20

list = [SpiderLegs(K.orbit)]


for i in 1:frames-1
    SL = list[end]
    push!(list,spider_map(SL,K))
end


@userplot SpiderPlot
@recipe function f(sp::SpiderPlot)
    list, i = sp.args
    xlim --> (-5,3)
    ylim --> (-4,4)
    aspect_ratio --> 1
    real(list[i].legs), imag(list[i].legs)
end

anim = @animate for i ∈ 1:frames
    spiderplot(list, i)
end

gif(anim, "anim_fps1.gif", fps = 1)

