@userplot SpiderPlot
@recipe function f(sp::SpiderPlot)
    list, i = sp.args
    xlim --> (-20,20)
    ylim --> (-20,20)
    legend --> false
    aspect_ratio --> 1
    real(list[i].legs), imag(list[i].legs)
    linewidth --> 3
    linecolor --> :red
    real(list[i].boundary), imag(list[i].boundary)

end

anim = @animate for i ∈ 1:frames
    spiderplot(list, i)
end

gif(anim, "anim_fps1.gif", fps = 1)