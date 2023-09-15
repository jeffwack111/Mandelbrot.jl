include("Spiders.jl")

θ = 3//202
info = SpiderInfo(θ)
frames = 50

list = [standardlegs(θ)]

for i in 1:frames-1
    SL = list[end]
    push!(list,spidermap(info,SL))
end




