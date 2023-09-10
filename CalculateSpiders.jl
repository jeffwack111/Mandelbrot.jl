include("Spiders.jl")

θ = 1//62
info = SpiderInfo(θ)
frames = 100

list = [standard_legs(θ)]

for i in 1:frames-1
    SL = list[end]
    push!(list,spider_map(info,SL))
end




