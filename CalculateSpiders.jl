include("Spiders.jl")

θ = 13//163
K = KneadingSequence(θ)
frames = 150

list = [SpiderLegs(K.orbit)]

for i in 1:frames-1
    SL = list[end]
    push!(list,spider_map(SL,K))
end




