include("Spiders.jl")

K = KneadingSequence(333//1604)
frames = 20

list = [SpiderLegs(K.orbit)]

for i in 1:frames-1
    SL = list[end]
    push!(list,spider_map(SL,K))
end




