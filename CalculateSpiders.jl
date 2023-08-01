include("Spiders.jl")

θ = 4//17
K = kneading_sequence(θ)
frames = 300

list = [SpiderLegs(Orbit(θ))]

for i in 1:frames-1
    SL = list[end]
    push!(list,spider_map(SL,K))
end




