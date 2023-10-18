using JSON

newentry = Dict(string(θ) => list[end].legs[2][1]/2.0)

open("RayDictionary.txt","r") do file
    d = JSON.parse(file)
end

merge!(d,newentry)

rm("RayDictionary.txt")

open("RayDictionary.txt","w") do file
    JSON.print(file,d,2)
end

