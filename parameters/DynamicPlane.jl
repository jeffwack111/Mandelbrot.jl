include("DynamicRays.jl")
include("RenderFractal.jl")
include("../spiders/Spiders.jl")
include("../trees/Lamination.jl")
include("../trees/ShowTree.jl")

function characteristicrays(theta::Rational) 

    OHT = OrientedHubbardTree(theta)

    OZ = labelonezero(OHT)

    anglelist = allanglesof(OZ,OHT)

    (EdgeList,Nodes) = adjlist(OHT.adj)
    
    c = spideriterate(theta,100)
    
    #critical orbit
    rays = dynamicrays(c,theta,100,10,20)
    merge!(rays,dynamicrays(c,conjugate(theta),100,10,20))

    #periodic branch points
    characteristicpoints = characteristicset(OHT)
    for point in characteristicpoints
        merge!(rays,dynamicrays(c,angleof(first(anglelist[point])),100,10,20))
    end

    merge!(rays,dynamicrays(c,1//2,100,10,20))
    merge!(rays,dynamicrays(c,0//1,100,10,20))

    #preperiod branch points. This is slow!
    for node in Nodes
        for angle in anglelist[node]
            if !(angle in keys(rays))
                theta = angleof(angle)
                push!(rays,Pair(theta,dynamicrays(c,theta,100,10,20)[theta]))
            end
        end
    end

    orbit = [0.0+0.0im]
    n = period(theta)
    for ii in 1:n-1
        push!(orbit,orbit[end]^2+c)
    end

    zvalues = []
    for node in Nodes
        if '*' in node.items #then it is in the critical orbit
            idx = findone(x->x=='*',node.items)
            push!(zvalues,orbit[mod1(n-idx+2,n)])
        else
            list = [rays[angleof(angle)][end] for angle in anglelist[node]]
            push!(zvalues,sum(list)/length(list))
        end
    end

    pos = Point.(real.(zvalues),imag.(zvalues))
 
    (fig,ax) = showtree(OHT.adj,(EdgeList,Nodes),OZ,pos)

    onlyrays = [pair[2] for pair in rays]
    #plotrays!(ax,onlyrays)

    julia = inverseiterate(c,22)
    #scatter!(ax,real(julia),imag(julia),markersize = 1,color = "black")
    
    #scatter!(ax,real(orbit),imag(orbit))

    return fig
end