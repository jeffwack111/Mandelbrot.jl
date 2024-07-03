using CairoMakie
include("HubbardTrees.jl")
include("OrientTrees.jl")
include("../spiders/Spiders.jl")
include("../parameters/DynamicRays.jl")

function generationposition(E,root)
    T = [[root],E[root]]
    
    nadded = 1 + length(T[2])
    n = length(E)

    while nadded < n
        parents = T[end]
        children = []
        for parent in parents
            #cycle E[parent] so that the grandparent is in the last position
            s = findone(x-> x in T[end-1], E[parent])
            for u in circshift(E[parent],-s)
                if !(u in T[end-1])
                    push!(children,u)
                end
            end
        end
        push!(T,children)
        nadded += length(children)
    end
    Ngens = length(T)

    X = zeros(Float64,n) 
    Y = zeros(Float64,n)
    
    for (gen,vertices) in enumerate(T)
        k = length(vertices)
        for (ii,u) in enumerate(vertices)
            X[u] = 2*ii/(k+1) -1
            Y[u] = 2*(1 - gen/(Ngens+1))-1
        end
    end

    return Point.(X,Y)
end

function plottree!(scene, OHT::OrientedHubbardTree)
    OZ = labelonezero(OHT)

    anglelist = allanglesof(OZ,OHT)

    (EdgeList,Nodes) = adjlist(OHT.adj)

    criticalorbit = orbit(OHT.zero)
    labels = []
    nodecolors = []
    for node in Nodes
        firstchar = node.items[1]
        if firstchar == '*'
            push!(nodecolors,"black")
        elseif firstchar == 'A' #we are fully in one of the 4 regions
            if OZ[node] == '0'
                push!(nodecolors,"blue")
            elseif OZ[node] == '*'
                push!(nodecolors,"turquoise")
            elseif OZ[node] == '1'
                push!(nodecolors,"green")
            end
        elseif firstchar == 'B'
            if OZ[node] == '0'
                push!(nodecolors,"red")
            elseif OZ[node] == '*'
                push!(nodecolors,"orangered2")
            elseif OZ[node] == '1'
                push!(nodecolors,"orange")
            end
        end
    end

    for node in Nodes
        idx = findall(x->x==node, criticalorbit.items)
        if isempty(idx)
            push!(labels,Pair(node,repr("text/plain",node)))
        else
            push!(labels,Pair(node,string(idx[1]-criticalorbit.preperiod-1)))
        end
    end

    theta = angleof(first(criticalanglesof(OZ,OHT)))
    c = spideriterate(theta,100)

    #critical orbit
    rays = Dict()
    #periodic branch points
    characteristicpoints = characteristicset(OHT)
    for point in characteristicpoints
        merge!(rays,dynamicrays(c,angleof(first(anglelist[point])),100,10,20))
    end

    #preperiod branch points. This is slow!
    for node in Nodes
        for angle in anglelist[node]
            if !(angle in keys(rays))
                theta = angleof(angle)
                push!(rays,Pair(theta,dynamicrays(c,theta,100,10,20)[theta]))
            end
        end
    end
    
    merge!(rays,dynamicrays(c,1//2,100,10,20))
    merge!(rays,dynamicrays(c,0//1,100,10,20))

    paramorbit = [0.0+0.0im]
    n = period(theta)
    for ii in 1:n-1
        push!(paramorbit,paramorbit[end]^2+c)
    end

    zvalues = []
    for node in Nodes
        if '*' in node.items #then it is in the critical orbit
            idx = findone(x->x=='*',node.items)
            push!(zvalues,paramorbit[mod1(n-idx+2,n)])
        else
            list = [rays[angleof(angle)][end] for angle in anglelist[node]]
            push!(zvalues,sum(list)/length(list))
        end
    end

    pos = Point.(real.(zvalues)/2,imag.(zvalues)/2) #divide by 2 here for the scene coordinate system

    colorsforinterior = ["red","blue","green","orange"]

    for (ii,p) in enumerate(EdgeList)
        for n in p
            if nodecolors[ii] in colorsforinterior 
                lines!(scene,[pos[ii],pos[n]],color = nodecolors[ii])
            elseif nodecolors[n] in colorsforinterior
                lines!(scene,[pos[ii],pos[n]],color = nodecolors[n])
            elseif nodecolors[ii] !== "black" 
                lines!(scene,[pos[ii],pos[n]],color = nodecolors[ii])
            elseif nodecolors[n] !== "black" 
                lines!(scene,[pos[ii],pos[n]],color = nodecolors[n])
            end
        end
    end

    scatter!(scene,pos,color = nodecolors)
    tex = [node[2] for node in labels]
    text!(scene,pos,text = tex)
    
    return scene
end

function plottree!(scene, H::HubbardTree)

    (E, nodes) = adjlist(H.adj)
    root = H.zero
    
    criticalorbit = orbit(root)
    
    labels = []
   
    for node in nodes
        idx = findall(x->x==node, criticalorbit.items)
        if isempty(idx)
            push!(labels,Pair(node,repr("text/plain",node)))
        else
            push!(labels,Pair(node,string(idx[1]-criticalorbit.preperiod-1)))
        end
    end

    rootindex = findall(x->x==root,nodes)[1]

    pos = generationposition(E,rootindex)

    for (ii,p) in enumerate(E)
        for n in p
            lines!(scene,[pos[ii],pos[n]],linewidth = 1)
        end
    end

    scatter!(scene,pos)
    tex = [node[2] for node in labels]
    text!(scene,pos,text = tex)
    return scene
end

function plottree!(scene,angle::Rational)
    return plottree!(scene,OrientedHubbardTree(angle))
end

function plottree(H)
    scene = Scene()
    return plottree!(scene,H)
end