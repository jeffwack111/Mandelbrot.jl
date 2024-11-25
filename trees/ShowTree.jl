using GLMakie
include("HubbardTrees.jl")
include("OrientTrees.jl")
include("EmbedTrees.jl")
include("../spiders/Spiders.jl")
include("../parameters/DynamicRays.jl")

### Convenience functions
function plottree!(scene,angle::Rational)
    return plottree!(scene,OrientedHubbardTree(angle))
end

function plottree!(scene,K::Sequence)
    return plottree!(scene,HubbardTree(K))
end

function plottree(H)
    fig = Figure()
    ax = Axis(fig[1,1])
    return (fig,plottree!(ax,H))
end
###


function plottree!(scene, OHT::OrientedHubbardTree,steps = 8,colors = [])
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

    zvaluedict = embednodes(OHT)
    zvalues = []
    for node in Nodes
        push!(zvalues,zvaluedict[node])
    end

    E = refinedtree(OHT,zvaluedict,steps)

    pos = Point.(real.(zvalues)/2,imag.(zvalues)/2) #divide by 2 here for the scene coordinate system

    colorsforinterior = ["red","blue","green","orange"]

    for (ii,p) in enumerate(EdgeList)
        for n in p
            cmplxedge = E[Set([Nodes[ii],Nodes[n]])][2]
            realedge = Point.(real.(cmplxedge)/2,imag.(cmplxedge)/2) 
            if nodecolors[ii] in colorsforinterior 
                col = nodecolors[ii]
            elseif nodecolors[n] in colorsforinterior
                col = nodecolors[n]
            elseif nodecolors[ii] !== "black" 
                col = nodecolors[ii]
            elseif nodecolors[n] !== "black" 
                col = nodecolors[n]
            end
            lines!(scene,realedge,color = col,linewidth = 1,transparency = true,overdraw = true)
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
