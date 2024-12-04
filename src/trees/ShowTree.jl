using GLMakie
using ColorSchemes
include("EmbedTrees.jl")
include("../spidermap/SpiderMap.jl")
include("../parameters/DynamicRays.jl")
include("../parameters/RenderFractal.jl")

@recipe(TreePlot,EHT) do scene #TODO add kwargs to toggle plotting the julia set and the external rays
    Theme()
end

function Makie.plot!(myplot::TreePlot)
    EHT = myplot.EHT[]
    (EdgeList,Nodes) = adjlist(EHT.adj)
    criticalorbit = orbit(EHT.criticalpoint)

    labels = []
    nodecolors = []
    for node in Nodes
        firstchar = node.items[1]
        if firstchar == star()
            push!(nodecolors,"black")
        elseif firstchar == A() #we are fully in one of the 4 regions
            if EHT.onezero[node] == Digit{2}(0)
                push!(nodecolors,"blue")
            elseif EHT.onezero[node] === nothing
                push!(nodecolors,"turquoise")
            elseif EHT.onezero[node] == Digit{2}(1)
                push!(nodecolors,"green")
            end
        elseif firstchar == B()
            if EHT.onezero[node] == Digit{2}(0)
                push!(nodecolors,"red")
            elseif EHT.onezero[node] === nothing
                push!(nodecolors,"orangered2")
            elseif EHT.onezero[node] == Digit{2}(1)
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

    zvalues = [EHT.vertices[node] for node in Nodes]

    pos = Point.(real.(zvalues),imag.(zvalues)) 

    colorsforinterior = ["red","blue","green","orange"]

    for (ii,p) in enumerate(EdgeList)
        for n in p
            cmplxedge = EHT.edges[Set([Nodes[ii],Nodes[n]])][2]
            realedge = Point.(real.(cmplxedge),imag.(cmplxedge)) 
            if nodecolors[ii] in colorsforinterior 
                col = nodecolors[ii]
            elseif nodecolors[n] in colorsforinterior
                col = nodecolors[n]
            elseif nodecolors[ii] !== "black" 
                col = nodecolors[ii]
            elseif nodecolors[n] !== "black" 
                col = nodecolors[n]
            end
            lines!(myplot,realedge,color = col,linewidth = 1,transparency = true,overdraw = true)
        end
    end

    julia = inverseiterate(EHT.parameter,20)
    scatter!(myplot,real(julia),imag(julia),markersize = 1,color = "black")


    rays = collect(values(EHT.rays))
    n = length(rays)
    for (j ,ray) in enumerate(rays)
        lines!(myplot,real(ray),imag(ray),color = get(ColorSchemes.rainbow, float(j)/float(n)))
    end

    scatter!(myplot,pos,color = nodecolors)
    tex = [node[2] for node in labels]
    text!(myplot,pos,text = tex)
    limits!(-2,2,-2,2)
    return myplot
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
