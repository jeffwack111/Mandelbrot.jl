using CairoMakie
include("HubbardTrees.jl")
include("OrientTrees.jl")

function showadjlist(E,root,labels,nodecolors=[],pos=[])

    if isempty(pos)
        pos = generationposition(E,root)
    end
    
    f = Figure()

    ax1 = Axis(f[1, 1],aspect = 1, limits = (0,1,0,1))
    hidedecorations!(ax1)
    hidespines!(ax1)

    colorsforinterior = ["red","blue","green","orange"]

    for (ii,p) in enumerate(E)
        for n in p
            if isempty(nodecolors)
                lines!([pos[ii],pos[n]])
            else
                if nodecolors[ii] in colorsforinterior 
                    lines!([pos[ii],pos[n]],color = nodecolors[ii])
                elseif nodecolors[n] in colorsforinterior
                    lines!([pos[ii],pos[n]],color = nodecolors[n])
                elseif nodecolors[ii] !== "black" 
                    lines!([pos[ii],pos[n]],color = nodecolors[ii])
                elseif nodecolors[n] !== "black" 
                    lines!([pos[ii],pos[n]],color = nodecolors[n])
                end
            end
        end
    end
    if isempty(nodecolors)
        scatter!(pos)
    else
        scatter!(pos,color = nodecolors)
    end
    tex = [node[2] for node in labels]
    text!(pos,text = tex)
    
    return f
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
            X[u] = ii/(k+1)
            Y[u] = 1 - gen/(Ngens+1)
        end
    end

    return Point.(X,Y)
end

function showtree(H::Dict)
    (E,nodes) = adjlist(H)
    return showtree(E,nodes)
end

function showtree(E,nodes) 
    
    root = filter(x->x.items[1]=='*',nodes)[1]
    
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
    return showadjlist(E,rootindex,labels)
end

function showtree(H::Dict,OZ)
    (E,nodes) = adjlist(H)
    return showtree(H,(E,nodes),OZ)
end

function showtree(H::Dict,(E,nodes),OZ,pos=[]) 
    
    #We need to make colorlist and edge color matrix here
    root = filter(x->x.items[1]=='*',collect(keys(H)))[1]
    
    criticalorbit = orbit(root)
    
    labels = []
    nodecolors = []
    for node in nodes
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

    for node in nodes
        idx = findall(x->x==node, criticalorbit.items)
        if isempty(idx)
            push!(labels,Pair(node,repr("text/plain",node)))
        else
            push!(labels,Pair(node,string(idx[1]-criticalorbit.preperiod-1)))
        end
    end

    rootindex = findall(x->x==root,nodes)[1]
    return showadjlist(E,rootindex,labels,nodecolors,pos)
end

function showtree(angle::Rational)
    (OH,b) = orientedtree(AngledInternalAddress(angle))
    OZ = labelonezero(OH,b)
    return showtree(OH,OZ)
end

function showtree(K::Sequence)
    return showtree(hubbardtree(K))
end

