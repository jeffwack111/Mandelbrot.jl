using CairoMakie
include("HubbardTrees.jl")
include("EmbedTrees.jl")

function showadjlist(E,root,labels,nodecolors=[])
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
    pos = Point.(X,Y)
    
    f = Figure()

    ax1 = Axis(f[1, 1],aspect = 1, limits = (0,1,0,1))
    hidedecorations!(ax1)
    hidespines!(ax1)

    for (ii,p) in enumerate(E)
        for n in p
            if isempty(nodecolors)
                lines!([pos[ii],pos[n]])
            else
                if nodecolors[ii] == "red" || nodecolors[n] == "red"
                    lines!([pos[ii],pos[n]],color = "red")
                elseif nodecolors[ii] == "blue" || nodecolors[n] == "blue"
                    lines!([pos[ii],pos[n]],color = "blue")
                elseif nodecolors[ii] == "orange" || nodecolors[n] == "orange"
                    lines!([pos[ii],pos[n]],color = "orange")
                elseif nodecolors[ii] == "green" || nodecolors[n] == "green"
                    lines!([pos[ii],pos[n]],color = "green")
                elseif nodecolors[ii] == "orangered2" || nodecolors[n] == "orangered2"
                    lines!([pos[ii],pos[n]],color = "orangered2")
                elseif nodecolors[ii] == "turquoise" || nodecolors[n] == "turquoise"
                    lines!([pos[ii],pos[n]],color = "turquoise")
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


function drawtree(H::Dict) 
    (E,nodes) = adjlist(H)
    root = filter(x->x.items[1]=='*',collect(keys(H)))[1]
    
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

function drawtree(H::Dict,OZ,boundary) 
    (E,nodes) = adjlist(H)
    #We need to make colorlist and edge color matrix here
    root = filter(x->x.items[1]=='*',collect(keys(H)))[1]
    
    criticalorbit = orbit(root)
    
    labels = []
    nodecolors = []
    zeroindex = findone(x->x.items[1] == '*',boundary)
    for node in nodes
        if node in boundary
            #determine if it is before, after, or equal to zero
            nodeindex = findone(x->x == node,boundary)
            if nodeindex < zeroindex #then we are in region B
                push!(nodecolors,"orangered2")
            elseif nodeindex == zeroindex #then this is zero
                push!(nodecolors,"black")
            else #then we are in region A
                push!(nodecolors,"turquoise")
            end
        else #we are fully in one of the 4 regions
            neighb = H[node][1]
            if node.items[1] == 'A'
                if OZ[node][neighb] == '0'
                    push!(nodecolors,"blue")
                else
                    push!(nodecolors,"green")
                end
            else
                if OZ[node][neighb] == '0'
                    push!(nodecolors,"red")
                else
                    push!(nodecolors,"orange")
                end
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
    return showadjlist(E,rootindex,labels,nodecolors)
end

function drawtree(angle::Rational)
    (OH,b) = embed(AngledInternalAddress(angle))
    OZ = labelonezero(OH,b)
    return drawtree(OH,OZ,b)
end

function drawtree(K::Sequence)
    return drawtree(hubbardtree(K))
end

