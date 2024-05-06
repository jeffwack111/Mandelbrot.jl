using CairoMakie
include("HubbardTrees.jl")

function showorientedtree(E,root,labels)
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
            lines!([pos[ii],pos[n]])
        end
    end

    scatter!(pos)
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
    return showorientedtree(E,rootindex,labels)
end

function drawtree(angle::Rational)
    return drawtree(embed(AngledInternalAddress(angle)))
end



