
function newtriple(X,(A,B,C))
    newA = 2*X*A+1
    newB = 2*X*B + A^2
    newC = 2*X*C + 2*A*B
    return (newA,newB,newC)
end

function coefs(X)
    list = []
    (A,B,C) = (1,0,0)
    for x in X
        push!(list,(A,B,C))
        (A,B,C) = newtriple(x,(A,B,C))
    end
    push!(list,(A,B,C))
    return list
end

function perturb(X,coefs,delta)
    delta2 = delta^2
    delta3 = delta^3
    for (ii,x) in enumerate(zip(X,coefs))
        (A,B,C) = x[2]
        Y = x[1]+delta*A + delta2*B + delta3*C
        if abs2(Y)>4
            return (ii,Y)
        end
    end
    return (NaN,Y)
end

pilot = [0.0+1.0im]
for ii in 1:100
    push!(pilot,pilot[end]^2+pilot[1])
end

C = coefs(pilot)

rightcenterdelta = 0.00001+0.000001im

patch = julia_patch(0.0+0.0im,rightcenterdelta)

PA = [(pilot,C,delta) for delta in patch]

results = [perturb(stuff...) for stuff in PA]

pic = [x[1] for x in results]

heatmap(pic,nan_color = RGBAf(0,0,0,1),colormap = :PRGn_9)
