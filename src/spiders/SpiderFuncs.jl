using IterTools  

function region(point::Complex,boundary::Vector{<:Complex})
    x2 = real(point)
    y2 = imag(point)

    intersections = 0 

    for pair in partition(boundary,2,1)
        x3 = real(pair[1])
        y3 = imag(pair[1])
        
        x4 = real(pair[2])
        y4 = imag(pair[2])

        tn = y3*(x3-x4) - x3*(y3-y4)
        un = x3*y2-y3*x2
        d = y2*(x3-x4) - x2*(y3-y4)

        if sign(tn)*sign(d) == -1
            intersections += 0
        elseif sign(un)*sign(d) == -1
            intersections += 0
        elseif abs(tn) > abs(d)
            intersections += 0
        elseif abs(un) > abs(d)
            intersections += 0
        else
            intersections += 1
        end
    end

    if intersections%2 == 0
        return 'A'
    else
        return 'B'
    end

end

function cross_cut(z1::Complex,z2::Complex)
    if sign(imag(z1)) == sign(imag(z2))
        return 1
    else
        if sign(real(z1)) == sign(real(z2))
            if real(z1) >= 0
                return 1
            else
                return -1
            end
        else
            #At this stage what we know is that both signs change
            #Now our points are sorted, so the line from a to b crosses the y-axis from left to right

            if real(z1) < real(z2)
                a = z1
                b = z2
            else
                a = z2
                b = z1
            end
            
            diff = b-a
            rise = imag(diff)
            run = real(diff)
            x = real(a)
            y = imag(b)

            #The line segment crosses the cut when the magnitude of the slope is larger than the slope of the line 
            #connecting the first point to the origin
            if abs(rise/run)>=abs(y/x)
                return -1
            else
                return 1
            end    
        end
    end
end

function lift_path(path::Vector{<:Complex},λ::Complex,branch)

    lift = [P_inv(path[1],λ,branch)]

    for pair in partition(path,2,1)
        #does the segment divided by lambda cross the cut?
        a = pair[1]/λ
        b = pair[2]/λ
        branch = branch * cross_cut(a,b)
        append!(lift,P_inv(pair[2],λ,branch))
    end
    
    return lift

end

function path_sqrt(path::Vector{ComplexF64})
    branch = 1
    newpath = [sqrt(path[1])]
    cut = (0.0 + 0.0im, -1000 + 0.0im)
    for segment in partition(path,2,1)
        if test_intersection(cut...,segment...)
            branch *= -1
        end
        append!(newpath,branch*sqrt(segment[2]))
    end
    return newpath
end

function test_intersection(z1::Complex,z2::Complex,w1::Complex,w2::Complex)
    
    x1 = real(z1)
    y1 = imag(z1)

    x2 = real(z2)
    y2 = imag(z2)

    x3 = real(w1)
    y3 = imag(w1)

    x4 = real(w2)
    y4 = imag(w2)

    tn = (x1 - x3)*(y3 - y4) - (y1 - y3)*(x3 - x4)
    un = (x1 - x3)*(y1 - y2) - (y1 - y3)*(x1 - x2)
    d = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)

    if d == 0
        return false
    elseif sign(tn)*sign(d) == -1
        return false
    elseif sign(un)*sign(d) == -1
        return false
    elseif abs2(tn) > abs2(d)
        return false
    elseif abs2(un) > abs2(d)
        return false
    else
        return true
    end

end

#function angle(z1::Complex,z2::Complex,z3::Complex)
    #calculate the angle at z2 formed by the segments (z1,z2) and (z2,z3)
#end

