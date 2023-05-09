function P_inv(z::Complex,λ::Complex,branch)
    return 2*(branch*sqrt(z/λ)-1)
end

function orbit(theta::Rational)
    orb = Rational{Int64}[]
    angle = theta
        while isempty(findall(x->x==angle,orb))
            push!(orb,angle)
            angle = angle*2
            angle = angle%1//1
        end
    joint = findall(x->x==angle,orb)[1]
    return orb,joint
end


function itinerary(theta::Real,orbit::Vector{<:Real})
    star_2 = theta/2
    star_1 = (theta+1)/2
    theta_itinerary = Any[]
    for angle in orbit
        if angle == star_1
            push!(theta_itinerary,1)
        elseif angle == star_2
            push!(theta_itinerary,2)
        elseif angle > star_2 && angle < star_1
            push!(theta_itinerary,"A")
        else
            push!(theta_itinerary,"B")
        end
    end
    return theta_itinerary
end

function region(point::Complex,boundary::Vector{<:Complex})
    x2 = real(point)
    y2 = imag(point)
    
    n_segments = length(boundary)-1

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
        return "A"
    else
        return "B"
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
