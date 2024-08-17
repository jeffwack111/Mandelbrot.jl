using CairoMakie
include("SpiderFuncs.jl")

A = 1.0-3.0im
B = -2.0+1.0im

C = -0.75-1.0im
D = sqrt(C)

seg = [A,B]
path = collect(LinRange(A,B,20))
segimage = path_sqrt(seg)
refinedsegimage = collect(LinRange(segimage[1],segimage[end],20))
segpreimage = [z*z for z in refinedsegimage]
pathimage = path_sqrt(path)

f = Figure()
ax = Axis(f[1,1])
lines!(ax,real(path),imag(path))
lines!(ax,real(pathimage),imag(pathimage))
lines!(ax,real(segimage),imag(segimage))
lines!(ax,real(segpreimage),imag(segpreimage))
scatter!(ax,real([C,D]),imag([C,D]))

f