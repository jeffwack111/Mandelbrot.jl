include("RenderFractal.jl")
#using CairoMakie

J = julia_patch(0.0+0.0im,2.0+0.0im)
f(z) = z*z - 1
PA = problem_array(J,f,escape(4),100)
C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
colors = fill(C,size(PA))
pic = escape_time.(PA,colors)

#heatmap(pic)