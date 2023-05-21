include("RenderFractal.jl")

J = julia_patch(0.0+0.0im,2.0+0.0im)
f(z) = z*z - 1
PA = problem_array(J,f,4,100)
ET = escape_time.(PA)
C = colorschemes[:glasbey_bw_minc_20_hue_150_280_n256].colors
pic = apply_color(ET,C)