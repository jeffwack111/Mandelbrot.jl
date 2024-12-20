module Mandelbrot

export RationalAngle,
       BinaryExpansion, 
       KneadingSequence, 
       InternalAddress, 
       HubbardTree, 
       AngledInternalAddress, 
       OrientedHubbardTree, 
       parameter,
       treeplot,
       spiderplot,
       brotplot

using ColorSchemes,
      Colors,
      GLMakie,
      IterTools,
      Primes

include("Sequences.jl")
include("angledoubling.jl")

include("Graphs.jl")
include("HubbardTrees.jl")

include("orienttrees.jl")
include("dynamicrays.jl")
include("embedtrees.jl")

include("interiorbinarydecomp.jl")
include("juliaset.jl")
include("lamination.jl")
include("mandelbrotset.jl")

include("spiderfuncs.jl")
include("spidermap.jl")

include("renderfractal.jl")
include("perturb.jl")

include("showrays.jl")
include("showspider.jl")
include("showtree.jl")

end
