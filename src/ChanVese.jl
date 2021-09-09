module ChanVese
using DistanceTransforms
using Images
using ImageFiltering

include("initializers.jl")
include("classical_chan_vese.jl")

export
    checkerboard,
    disk
    classical_chan_vese

end
