module ChanVese
using DistanceTransforms
using Images
using ImageFiltering

include("initializers.jl")
include("classical_chan_vese.jl")

export
    # initializers.jl
    checkerboard,
    disk

    # classical_chan_vese.jl
    classical_chan_vese

end
