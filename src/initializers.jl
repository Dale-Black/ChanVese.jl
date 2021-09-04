"""
    checkerboard(shape)

Generates a checkerboard level set function.
According to Pascal Getreuer, such a level set function has fast
convergence.

# Arguments
- shape: 2D shape of the output array
"""
function checkerboard(shape)
    x₀ = reshape(collect(0:(shape[1] - 1)), shape[1], 1)
    y₀ = reshape(collect(0:(shape[2] - 1)), 1, shape[2])
    return 𝚽₀ = sin.(pi / 5 .* x₀) .* sin.(pi / 5 .* y₀)
end

"""
    disk(shape)

Generates a disk level set function.
The disk covers the whole image along its smallest dimension
or `factor` can be given to decrease the size of the disk.

# Arguments
- shape: 2D shape of the output array
- factor: size to decrease the disk, default is
    1 which implies no size change
"""
function disk(shape, factor=1)
    res = ones(shape)
    centerY = Int(round((shape[1] - 1) / 2))
    centerX = Int(round((shape[2] - 1) / 2))
    res[centerY, centerX] = 0.0
    radius = Float32(min(centerX, centerY)) / factor
    return disk = (radius .- DistanceTransforms.euclidean(res)) / radius
end
