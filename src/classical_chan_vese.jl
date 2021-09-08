"""
    c1(img, HΦ)

Calculate the average value inside the function

# Arguments
- img: input image
- HΦ: Heaviside of the level set function `Φ`
"""
function c1(img, HΦ)
    ∫H𝚽 = sum(HΦ)
    ∫u₀H𝚽 = sum(img .* HΦ)
    return c₁ = ∫u₀H𝚽 / ∫H𝚽
end

"""
    c2(img, HΦ)

Calculate the average value outside the function

# Arguments
- img: input image
- HΦ: Heaviside of the level set function `Φ`
"""
function c2(img, HΦ)
    H𝚽ⁱ = 1 .- HΦ
    ∫H𝚽ⁱ = sum(H𝚽ⁱ)
    ∫u₀H𝚽ⁱ = sum(img .* H𝚽ⁱ)
    return c₂ = ∫u₀H𝚽ⁱ / ∫H𝚽ⁱ
end

"""
    δₕ(x, h=1)

Modified Dirac Delta function

# Arguments
- x: input
"""
δₕ(x, h=1) = h ./ (h .^ 2 .+ x .^ 2)

raw"""
	calculate_variation(img, 𝚽ⁿ, μ, λ₁, λ₂, Δt)

Calculates the variation of level set `𝚽ⁿ` based on algorithm parameters.
Solves the set of equations from the given paper below:

```math
\begin{aligned}
	\phi_{i+} 	&= \phi_{i, j+1}^n - \phi_{i, j}^n \\
	\phi_{i-} 	&= \phi_{i, j}^n - \phi_{i, j - 1}^n \\
	\phi_{i} 	&= \frac{\phi_{i, j+1}^n - \phi_{i, j-1}^n}{2} \\
	\phi_{j+} 	&= \phi_{i+1, j}^n - \phi_{i, j}^n \\
	\phi_{j-} 	&= \phi_{i, j}^n - \phi_{i-1, j}^n \\
	\phi_{j} 	&= \frac{\phi_{i+1, j}^n - \phi_{i-1, j}^n}{2} \\
\end{aligned}
\tag{9a}
```
```math
\begin{aligned}
	C_1 &= \frac{1}{\sqrt{\phi_{i+}^2 - \phi_{j}^2}} \\
	C_2 &= \frac{1}{\sqrt{\phi_{i-}^2 - \phi_{j}^2}} \\
	C_3 &= \frac{1}{\sqrt{\phi_{i}^2 - \phi_{j+}}} \\
	C_4 &= \frac{1}{\sqrt{\phi_{i}^2 - \phi_{j+}}}
\end{aligned}
\tag{9b}
```
```math
	C = C_1 * \phi_{i, j+1}^n + C_2 * \phi_{i, j-1}^n + C_3 * \phi_{i+1, j}^n + C_4 * \phi_{i-1, j}^n
\tag{9c}
```
```math
\begin{aligned}
D = -\lambda_{1} (u_{0_{i,j}} - c_1*(\phi_{i,j}^n))^2 + \lambda_{2} (u_{0_{i,j}} - c_2 * (\phi_{i,j}^n))^2
\end{aligned}
\tag{9d}
```
```math
\begin{aligned}
	\phi_{i,j}^{n+1} = \frac{\phi_{i,j}^n + \Delta t * \delta(\phi_{i,j}^n) * (\mu * C + D)}{ 1 + \mu * \Delta t * \delta(\phi_{i,j}^n) * (C_1 + C_2 + C_3 + C_4)}
\end{aligned}
\tag{9e}
```

# Arguments
- img: input image
- 𝚽ⁿ: level set function
- μ: "edge length" weight parameter. Higher `μ` values will
	produce a 'round' edge, while values closer to zero will
	detect smaller objects.
- λ₁: "difference from average" weight parameter for the output
	region with value 'true'. If it is lower than `λ₂`, this
	region will have a larger range of values than the other.
- λ₂: "difference from average" weight parameter for the output
	region with value 'false'. If it is lower than `λ₁`, this
	region will have a larger range of values than the other.
- Δt: A multiplication factor applied at calculations for each step,
	serves to accelerate the algorithm. While higher values may
	speed up the algorithm, they may also lead to convergence
	problems.

# Citation
'The Chan-Vese Algorithm' [Rami Cohen](https://arxiv.org/pdf/1107.2782.pdf)
"""
function calculate_variation(img, 𝚽ⁿ, μ, λ₁, λ₂, Δt)
    ϵ = 1e-16
    𝚽⁺ = padarray(𝚽ⁿ, Pad(1, 1))

    # 9a
    𝚽ᵢ₊ = 𝚽⁺[1:(end - 1), 2:end] .- 𝚽⁺[1:(end - 1), 1:(end - 1)]
    𝚽ᵢ₋ = 𝚽⁺[1:(end - 1), 1:(end - 1)] .- 𝚽⁺[1:(end - 1), 0:(end - 2)]
    𝚽ᵢ = (𝚽⁺[1:(end - 1), 2:end] .- 𝚽⁺[1:(end - 1), 0:(end - 2)]) / 2

    𝚽ⱼ₊ = 𝚽⁺[2:end, 1:(end - 1)] - 𝚽⁺[1:(end - 1), 1:(end - 1)]
    𝚽ⱼ₋ = 𝚽⁺[1:(end - 1), 1:(end - 1)] - 𝚽⁺[0:(end - 2), 1:(end - 1)]
    𝚽ⱼ = (𝚽⁺[2:end, 1:(end - 1)] - 𝚽⁺[0:(end - 2), 1:(end - 1)]) / 2

    # 9b
    C₁ = 1 ./ sqrt.(ϵ .+ 𝚽ᵢ₊ .^ 2 .+ 𝚽ⱼ .^ 2)
    C₂ = 1 ./ sqrt.(ϵ .+ 𝚽ᵢ₋ .^ 2 .+ 𝚽ⱼ .^ 2)
    C₃ = 1 ./ sqrt.(ϵ .+ 𝚽ᵢ .^ 2 .+ 𝚽ⱼ₊ .^ 2)
    C₄ = 1 ./ sqrt.(ϵ .+ 𝚽ᵢ .^ 2 .+ 𝚽ⱼ₋ .^ 2)

    # 9c
    C =
        𝚽⁺[1:(end - 1), 2:end] .* C₁ .+ 𝚽⁺[1:(end - 1), 0:(end - 2)] .* C₂ .+
        𝚽⁺[2:end, 1:(end - 1)] .* C₃ .+ 𝚽⁺[0:(end - 2), 1:(end - 1)] .* C₄

    H𝚽 = 1 .* (𝚽ⁿ .> 0)
    c₁, c₂ = c1(img, H𝚽), c2(img, H𝚽)

    # 9d
    D = -λ₁ .* (img .- c₁) .^ 2 .+ λ₂ .* (img .- c₂) .^ 2

    # 9
    return 𝚽ⁿ⁺¹ =
        (𝚽ⁿ .+ Δt .* δₕ(𝚽ⁿ) .* (μ * C .+ D)) ./
        (1 .+ μ .* Δt .* δₕ(𝚽ⁿ) .* (C₁ .+ C₂ .+ C₃ .+ C₄))
end
raw"""
	calculate_reinitial(𝚽, Δt)

Calculates the way in which to reinitialize the level set `𝚽` based on 
algorithm parameters. Solves the set of equations from the given paper below:

```math
\begin{aligned}
	a &= (\Psi_{i,j} - \Psi_{i-1,j}) \\
	b &= (\Psi_{i,j} - \Psi_{i+1,j}) \\
	c &= (\Psi_{i,j} - \Psi_{i,j-1}) \\
	d &= (\Psi_{i,j} - \Psi_{i,j+1})
\end{aligned}
\tag{12a}
```
```math
\begin{aligned}
	G(\Psi_{i,j}^{n}) &=
    \begin{cases}
		\sqrt{max((a^+)^2, (b^-)^2) + max((c^+)^2, (d^-)^2)} - 1, 
			\ &\Phi(x_i, j_i, t) > 0 \\
		\sqrt{max((a^-)^2, (b^+)^2) + max((c^-)^2, (d^+)^2)} - 1, 
			\ &\Phi(x_i, j_i, t) < 0 \\
      0 \ &\text{otherwise}
    \end{cases} 
\end{aligned}
\tag{12b}
```
```math
\begin{aligned}
	sign(\Phi(x,y,t)) &= \frac{\Phi(x,y,t)}{\sqrt{\Phi(x,y,t)^2}}
\end{aligned}
\tag{11a}
```
```math
\begin{aligned}

	\Psi_{i,j}^{n+1} &= \Psi_{i,j}^{n} - \Delta \tau * sign(\Phi(x,y,t) * G(\Psi_{i,j}^{n})
\end{aligned}
\tag{11b}
```

# Arguments
- 𝚽: level set function
- Δt: A multiplication factor applied at calculations for each step,
	serves to accelerate the algorithm. While higher values may
	speed up the algorithm, they may also lead to convergence
	problems.

# Citation
'The Chan-Vese Algorithm' [Rami Cohen](https://arxiv.org/pdf/1107.2782.pdf)
"""
function calculate_reinitial(𝚽, Δt)
    ϵ = 1e-8
    𝚽⁺ = padarray(𝚽, Pad(1, 1))

    # 12a
    a = 𝚽⁺[1:(end - 1), 1:(end - 1)] - 𝚽⁺[1:(end - 1), 0:(end - 2)]
    b = 𝚽⁺[1:(end - 1), 2:end] - 𝚽⁺[1:(end - 1), 1:(end - 1)]
    c = 𝚽⁺[1:(end - 1), 1:(end - 1)] - 𝚽⁺[0:(end - 2), 1:(end - 1)]
    d = 𝚽⁺[2:end, 1:(end - 1)] - 𝚽⁺[1:(end - 1), 1:(end - 1)]

    a⁺ = max(a, 0)
    a⁻ = min(a, 0)
    b⁺ = max(b, 0)
    b⁻ = min(b, 0)
    c⁺ = max(c, 0)
    c⁻ = min(c, 0)
    d⁺ = max(d, 0)
    d⁻ = min(d, 0)

    G = zeros(size(𝚽))
    index⁺ = 𝚽 .> 0
    index⁻ = 𝚽 .< 0

    # 12b
    G =
        (sqrt.(max.(a⁺ .^ 2, b⁻ .^ 2) + max.(c⁺ .^ 2, d⁻ .^ 2)) .- 1) .* index⁺ .+
        (sqrt.(max.(a⁻ .^ 2, b⁺ .^ 2) .+ max.(c⁻ .^ 2, d⁺ .^ 2)) .- 1) .* index⁻

    # 11a
    sign𝚽 = 𝚽 ./ sqrt.(𝚽 .^ 2 .+ ϵ)

    # 11b
    return 𝚿 = 𝚽 .- Δt .* sign𝚽 .* G
end

function reinitialize(𝚽, Δt, max_reiter=5)
    iter = 0
    while iter < max_reiter
        𝚽 = calculate_reinitial(𝚽, Δt)
        iter += 1
    end

    return 𝚽
end

raw"""
	classical_chan_vese(
		img; μ=0.25, λ₁=1.0, λ₂=1.0, tol=1e-3, max_iter=500, Δt=0.5, reinitial_flag=false
	)

Chan-Vese segmentation algorithm. Active contour model by evolving a level set. 
Can be used to segment objects without clearly defined boundaries.

# Arguments
- img: input image
- μ: "edge length" weight parameter. Higher `μ` values will
	produce a 'round' edge, while values closer to zero will
	detect smaller objects.
- λ₁: "difference from average" weight parameter for the output
	region with value 'true'. If it is lower than `λ₂`, this
	region will have a larger range of values than the other.
- λ₂: "difference from average" weight parameter for the output
	region with value 'false'. If it is lower than `λ₁`, this
	region will have a larger range of values than the other.
- tol: Level set variation tolerance between iterations. If the
	L2 norm difference between the level sets of successive
	iterations normalized by the area of the image is below this
	value, the algorithm will assume that the solution was
	reached.
- max_iter: Maximum number of iterations allowed before the algorithm
	interrupts itself.
- Δt: A multiplication factor applied at calculations for each step,
	serves to accelerate the algorithm. While higher values may
	speed up the algorithm, they may also lead to convergence
	problems.
- reinitial_flag: not sure what this does yet

# Citation
'Active contours without edges' [Chan, Vese](10.1109/83.902291)

'Chan–Vese Segmentation' [Pascal Getreuer](https://www.ipol.im/pub/art/2012/g-cv/article.pdf)

'The Chan-Vese Algorithm' [Rami Cohen](https://arxiv.org/pdf/1107.2782.pdf)
"""
function classical_chan_vese(
    img; μ=0.25, λ₁=1.0, λ₂=1.0, tol=1e-3, max_iter=500, Δt=0.5, reinitial_flag=false
)
    iter = 0
    D = ndims(img)
    if D == 3
        img = PermutedDimsArray(img, (2, 3, 1))
    end
    𝚽ⁿ = checkerboard(size(img))
    δ = tol + 1
    img .= img .- minimum(img)

    if maximum(img) != 0
        img = img ./ maximum(img)
    end

    while (δ > tol) & (iter < max_iter)
        𝚽ⁿ⁺¹ = calculate_variation(img, 𝚽ⁿ, μ, λ₁, λ₂, Δt)
        δ = sqrt(meanfinite((𝚽ⁿ⁺¹ .- 𝚽ⁿ) .^ 2, (1, 2))[1])
        if reinitial_flag
            𝚽ⁿ .= reinitialize(𝚽ⁿ⁺¹, Δt)
        else
            if D == 3
                r = axes(𝚽ⁿ⁺¹)
                𝚽ⁿ .= 𝚽ⁿ⁺¹[:, :, first(r[3])]
            else
                𝚽ⁿ .= 𝚽ⁿ⁺¹
            end
        end

        iter += 1
    end

    return 𝚽ⁿ, iter
end
