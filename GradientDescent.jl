### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ dfaa2c0e-882f-11eb-2d6b-715784168557
md"""
Ordinary Least Squares (OLS) problem:

$\arg \min _{\mathrm{{\boldsymbol w}}} \sum _{i=1}^ N \left(Y_ i - \mathrm{{\boldsymbol x}}_ i^{\intercal }\mathrm{{\boldsymbol w}}\right)^2,$
where $\mathrm{{\boldsymbol x}}_ i^{\intercal }\mathrm{{\boldsymbol w}}$ is the $\hat{Y_i}$. 
"""

# ╔═╡ f5019df6-8831-11eb-0a74-f7edd872d977
md"""
Function $f$ is convex if at each point, the gradient gives a linear lower bound, i.e., for all $u, w$:

```math
f(u) \geq f(w) + \langle \nabla f(w), u-w\rangle.
```
When $\nabla f(w) = 0$, the inner product is zero and for every point of $u$ and $w$, we also have

$f(u) \geq f(w),$
which means that the function has the global minimum value at $w$.

Chord across bowl: f is convex if for all $u, w$ and $0 \leq \lambda \leq 1$:

```math
f\left(\lambda w + (1-\lambda)u\right) \leq \lambda f(w) + (1-\lambda)f(u)
```

Jensen’s inequality: if $g$ is convex, then
```math
g(\mathbb{E}[X]) \leq \mathbb{E}[g(X)]
```
"""

# ╔═╡ 02d59a60-8835-11eb-08d1-83c437986d9c


# ╔═╡ Cell order:
# ╠═dfaa2c0e-882f-11eb-2d6b-715784168557
# ╠═f5019df6-8831-11eb-0a74-f7edd872d977
# ╠═02d59a60-8835-11eb-08d1-83c437986d9c
