### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 7fecf5ea-869e-11eb-1b61-bfec1817f730
begin
	using StatsPlots
	X = [0.0339, 0.0423, 0.213, 0.257, 0.273, 0.273, 0.450, 0.503, 0.503,
		0.637, 0.805, 0.904, 0.904, 0.910, 0.910, 1.02, 1.11, 1.11, 1.41,
		1.72, 2.03, 2.02, 2.02, 2.02]
	Y = [-19.3, 30.4, 38.7, 5.52, -33.1, -77.3, 398.0, 406.0, 436.0, 320.0, 373.0,
		93.9, 210.0, 423.0, 594.0, 829.0, 718.0, 561.0, 608.0, 1.04E3, 1.10E3,
		840.0, 801.0, 519.0]
	N = 24
	plot(X, Y, seriestype = :scatter, title = "Hubble", xlabel = "Distance (Mpc)", ylabel = "Apparent velocity(km/s)", label = "galaxies")
end

# ╔═╡ 6eed7c20-869d-11eb-3452-2ded716d82ae
md"""
# 1. Revise
### Mean
$\bar{X} = \frac{1}{N} \sum _{i=1}^ N X_ i$
### Standard Deviation
$s_ x = \sqrt{N^{-1} \sum _{i=1}^ N (x_ i - \bar{x})^2}$
### Covariance
Sample covariance: 

$s^2_{X,Y} = \frac{1}{N-1} \sum _{i=1}^ N (X_ i - \bar{X}) (Y_ i - \bar{Y})$
The extreme values:
```math
\begin{split}
|\text{Cov}(X,Y)|^2 &= |\mathbb{E}((X-\mu)(Y-\nu))|^2 \\
& = |\langle X-\mu, Y-\nu \rangle|^2 \\
& \leq \langle X-\mu, X-\mu\rangle \langle Y-\nu, Y-\nu\rangle \\
& = \mathbb{E}((X-\mu)^2)\mathbb{E}((Y-\nu)^2) \\
& = \text{Var}(X)\text{Var}(Y) \sim s_X s_Y
\end{split}
```
### Correlation coefficient
The sample correlation coefficient:

$r_{X,Y} = \frac{s^2_{X,Y}}{s_ X s_ Y} = \frac{1}{N-1} \sum _{i=1}^ N \left(\frac{X_ i - \bar{X}}{s_ X}\right) \left(\frac{Y_ i - \bar{Y}}{s_ Y}\right)$
The population correlation coeeficient:

$\rho_{X,Y} =\frac{\text{Cov}(X,Y)}{\sigma_X\sigma_Y} =\frac{\mathbb{E}[(X-\mu_X)(Y-\mu_Y)]}{\sigma_X\sigma_Y}$
### Data Set
The $Y$-axis is the apparent velocity, measured in kilometers per second. Positive velocities are galaxies moving away from us, negative velocities are galaxies that are moving towards us. The $X$-axis is the distance of the galaxy from us, measured in mega-parsecs (Mpc); one parsec is 3.26 light-years, or 30.9 trillion kilometers.
"""

# ╔═╡ Cell order:
# ╟─6eed7c20-869d-11eb-3452-2ded716d82ae
# ╟─7fecf5ea-869e-11eb-1b61-bfec1817f730
