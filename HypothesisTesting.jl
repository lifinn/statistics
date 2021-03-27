### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 287e7c2e-768c-11eb-0916-13431a504633
begin
	using Interact, InteractiveUtils, PlutoUI
	using Distributions
	using StatsPlots
	H1 = Binomial(31000, 0.0013)
	H0 = Binomial(31000, 0.002)
	@bind α Slider(0.01:0.01:0.1, show_value=true, default=0.05)
end

# ╔═╡ 8db0b01e-7b76-11eb-0606-15cc778e4879
begin
	using HypothesisTests
	drug =    [6.1, 7.0, 8.2, 7.6, 6.5, 7.8, 6.9, 6.7, 7.4, 5.8]
	placebo = [5.2, 7.9, 3.9, 4.7, 5.3, 4.8, 4.2, 6.1, 3.8, 6.3]
	X = drug - placebo
	OneSampleTTest(X)
end

# ╔═╡ 66ef8d9c-7617-11eb-0204-7fcfc3745ca3
md"""
# Hypothesis Testing
In hypothesis testing, we start with some default theory -- called a null hypothesis -- and we ask if the data provide sufficient evidence to reject the theory. If not we retain the null hypothesis. However, a more accurate term should be "failing to reject the null", because we are in the process looking for evidence against $H_0$ in the data and we decide whether to reject $H_0$, rather than looking for evidence that in favour of $H_0$. Thus, we can never prove $H_0$ is true, because we do not prove it from the beginning. The real starting point is that we assume $H_0$ is true and try to find anything to prove its faultiness. Even if we do not find any evidence against $H_0$, this does not mean that $H_0$ is true.

In this setting, the $H_0$ and $H_1$ are not in the sysmetric role and the data collected is only used for disapprove $H_0$. The process of selecting $H_0$ and $H_1$, therefore, has a critical effect on the scientific study. In the court, the verdict of jury is "not guilty" instead of "innocent".
"""

# ╔═╡ 4aba30da-7696-11eb-3786-1b7e4dd442a4
md"""
# 1. Test statistic
The decision whether to reject the null hypothesis is based on a test statistic. The test statistic is a function of the random variables modelling the data. Hence it is a random variable itself, and its distribution depends on the parameters defining the model. For any specific hypothesis test, the test statistic that we choose needs to distinguish between the null and the alternative hypotheses, and have a distribution that is known and computable.
"""

# ╔═╡ 3c4116a6-7617-11eb-3e16-998daaae41ba
md"""
# 2. Errors
- A test is a statistic $\psi \in \{0, 1\}$ such that:

1. if $\psi = 0$, $H_0$ is not rejected
2. if $\psi = 1$, $H_0$ is rejected

- Reject region of the test $\psi$:

$R_{\psi} = \{x \in E^n: \psi(x) = 1\},$

which is the subset of samples that lead to a rejection.

- Type I error of a test $\psi$ (rejecting $H_0$ when it is actually true)

Function $\alpha_{\psi}$ is from $\Theta_0$ to $\mathbb{R}$ that maps $\theta$ in $\Theta_0$ to a probability $\mathbb{P}_{\theta}[\psi = 1]$. Thus, $\mathbb{R} \in [0, 1]$.

- Type II error of a test $\psi$ (not rejecting $H_0$ when $H_1$ is actually true)

Function $\beta_{\psi}$ is from $\Theta_1$ to $\mathbb{R}$ that maps $\theta$ in $\Theta_1$ to a probability $\mathbb{P}_{\theta}[\psi = 0]$. Thus, $\mathbb{R} \in [0, 1]$.

- The power of a test is defined as the probability of rejecting $H_0$ when $H_1$ is true

$\pi_{\psi} = \inf_{\theta \in \Theta_1}(1 - \mathbb{P}(\beta_{\psi})) \in [0, 1]$

The worst possible value you could take over all possible realities is H_1.
"""

# ╔═╡ 14af85b8-769a-11eb-0a79-fb682c101e0e
md"""
- Null Hypothesis $H_0$: $\pi=0.002$  (blue distribution on the right)
- Alternative Hypothesis $H_1$: $\pi=0.0013$  (orange distribution on the left

The distribution of the test statistic under the null and the alternative hypotheses are shown in the plot below. The blue distribution on the right corresponds to the null hypothesis $H_0$, while the orange distribution on the left corresponds to the alternative hypothesis $H_1$.
We reject $H_0$ if our sample has a $p$-value less than $\alpha=0.05$, i.e. any value on the $x$-axis to the left of the vertical line separating the green from the pink regions. The power is then the area to the left of this line that is under the curve given $H_1$(the left curve).
"""

# ╔═╡ 8d9b8ad6-76ad-11eb-1368-19bac4d4e824
md"""
## 3. The power of a test
"""

# ╔═╡ 95d86d8a-76b3-11eb-3d95-0f227c556cd6
md"""
The test $\psi$ has level $α$ that

$α_{\psi}(\theta) \leq α, \forall\theta \in \Theta_0$
"""

# ╔═╡ 9c4b9540-76b1-11eb-07a4-4134a1d33132
md"""
The power of the test for $H_0$ Binomial(31000, 0.002) and $H_1$ Binomial(31000, 0.0013)"""

# ╔═╡ 4b60dc5e-7688-11eb-3e70-5f19fca25155
begin
	plot(0:1:100, pdf.(H0, 0:1:100), label="\$H_0\$")
	plot!(0:1:100, pdf.(H1, 0:1:100), label="\$H_1\$")
	plot!([quantile(H0, 1E-10):0.01:quantile(H0, α)], seriestype="vspan", label="Power", xticks=([quantile(H0, α)], [quantile(H0, α)]))
end

# ╔═╡ accdb8ca-76ad-11eb-3061-ddc6c68ad5ee
begin
	plot(cdf.(H1, 1:1:80), label="\$H_1\$")
	plot!([quantile(H0, α)], seriestype="vline", xticks=([quantile(H0, α)], [quantile(H0, α)]), label="Critical")
	scatter!([quantile(H0, α)], [cdf(H1, quantile(H0, α))], label="Power")
end


# ╔═╡ 6ee2a066-76ab-11eb-3faf-258b06fc9e0c
md"""
The type I error is entirely controllable because we decide on the boundary for the test statistic as to whether we will reject the null hypothesis.

In a simple hypothesis test setting, the significance level $\alpha$ and the test statistic together define the rejection region. We can then compute the probability of a type II error as the probability that the test statistic takes a value outside the rejection region, under the alternative hypothesis.

If we use a different test statistic, then its distribution under the null hypothesis and the alternative hypothesis can also change. Hence, even under the same significance level, the type II error can be different.

There is a direct tradeoff between reducing the type I and type II errors. In nearly all cases, with a fixed test statistic, reducing the type I error results into increasing the type II error, and vice versa.

Keeping the significance level constant, a one-sided hypothesis test has a higher power than the corresponding two-sided hypothesis test. An exception is when the distribution of the observation under the alternative hypothesis is bimodal.
"""

# ╔═╡ 0651542e-7b7c-11eb-2c2c-336b9eeedec2
md"
# 4. T-test
- T-statistic:

$T_n:=\frac{\bar{X}_n -\mu}{\sqrt{\hat{\sigma^2}/n}},$

where

$\bar{X}_n:=\frac{1}{n}\sum_{i=1}^n X_i,$

$\sigma^2 = \frac{1}{n-1}\sum_{i=1}^n (X_i - \bar{X}_n)^2$ and

$X_i, \dots, X_n \stackrel{i.i.d}{\sim} \mathcal{N}(0,1).$
-  $\chi^2$ distribution with $n$ degree of freedom is defined as

$\sum_{i=1}^n y_i^2 \sim \chi^2_n \sim \sum_{i=1}^n z_i,$ where

$y_i \sim \mathcal{N}(0,1)$
$z_i \sim \chi^2_1$

- t distribution with $n$ degree of freedom

$\frac{y}{\sqrt{z/n}} \sim t_n,$

where

$y \sim \mathcal{N}(0,1)$
$z \sim \chi^2_n$
"

# ╔═╡ 555c2b12-7b85-11eb-1f26-6b972a853c9d
md"
## 4.1 We claim $T_n \sim t_{n-1}$

By the definition of T-statistic,

$T_n:=\frac{\bar{X}_n -\mu}{\hat{\sigma}/\sqrt{n}},$

Replace $\hat{\sigma}$ with $\sigma$,

$T_n:=\frac{\frac{\bar{X}_n -\mu}{\sigma/\sqrt{n}}} {{\sqrt{\hat{\sigma}^2/\sigma^2}}},$

where the nominator $\sim \mathcal{N}(0,1)$, and the denominator:

$\sqrt{\frac{1}{n-1}}\sqrt{\frac{1}{\sigma^2}\sum_{i=1}^n(X_i - \bar{X}_n)^2}.$

By the definition of t distribution with $n$ degree of freedom, we have $y$ and $n-1$. The rest of work is to prove that the $\frac{1}{\sigma^2}\sum_{i=1}^n(X_i - \bar{X}_n)^2$ is $\chi_{n-1}^2$.

$\chi^2_n = \sum_{i=1}^n (\frac{x_i - \mu}{\sigma})^2 = \frac{1}{\sigma^2}\sum_{i=1}^n(X_i - \bar{X}_n + \bar{X}_n - \mu)^2$

Collect terms and then we have,

$\frac{1}{\sigma^2}\sum_{i=1}^n(X_i-\bar{X}_n)^2 + \frac{n}{\sigma^2}(X_n-\mu)^2 + 2\frac{1}{\sigma^2}(\bar{X}_n-\mu)\sum_{i=1}^n(X_i-\bar{X}_n),$

where $\sum_{i=1}^n(X_i-\bar{X}_n) = 0.$

$\frac{n}{\sigma^2}(X_n-\mu)^2 = (\frac{\bar{X}_n - \mu}{\sigma/\sqrt{n}})^2 \sim \chi_1^2,$

thus 

$\chi_n^2 = \frac{1}{\sigma^2}\sum_{i=1}^n(X_i-\bar{X}_n)^2 + \chi_1^2$


$\frac{1}{\sigma^2}\sum_{i=1}^n(X_i-\bar{X}_n)^2 = \chi_{n-1}^2$

In summary T-statistic follows a t-distribution with $n-1$ degrees of freedom.

"

# ╔═╡ 97d7fdb2-7f55-11eb-2529-3daac74edaf6
md"""## 4.2 Case Study
We have two datasets
```julia
drug =    [6.1, 7.0, 8.2, 7.6, 6.5, 7.8, 6.9, 6.7, 7.4, 5.8]
placebo = [5.2, 7.9, 3.9, 4.7, 5.3, 4.8, 4.2, 6.1, 3.8, 6.3]
```
 $X$ is the difference between the sleep time of patients who takes *drug* and that of those who take *placebo*. The size of $X$ is 10. 
### 1. State the null and alternative hypotheses.
Let $\mu$ denote the mean value of $X$. The null hypothesis is that the sleep time of both groups of patients have no difference, that is, the mean equals to 0, while the alternative hypothesis is that the drug extends the sleep time of patients, that is, the mean is greater than 0.

$H_0: \mu = 0$
$H_1: \mu > 0$

### 2. Discuss the logic of this hypothesis test.
Basically, the logic of this hypothesis test is as follows: If the null hypothesis
is true, then the mean difference, $X$, of the sample of patients should approximately equal 0. We say “approximately equal” because we cannot
expect a sample mean to equal exactly the population mean; some sampling error
is anticipated. However, if the sample mean difference of sleep time is “too much
greater” than 0, we would be inclined to reject the null hypothesis in
favor of the alternative hypothesis.

### 3. Obtain a precise criterion for deciding whether to reject the null hypothesis in favor of the alternative hypothesis.
$T_n:=\frac{\bar{X}_n -\mu}{\hat{\sigma}/\sqrt{n}} = \frac{1.78 - 0}{1.7681/\sqrt{10}} = 3.1835$
"""


# ╔═╡ b8f55364-7b77-11eb-2d65-d5b556b5af4b
begin
	plot(-4:0.1:4, pdf.(TDist(9), -4:0.1:4), label="\$t(9)\$")
	plot!(-4:0.1:4, pdf.(Normal(0, 1), -4:0.1:4), label="\$N(0,1)\$")
	plot!([3.18], seriestype="vline", xticks=([3.18], [3.18]), label="Critical")
end

# ╔═╡ ddd94a7a-7f58-11eb-1fb0-f9657b31d15e
md"""
### 4. Apply the criterion to the sample data and state the conclusion.

We see that the value of the test statistic, −2.65, is less than the cutoff point of −1.645 and, hence, we reject $H_0$. We calculated the t-statistic, 3.1835, is greater than the cutoff point of the critical value of 95% level, $(round(quantile(TDist(9), 0.95), digits=4)), and, hence, we reject the null hypothesis. 
"""

# ╔═╡ Cell order:
# ╟─66ef8d9c-7617-11eb-0204-7fcfc3745ca3
# ╟─4aba30da-7696-11eb-3786-1b7e4dd442a4
# ╟─3c4116a6-7617-11eb-3e16-998daaae41ba
# ╟─14af85b8-769a-11eb-0a79-fb682c101e0e
# ╟─8d9b8ad6-76ad-11eb-1368-19bac4d4e824
# ╟─95d86d8a-76b3-11eb-3d95-0f227c556cd6
# ╟─287e7c2e-768c-11eb-0916-13431a504633
# ╟─9c4b9540-76b1-11eb-07a4-4134a1d33132
# ╟─4b60dc5e-7688-11eb-3e70-5f19fca25155
# ╟─accdb8ca-76ad-11eb-3061-ddc6c68ad5ee
# ╟─6ee2a066-76ab-11eb-3faf-258b06fc9e0c
# ╟─0651542e-7b7c-11eb-2c2c-336b9eeedec2
# ╟─555c2b12-7b85-11eb-1f26-6b972a853c9d
# ╟─97d7fdb2-7f55-11eb-2529-3daac74edaf6
# ╠═8db0b01e-7b76-11eb-0606-15cc778e4879
# ╟─b8f55364-7b77-11eb-2d65-d5b556b5af4b
# ╠═ddd94a7a-7f58-11eb-1fb0-f9657b31d15e
