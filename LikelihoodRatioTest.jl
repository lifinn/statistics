### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 065fc4aa-834c-11eb-31bd-df1b0297bb3b
begin
	using Distributions, StatsPlots
	plot(0:0.05:8, pdf.(Chisq(3), 0:0.05:8), label = "\$ x^2_3\$")
	plot!([0.48], seriestype="vline", xticks=([0.48], [0.48]), label="Critical")
	plot!([quantile(Chisq(3), 0.95)], seriestype="vline", xticks=([quantile(Chisq(3), 0.95)], [quantile(Chisq(3), 0.95)]), label="95%")

end

# ╔═╡ 0997866e-82af-11eb-025c-c97f848f2a99
md"""
# Likelihood Ratio Test
## 1. Likelihood Ratio Test Definitions
Consider testing

$H_0: \theta \in \Theta _0, \, H_ A: \theta \in \Theta _ A.$ where $\Theta _0$ and $\Theta _A$ are disjoint subsets $(\Theta _0 \cap \Theta _ A =\emptyset)$  of a parameter space $\Theta =\Theta _0\cup \Theta _ A.$

The likelihood ratio test statistic $\Lambda (x)$ is defined as (negative twice) the logarithm of the likelihood ratio $L(x)$

$\Lambda (x) = -2\log (L(x)) \qquad \text {where } L(x)\, =\, \frac{\max _{\theta \in {\color{blue}{\Theta _0}} }p(x;\theta )}{\max _{\theta \in {\color{blue}{\Theta }} }p(x;\theta )}$

Equivalently, in the language of maximum likelihood estimators,

$L(x) = \frac{ p\left(x;\hat{\theta }_{\text {MLE}}^{\text {constrained}}\right)}{p\left(x;\hat{\theta }_{\text {MLE}}\right)}$

where $\hat{\theta }_{\text {MLE}}$ is the maximum likelihood estimator of $\theta$ and $\hat{\theta }_{\text {MLE}}^{\text {constrained}}$ is the constrained maximum likelihood estimator of $\theta$ within $\theta_0$.) In practice, replacing $\Theta _A$ with $\Theta$ has little effect on the test statistic. Moreover, the theoretical properties of $\lambda$ are much simpler in this way.

> Given the definition of the likelihood ratio test, if we expect the null hypothesis to be false, $L(x) \lt \lt 1$ and $\Lambda(x) \rightarrow \infty$.
"""

# ╔═╡ ea031e44-82b3-11eb-382f-31eff90163ae
md"""
## 2. The Distribution of Likelihood Ratio test statistics
Wilk's Theorem states that when the sample size is large, the distribution of $\Lambda$ under $H _0$ approaches a $\chi^2$ distribution:

$\Lambda \overset {n\to \infty }{\longrightarrow } \chi _ d^2 \qquad \text {where } d = \text {dim}(\Theta ) - \text {dim}(\Theta _0),$

where $d$ is the degree of freedom of the $\chi ^2$ distribution.
"""

# ╔═╡ 944d9212-82b4-11eb-2aa4-dfae6b54fb56
md"""
## 3. Power of Likelihood Ratio Test
The Neyman–Pearson lemma states that among all tests that test for the simple hypotheses $H_0: \theta =\theta _0\, ;\, H_ A:\theta =\theta _ A$ at significance level $\alpha$, the likelihood ratio test is the most powerful. That is, among all tests testing the same simple hypotheses and at the same significance level, the likelihood ratio test gives the largest probability of rejecting the null when indeed the alternate is true.
"""

# ╔═╡ 94b1f7fe-82b7-11eb-1cc8-45312d47854b
md"""
## 4. A case study
There are 4 types of progeny: round yellow, wrinkled yellow, round green, wrinkled green. The number of each type is multinomial with probability $p = (p_1, p_2, p_3, p_4)$. We guess that $p$ is equal to

$p_0 = (\frac{9}{16},\frac{3}{16},\frac{3}{16},\frac{1}{16})$

In $n=556$ trails he observed $X = (315, 101, 108, 32)$. We want to test 

$H_0: p = p_0 \text{ vs } H_1: p \neq p_0$

### 4.1 Log-Likelihood of multinomial distribution

Multinomial Distribution is a multivariate analog of the binomial distribution. Real vector p in the parameter space

${p \in \mathbb{R}^k: 0 \leq p_i, i = 1, \dots, k, \text{ and } \sum_{i=1}^k p_i = 1}$

Its probability Mass Function

$f(x;p) = \frac{n!}{\prod_{i=1}^k x_i!}\prod_{i=1}^k p_i^{x_i}, x \in S,$

where sample space

$S = {x \in \mathbb{Z}^k: 0 \leq x_i, i = 1, \dots, k, \text{ and } \sum_{i=1}^k x_i = n}$

The likelihood of $x_i:$

$L(x;p) = \prod_{i=1}^k p_i^{x_i}$
"""

# ╔═╡ 9ee8d046-8349-11eb-21c8-eba0d9f37a23
begin
	p₀ = [9/16, 3/16, 3/16, 1/16]
	X = [315, 101, 108, 32]
	n = 556
	p = X ./ n
	log_values = log.(p₀ ./ p)
	-2 * sum(X .* log_values)
end

# ╔═╡ 003a40a0-8347-11eb-045c-4df433a0b35b
md"""
### 4.2 Likelihood Ratio Test Statistic

From the definition we have

$\log{\frac{L(x;p_0)}{L(x;p)}} = \log{\frac{\prod_{i=1}^k p_{0i}^{x_{i}}}{\prod_{i=1}^k p_i^{x_i}}} = \sum_{i=1}^k x_i \log{\frac{p_{0i}}{p_i}}.$

Substitute the values we have test statistic $(round(-2 * sum(X .* log_values), digits=2))
"""

# ╔═╡ 01d1b83e-8352-11eb-1f66-111a01409ee0
begin
	using HypothesisTests
	MultinomialLRTest(X, p₀)
end

# ╔═╡ edb0ce64-834a-11eb-100f-e1a1cdc76492
md"""
### 4.3 Compare with $\chi^2$ distribution
Parameter $p$ has 4 values, but the sum of all parameters should be 1, thus we lost one degree of freedom and it lives in dimention 3. Parameter $p_0$ has only point, which is dimention of 0. In summary, the $\chi^2$ distribution has degree of freedom 3.
"""

# ╔═╡ 749cbf7a-834e-11eb-3794-27be6f7ea1d6
md"""
Therefore, the limiting distribution of $p$ under $H_0$ is $\chi^2_3$ and the $p$-value is $\mathbb{P}(\chi^2_3 > 0.48)$ = $(round(1 - cdf(Chisq(3), 0.48), digits = 2))

It is not evidence against $H_0$.
"""

# ╔═╡ cfa1aa10-8352-11eb-01cf-13942fb3ed55
md"""
### 4.4 Using Julia to do $\chi^2$ ratio test
```julia
MultinomialLRTest(observed value :> Integers, p_0 :> Float64)
```
"""

# ╔═╡ 674692a6-8510-11eb-1726-cdc17553ac81
md"""
### 4.5 Binomial and Multinomial test
In terms of the breast cancer death rates, $\pi_y$ in the treatment group and $\pi_c$ in the control group, are

$\pi _ T=\pi _ C,$
$\pi _ T\neq\pi _ C$

Let $Y_t$ and $Y_c$ be the numbers of cancer deaths in the treatment and control groups respectively. Assuming these are independent from each other, the probability of having $y_t$ breast cancer deaths in the treatment group and $y_c$  breast cancer deaths in the control group is the product.

```math
\begin{bmatrix}
          & death & alive  & total  \\
treatment & 39    & 30,961 & 31,000 \\
control   & 63    & 30,937 & 31,000 \\
total     & 102   & 61,898 & 62,000
\end{bmatrix}
```
---
This question can be viewed as a multinomial problem that the observed values of $\texttt{td}$, $\texttt{ta}$, $\texttt{cd}$ and $\texttt{cd}$ are $[39, 30961, 63, 30937]$ respectively, while the probability values under $h_0$ are $[0.0008225, 0.499177, 0.0008225, 0.499178]$. However, we only have 1 degree of freedom for $\Theta$ in this problem, because the the death rate and ratio of total death can generate the rest data.
"""

# ╔═╡ 65ea43fe-850f-11eb-3827-9360d1f3026a
begin
	π₀ = 102/62000 / 2 
	rπ₀ = 61898/62000 / 2
	MultinomialLRTest([39,30961,63,30937], [0.0008225,0.499177,0.0008225,0.499178])
end

# ╔═╡ 1a84e07c-8514-11eb-045c-d52a5de159de
begin
	plot(0:0.05:8, pdf.(Chisq(1), 0:0.05:8), label = "\$ x^2_3\$")
	plot!([5.71], seriestype="vline", xticks=([5.71], [5.71]), label="Critical")
	plot!([quantile(Chisq(1), 0.95)], seriestype="vline", xticks=([quantile(Chisq(1), 0.95)], [quantile(Chisq(1), 0.95)]), label="95%")

end

# ╔═╡ 86781346-8585-11eb-077e-55197182c28b
md"""
# 5. Problems in multiple testing
As we conduct the hypothesis testing based on $\alpha$% level, if the testing number is large enough, it is almost sure that we have $\alpha$% false discoveries in total.

Consider m hypothesis tests:

$H_{0i} \text{ vs. } H_{1i}, i = 1, \dots, m$ and let $P_1, \dots, P_m$ denote the $m$ $p$-value for these tests.

Two following metrics for false significance:
- Family-wise error rate (FWER) : the probability of making at least one false discovery, or type I error;
- False discovery rate (FDR) : the expected fraction of false significance results among all false significance results.

## 5.1 FWER
$\text {FWER}= \mathbf{P}({\color{blue}{V}}  \geq 1),$ where V is the total number of type I error.

When we have FWER with no corrections:

$\text {FWER} = \mathbf{P}(V\geq 1) = 1-\mathbf{P}(V=0)= 1 -(1- \alpha )^ m \approx 1 \qquad \text {for large } m.$

### Bonferroni Correction
To control FWER $\lt \alpha$, we can “correct" the $p$-value $p_i$ of each test to $mp_i$ and reject $H^i_0$  when $mp_i\lt\alpha$, where $m$ is the number of tests. This correction is called the Bonferroni correction.

$p_i \lt \frac{\alpha}{m}$

### Example for the Bonferroni Method
If $\alpha$ = .05 and the number of tests, m, are 2,638, $p_i = 0.00001895$. Hence, for any test with $p$-value $\lt 0.00001895$, we declare that there is a significant difference.
"""

# ╔═╡ 74ff1d0e-85a7-11eb-27a4-9977aded8672
md"""
## 5.2 FDR
Suppose we reject all null hypotheses whose $p$-value fall below some threshold. Let $m_0$ be the number of null hypotheses that are true and let $m_1 = m - m_0$. The test can be categorised in a 2 * 2 table. The definition of the false discovery proportion
```math 
\text{FDP} = 
\begin{cases} 
      V/R & R\gt 0 \\
      0 & R = 0, 
\end{cases}
```
where $V$ is the number of false discoveries and $R$ is the total number of cases that reject $H_0$.

### The Holm-Bonferroni Correction
1. Sort the $m$ $p$-values in increasing order $p^{(1)}\leq p^{(2)}\leq \ldots \leq p^{(i)}\leq \ldots \leq p^{(m)}$
2. Start with $i=1$, i.e. the minimum $p$-value. If $p^{(i)}\leq \frac{\alpha }{m-(i-1)}$, reject $H_0^{(i)}$
3. As soon as a hypothesis, say $H_0^{(k)}$, is not rejected, stop the procedure and simply do not reject any of $H_0^{(k)}, H_0^{(k+1)}\cdots , H_0^{(m)}$

The Holm-Bonferroni method is more powerful than the Bonferroni correction, since it increases the chance of rejecting the null hypothesis (accepting discovery), and hence reduces the chance of making type II error.

### The Benjamini-Hochberg Correction
The Benjamini-Hochberg method guarantees $\text {FDR}<\alpha$ for a series of $m$ independent tests. The procedure is as follows:
1. Sort the $m$ $p$-values in increasing order: $p^{(1)}\leq p^{(2)}\leq \ldots \leq p^{(i)}\leq \ldots \leq p^{(m)}$
2. Find the maximum $k$ such that
$p^{(k)}\leq \frac{k}{m}\alpha.$
3. Reject all of $H_0^{(1)}, H_0^{(2)}\cdots , H_0^{(k)}$

### Example for Benjamini-Hochberg Correction Method
Suppose that 10 independent hypothesis tests are carried leading to the following ordered $p$-values: $[0.00017, 0.00448, 0.00671, 0.00907, 0.01220, 0.33626, 0.39341, 0.53882, 0.58125, 0.98617]$
"""


# ╔═╡ a32bfd0c-85aa-11eb-0758-17363093ec94
begin
	pvalues = [0.00017, 0.00448, 0.00671, 0.00907, 0.01220, 
		       0.33626, 0.39341, 0.53882, 0.58125, 0.98617]
	plot(1:10, pvalues, seriestype = :scatter, title = "Three multiple testing methods", label = "p-values")
	α = 0.05
	plot!([0.05/10], seriestype = :hline, label = "B")
	HB(i) = α / (10 - i + 1)
	plot!(1:10, HB.(1:10), label = "HB")
	BH(i) = i * α / 10
	plot!(1:10, BH.(1:10), label = "BH")
end


# ╔═╡ 2c456e62-85d2-11eb-3fab-45c318292605
md"""
## 5.3 FDR versus FWER
Compared to $\text {FWER}$, $\text {FDR}$ has higher power. Put another way, $\text {FWER}$ is stricter than $\text {FDR}$. $\text {FDR}$ is easier to control than $\text {FWER}$. That is, we do need to apply as large a correction factor for the  $\text {FDR}$ than we do to the $\text {FWER}$ to get it under our significance threshold $\alpha$. Therefore, the power of a series of tests with $\text {FDR}$ controlled will be larger that the power of the series with $\text {FWER}$ controlled.
- When generating hypotheses, it is best to report the original p-values, and apply no correction, so that the power of the test is high and we don't reject a hypothesis too early. We must also report the number of hypotheses we tested, so that other scientists can accurately judge if the hypothesis merits further testing. 
- When approving pharmaceuticals, we do not wish to tolerate even a single false positive. A false positive could result in a patient being given a treatment that does not work, when there are other drugs available that are effective. Hence, we want $\text {FWER} \leq 5\%$.
- When screening, we may be testing a large number of hypotheses with the goal of identifying just a few candidates that are promising for further study. If the number of hypotheses is very large, the $\text {FWER}$ correction factor will be extremely strict, as it controls for even one false positive. In comparison, as the $\text {FDR}$ is based on a ratio, and so it scales with the number of hypotheses performed. We will also get proportionately more false positives, but this is acceptable as such errors will be corrected by further studies (such as clinical trials). Hence, we want $\text {FDR} \leq 10\%$
"""

# ╔═╡ Cell order:
# ╟─0997866e-82af-11eb-025c-c97f848f2a99
# ╟─ea031e44-82b3-11eb-382f-31eff90163ae
# ╟─944d9212-82b4-11eb-2aa4-dfae6b54fb56
# ╟─94b1f7fe-82b7-11eb-1cc8-45312d47854b
# ╟─003a40a0-8347-11eb-045c-4df433a0b35b
# ╠═9ee8d046-8349-11eb-21c8-eba0d9f37a23
# ╟─edb0ce64-834a-11eb-100f-e1a1cdc76492
# ╟─065fc4aa-834c-11eb-31bd-df1b0297bb3b
# ╟─749cbf7a-834e-11eb-3794-27be6f7ea1d6
# ╟─cfa1aa10-8352-11eb-01cf-13942fb3ed55
# ╟─01d1b83e-8352-11eb-1f66-111a01409ee0
# ╟─674692a6-8510-11eb-1726-cdc17553ac81
# ╟─65ea43fe-850f-11eb-3827-9360d1f3026a
# ╠═1a84e07c-8514-11eb-045c-d52a5de159de
# ╟─86781346-8585-11eb-077e-55197182c28b
# ╟─74ff1d0e-85a7-11eb-27a4-9977aded8672
# ╟─a32bfd0c-85aa-11eb-0758-17363093ec94
# ╟─2c456e62-85d2-11eb-3fab-45c318292605
