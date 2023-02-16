<h1> 
Technical Reference 
</h1>

## Random Dynamic System

The steady state response of a Linear Time-Invariant (LTI) structural dynamic system can be found by solving the discretized equation of motion in the frequency domain. 
The equation of motion in time and frequency domain are given below. 

<br>

$$
\mathbf{M} \ddot{\mathbf{u}} + 
\mathbf{C} \dot{\mathbf{u}} + 
\mathbf{K} \mathbf{u} =
\mathbf{f} 
\phantom{x}
\xrightarrow{
    \phantom{x}
    \mathscr{F}
    \phantom{x}
}
\phantom{x}
\left(
    \mathbf{K} - \omega^{2}\mathbf{M} + 
    i\omega \mathbf{C}
 \right)
 \mathbf{U} = 
 \mathbf{F}
$$

<br>

However, the equation depends on structural and load parameters, which are subject to uncertainties.
To account for this, the problem can be stated in a probabilistic settings, and the equation is reformulated as a function of random variables $\mathbf{x}$:

<br>

$$
\left[
    \mathbf{K}(\mathbf{x}) - \omega^{2}\mathbf{M}(\mathbf{x}) + 
    i\omega \mathbf{C}(\mathbf{x})
\right]
\mathbf{U}(\mathbf{x}) = 
\mathbf{F}(\mathbf{x})
$$

<br>

One can obtain distribution of $\mathbf{U}(\mathbf{x})$ using Monte Carlo Simulations (MCS) i.e. to solve the equations multiple times, each time using a different sample of $\mathbf{x}$. 
For large and complex system, the computational cost is too high, so one usually create a meta/surrogate model that can approximate the response while being cheaper to evaluate. 

<br>

## Polynomial Chaos Expansion (PCE)

In this project, two approaches are used; PCE and RPCE. These approximate response as a linear combination of basis functions and a ratio between two linear combinations of basis functions respectively:

<br>

$$
\mathbf{U}(\mathbf{x})
\approx 
\sum_{p=0}^{P-1} \Psi_{p}(\mathbf{x}) \mathbf{Y}_{p}
\phantom{x}
\text{or}
\phantom{x}
{U}_{i}(\mathbf{x})
\approx 
\frac{
    \sum_{p=0}^{P-1} \Psi_{p}(\mathbf{x}) {Y}_{pi}
}{
    \sum_{q=0}^{Q-1} \Psi_{q}(\mathbf{x}) {Y}_{qi}
}
$$

<br>

The basis functions are products of polynomial functions. The polynomials are chosen such that the expected value of product of two different polynomials is zero:

<br>

$$
\Psi_{p}(\mathbf{x}) =
\prod_{\alpha=0}^{d-1} \psi_{p_{\alpha}} \left(x_{\alpha}\right)
\phantom{x}
\text{where}
\phantom{x}
\int_{x\in\Omega_{x}}
\psi_{i}(x) \psi_{j}(x)
\cdot 
\mathscr{P}(x) \phantom{.} dx
= 0
\phantom{x}
\forall 
\phantom{x}
i\neq j
$$

<br>

For independent standard normal random variable, the polynomials are given by probabilist Hermite polynomials.
In this project, the polynomials are normalized:

<br>

$$
\int_{x\in\Omega_{x}}
\psi_{i}(x) \psi_{i}(x)
\cdot 
\mathscr{P}(x) \phantom{.} dx
= 1
\phantom{x}
\forall 
\phantom{x}
i\in \mathbb{Z}^{+}
$$

<br>

## Available Surrogate Models 

The two approaches are further classified into intrusive and non-intrusive. The former requires modification to the equations solver while the later is based on collocation method and can be paired with black-box solver. 
Further information on each model can be found in the following links:

<ul>
<li><a href="./IntrusivePCE.md">Intrusive PCE </a></li>
<li> Intrusive RPCE </li>
<li> Non - Intrusive PCE </li>
<li> Non - Intrusive RPCE </li>
</ul>