<h1>
    Intrusive PCE
</h1>

## Approximation to The Equation of Motion

The matrices and vectors are expanded in terms of basis functions, where each basis function is a product of several Hermite polynomials.
The expansion is given by the following expression. 

<br>

$$
\mathbf{Z}(\mathbf{x}) =
\sum_{p=0}^{P-1} \Psi_{p}(\mathbf{x}) \cdot \mathbf{Z}_{p}
$$

<br>

Applying the expansion to mass, damping, stiffness matrices, as well as force and response vector, an approximation to the equation of motion is obtained:

<br>

$$
\sum_{k=0}^{P-1} 
\sum_{q=0}^{P-1} 
\Psi_{k}(\mathbf{x}) 
\Psi_{q}(\mathbf{x}) 
\cdot 
\left[
    \mathbf{K}_{k} - \omega^{2}\mathbf{M}_{k} +
    i\omega \mathbf{C}_{k}
\right]
\mathbf{U}_{q}
\approx
\sum_{q=0}^{P-1} \Psi_{q}(\mathbf{x}) \mathbf{F}_{q}
$$

<br>

Mass, damping, stiffness matrix, and load vector come from the problem formulation itself.
Therefore, the main task is to obtain the "best" $\mathbf{U}_{q}$ to fulfil the approximation.

<br>

## Weak Formulation 

Since solution to the strong formulation is not known (may not even exist), the weak formulation is used instead.
Using the orthogonal polynomial function and expectation, the strong form is converted to the weak form as in the following:

<br>

$$
a(\mathbf{x}) \approx b(\mathbf{x})
\phantom{x}
\longrightarrow
\phantom{x}
E\left[a(\mathbf{x})\cdot \Psi_{p}(\mathbf{x})\right]
=
E\left[b(\mathbf{x})\cdot \Psi_{p}(\mathbf{x})\right]
\phantom{x}
\forall
\phantom{x}
p \in [0,1,...,P-1]
$$

<br>

Applying the same conversion to the approximation of equation of motion, a new system of equations are obtained:

<br>

$$
\sum_{k=0}^{P-1} 
\sum_{q=0}^{P-1} 
E\left[
\Psi_{k}(\mathbf{x}) 
\Psi_{p}(\mathbf{x}) 
\Psi_{q}(\mathbf{x}) 
\right]
\cdot 
\left[
    \mathbf{K}_{k} - \omega^{2}\mathbf{M}_{k} +
    i\omega \mathbf{C}_{k}
\right]
\mathbf{U}_{q}
=
\sum_{q=0}^{P-1} 
E\left[
\Psi_{p}(\mathbf{x}) 
\Psi_{q}(\mathbf{x}) 
\right]
\mathbf{F}_{q}
$$

<br>

Taking advantage of polynomials' orthonormality, the equation is given in matrix notation below:

<br>

$$
\begin{aligned}
\bar{\mathbf{F}}
=
\bar{\mathbf{H}} 
\bar{\mathbf{U}},
\phantom{xxxxx}
\bar{\mathbf{H}} &= 
\sum_{k=0}^{P-1} 
\mathbf{A}_{k} 
\otimes 
\left[
    \mathbf{K}_{k} - \omega^{2}\mathbf{M}_{k} +
    i\omega \mathbf{C}_{k}
\right]
\\
\bar{\mathbf{U}} &=
\left[
    \mathbf{U}_{0}^{T},
    \mathbf{U}_{1}^{T}, 
    ..., 
    \mathbf{U}_{P-1}^{T}, 
\right]^{T}
\\
\bar{\mathbf{F}} &=
\left[
    \mathbf{F}_{0}^{T},
    \mathbf{F}_{1}^{T}, 
    ..., 
    \mathbf{F}_{P-1}^{T}, 
\right]^{T}
\end{aligned}
$$

<br>

The term $\mathbf{A}_{k}$ is given by the following:

<br>

$$
\begin{aligned}
\left[\mathbf{A}_{k}\right]_{pq} &=
E\left[
\Psi_{k}(\mathbf{x}) 
\Psi_{p}(\mathbf{x}) 
\Psi_{q}(\mathbf{x}) 
\right] 
\\
&=
\prod_{\alpha=0}^{d-1}
E\left[
    \psi_{k_{\alpha}} (x_{\alpha})
    \psi_{p_{\alpha}} (x_{\alpha})
    \psi_{q_{\alpha}} (x_{\alpha})
\right]
\end{aligned}
$$

<br>

If $s=0.5(k_{\alpha}+p_{\alpha}+q_{\alpha})$ is an integer greater than $\max{(k_{\alpha},p_{\alpha},q_{\alpha})}$, the expectation term evaluates to the following expression. Otherwise, it evaluates to zero.

<br>

$$
E\left[
    \psi_{k_{\alpha}} (x_{\alpha})
    \psi_{p_{\alpha}} (x_{\alpha})
    \psi_{q_{\alpha}} (x_{\alpha})
\right]
= 
\frac{
    \sqrt{ k_{\alpha}! p_{\alpha}! q_{\alpha}!}
}{
    (s-k_{\alpha})! (s-p_{\alpha})! (s-q_{\alpha})!
}
$$
