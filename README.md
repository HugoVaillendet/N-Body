# N-Body

## Problem presentation

N-body problems are a class of gravitationnal problems where the number of bodies $N \geq 3$.
In these problems we consider only the gravitaionnal force applied by a body $j$ on a body $i$ :

$$\vec{F}_{i} = \sum_{i=1}^{N}G \frac{m_i m_j}{\lVert \vec{r}_{ij} \rVert^3}\vec{r}_{ij}$$

with m<sub>i</sub>, m<sub>j</sub> the respective masses, ||r<sub>ij</sub>|| the distance between i and j, r<sub>ij</sub> the distance vector and G the gravitational constant.

N-body problems difficulty lie in the absence of general analytical solutions for $N \geq 3$. To circumvent this we will use numerical methods.

## Numerical integration method

We need first to find a suitable integration method. One issue with common integrators like Euler or RK4 is that they are not symplectic (meaning they don't conserve energy over time).
We will instead use the leapfrog integrator (a variant of the verlet integration method). The algorithm works as follows :
