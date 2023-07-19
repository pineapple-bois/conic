# `Conics`
### A Python class to classify, manipulate and visualise conic sections.

----

Conic types are instantiated from a bi-variate polynomial equation by the [factory method](https://en.wikipedia.org/wiki/Factory_method_pattern) which is a design pattern in object-oriented programming. 

We create objects to represent conic sections without having to specify the type of conic section to create. 
In this sense, `Conics` is a classifier of conic sections. 

----

### Definition

A conic section is an algebraic curve of degree two whose coordinates satisfy a quadratic equation in two variables;

$$
Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
$$

with all coefficients $\in \mathbb{R}$ and $A, B, C$ not all zero. This is the general form equation of the conic

The above equation can be written in matrix notation as

$$
\left(\begin{array}{ll}
x & y
\end{array}\right)\left(\begin{array}{cc}
A & B / 2 \\
B / 2 & C
\end{array}\right)\left(\begin{array}{l}
x \\
y
\end{array}\right)+\left(\begin{array}{ll}
D & E
\end{array}\right)\left(\begin{array}{l}
x \\
y
\end{array}\right)+F=0 .
$$

Conic sections described by this equation can be classified in terms of the discriminant $\Delta = B^2 -4AC$.

The discriminant is $-4\Delta$ where $\Delta$ is the quadratic matrix determinant 

$$
\mathrm{det} \textbf{M} = \left|\begin{array}{cc}A & B / 2 \\ B / 2 & C\end{array}\right|
$$.

If the conic is [non-degenerate](https://en.wikipedia.org/wiki/Degenerate_conic) then, 

- if $B^2-4 A C<0$, the equation represents an ellipse.
- - if $A=C$ and $B=0$, the equation represents a circle, a special case of an ellipse.
- if $B^2-4 A C=0$, the equation represents a parabola.
- if $B^2-4 A C>0$, the equation represents a hyperbola.
- if $\tau=A+C=0$, the equation represents a rectangular hyperbola.

----

### Input Parsing and Processing: 
`Conic` takes a string representation of polynomial equation as input. Initial parsing relies on the [SymPy](https://www.sympy.org/en/index.html) library to eliminate fractions and multiply the equation by the LCM thus reducing it to the general form of a conic section.

The `poly_dictionary` program then decomposes the equation into a dictionary of the form; {(degree of x, degree of y): coefficient, ...} with integer coefficients.

----

### Example usage;

```python
mystery_conic = "x^2/4 + x*y/3 + y^2/9 -2*x + 3*y -1"
conic = Conic.create(mystery_conic)

# Type
In [1]: conic.type
Out[1]: Parabola

# Coefficients
In [2]: conic.coeff
Out[2]: {(2, 0): 9, (1, 1): 12, (1, 0): -72, (0, 2): 4, (0, 1): 108, (0, 0): -36}
```
