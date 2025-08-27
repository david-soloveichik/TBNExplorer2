Link to ChatGPT (GPT 5) chat: https://chatgpt.com/c/68aeb3c1-cad8-8323-9adb-8e6ddd1229f2?model=gpt-5

Question: 
My overarching goal is as follows, and I am happy to achieve it any way that I can (efficiently, without "post-filtering" since these are going to be very large problems). I have C, the space of all solutions to the system:
A.x = 0 
S.x >= 0 
x_i > 0
Let C_int be the integer points in C. In other words, C_int = C \intersect Z^n.
I want to find H(C) = {h \in (C_int \ {0}) : h != x + y for any x, y \in (C_int \ {0})}


Answer:
Got it. You want the set of **indecomposable** integer points

$$
H(C)=\{\,h\in\mathbb Z^d\cap C\setminus\{0\}\ :\ h\neq x+y\ \text{for any}\ x,y\in\mathbb Z^d\cap C\setminus\{0\}\,\},
$$

where

$$
C=\{x:\ A x=0,\ Sx\ge 0,\ x_i>0\}.
$$

The most direct (and scalable) way to get exactly that with **Normaliz** is to input your system *as written* and use the built-in “strict” machinery. Here’s the recipe and why it gives you $H(C)$ without any after-the-fact pruning.

---

# How to encode $C$ in Normaliz

### Option A — Work in the original variables $x$

Use **equations**, **inequalities**, and a **strict sign** on coordinate $i$:

```
amb_space d

equations rA
<rows of A>

inequalities rS
<rows of S>

strict_signs
[0 ... 0 1 0 ... 0]   # put 1 in position i

HilbertBasis           # or run: normaliz -N file.in
```

* `equations` encodes $A x=0$; `inequalities` encodes $Sx\ge 0$.
* `strict_signs` with a `1` at position $i$ encodes the **integer** strict condition $x_i>0$ (interpreted as $x_i\ge 1$ on lattice points).&#x20;

### Option B — Reduce dimension first (often faster)

Project to the integer kernel of $A$. Compute an integer basis $U\in\mathbb Z^{d\times k}$ of $\ker(A)\cap\mathbb Z^d$ and set $x=Uy$. Then run Normaliz in the $y$–space:

```
amb_space k

inequalities rS
<rows of (S*U)>

strict_inequalities 1
<row = (e_i^T * U)>        # encodes (e_i U)·y >= 1

HilbertBasis
```

Here `strict_inequalities` encodes linear forms $\alpha\cdot y\ge 1$.&#x20;

---

# What Normaliz returns and how it equals $H(C)$

With any *inhomogeneous/strict* input (as above), the goal **HilbertBasis** makes Normaliz compute **two** things:

1. **Hilbert basis of the recession monoid** (i.e., for the cone $A x=0,\ Sx\ge 0,\ x_i\ge 0$), and
2. **Module generators** that shift you into the strict layer $x_i\ge 1$.&#x20;

From these, your indecomposables are obtained with **no extra computation**:

* Every **module generator** has $x_i=1$, hence cannot split as a sum of two points that each have $x_i\ge 1$. So all module generators lie in $H(C)$.
* Any element of the **recession Hilbert basis** with $x_i>0$ is indecomposable in the (larger) non-strict monoid, hence also indecomposable in the strict one. These are exactly the remaining elements of $H(C)$.

So

$$
H(C)
= \{\text{module generators}\}\ \cup\ \{\,h\in \text{HB(recession)}\ :\ h_i>0\,\}.
$$

Normaliz prints these two lists in separate sections; you just read them off. (By design, in an inhomogeneous run **HilbertBasis** means “Hilbert basis of the recession monoid **and** module generators”.)&#x20;

> Note: If your recession cone has directions with $x_i=0$, Normaliz will still list those rays in the **recession** Hilbert basis—this is necessary to describe the arithmetic of the strict layer—but they are **not** candidates for $H(C)$ since they violate $x_i>0$.

---

# Performance tips (for large instances)

* **Reduce dimension** via the integer kernel $U$ (Option B) when $A$ has sizable rank; it can be a big win.
* Pick algorithm: **DualMode** (`-d`) is often faster when you have many inequalities; **PrimalMode** (`-P`) can win when the cone has few extreme rays. You can try both.&#x20;
* Make sure the **recession cone** is as pointed as possible after applying $A$ (Normaliz handles nonpointed cones but will also report a basis for the lineality space separately). If lineality is present, the output still splits cleanly into recession HB + module + lineality basis (as many interfaces document). ([doc.sagemath.org][1])

---

[1]: https://doc.sagemath.org/html/en/reference/discrete_geometry/sage/geometry/polyhedron/backend_normaliz.html "The Normaliz backend for polyhedral computations - Combinatorial and Discrete Geometry"
