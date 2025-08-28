4ti2 setup:

## Problem data (with dimensions)

* $B \in \mathbb{Z}^{m\times n}$
* $S \in \mathbb{Z}^{q\times n}$
* $r \in \mathbb{Z}^n$
* Cone: $C=\{r\in\mathbb{R}^n: Br=0,\; Sr\ge 0\}$
* Slice (to force positivity on coordinate $i$): add $r_i \ge 1$.

Let the total number of rows be

$$
d \;=\; m \;+\; q \;+\; 1.
$$

You’ll create four 4ti2 files for a project name, say `slice`. 

---

## Files for `zsolve`

### 1) `slice.mat`  — the coefficient matrix (size $d\times n$)

Header: `d n`
Rows (in this order):

1. The $m$ rows of $B$.
2. The $q$ rows of $S$.
3. One final row $e_i^\top$ with a 1 in column $i$ and 0 elsewhere (to encode $r_i\ge 1$).

```
d n
# -- B block: m rows (equalities)
b11 b12 ... b1n
b21 b22 ... b2n
...
bm1 bm2 ... bmn

# -- S block: q rows (inequalities)
s11 s12 ... s1n
s21 s22 ... s2n
...
sq1 sq2 ... sqn

# -- slice row: r_i >= 1
0 ... 0 1 0 ... 0   # 1 in the i-th position
```

### 2) `slice.rel`  — row-wise relation types (length $d$)

Header: `1 d`
Use:

* `=` for each of the $m$ rows of $B$ (equalities),
* `>` for each of the $q$ rows of $S$ (meaning “$\ge$”),
* `>` for the last slice row $r_i\ge 1$.

You can put them separated by spaces or newlines:

```
1 d
= = = ... =      # m symbols
> > > ... >      # q symbols
>                # 1 symbol for the slice row
```

### 3) `slice.rhs`  — right-hand sides (length $d$)

Header: `1 d`

* Zeros for the $m$ equality rows and the $q$ inequality rows from $S$,
* `1` for the slice row.

```
1 d
0
0
...      # m zeros
0
0
...      # q zeros
1        # for r_i >= 1
```

### 4) `slice.sign`  — variable sign restrictions (length $n$)

Header: `1 n`
Use `0` for each variable to make it **free** (no sign restriction).
(4ti2 uses: `1` = variable $\ge 0$, `-1` = $\le 0$, `0` = free.)

```
1 n
0 0 0 ... 0
```

---

## Run 4ti2

```bash
zsolve slice
```

## What you get back

* `slice.zinhom`: a **finite list** of vectors in $\mathbb{Z}^n$.
  We call this set $T_i$. These are the **minimal inhomogeneous solutions** for the slice; in our terminology, they are exactly the *module generators over the original monoid* for

  $$
  P=\{r:\; Br=0,\; Sr\ge 0,\; r_i\ge 1\}.
  $$
* `slice.zhom`: a (finite) set of **homogeneous generators** $H\subset\mathbb{Z}^n$ for the recession cone (the solutions of $Br=0,\; Sr\ge 0$).
