# Examples

A number of examples with use cases is provided.

## LSM examples

In this subsection all example code use **Least-Squares Method (LSM)** for reconstruction.

### Example_FHN
This example illustrates reconstructing FitzHugh-Nagumo system from 20 random data points. Running code as is, we obtain the following problem solution:
```
f1 =  x2
f2 = - x1 + 4*x2 - x2^3
```
This result is obtained when a derivative point set $W$ is obtained analytically:
```matlab
w = transpose(FHN(0,transpose(y)));
```
You can change 1 to 0 in the following code line:
```matlab
derivativecalc = 1; % Set 1 to find derivatives analytically, set 0 to find derivatives numerically
```
After that, the derivative will be found numerically with a 4-th order finite difference derivative:
```matlab
w = diff4(y,t); 
```
The result will lose its accuracy showing high sensitivity of the algorithm to numerical and truncation errors:
```
f1 = 0.99996*x2
f2 = -0.99811*x1 + 3.9906*x2 - 0.99748*x2^3
```
Actual phase portrait looks like

![FHNPhasePortrait](https://github.com/worriedlemon/oderecon/blob/main/GITHUB_GRAPHICS/FHN.png)

Adding strict accuracy options for the ode solver and replacing the `ode45` solver to the more accurate `ode78` or `ode113` improves this result, please, check it out.

### Example_Linear2

This code illustrates reconstruction of a stable 2-dimensional linear system with exponential decay. Derivatives are found analytically. The output of the code is
```
f1 =  x2
f2 = - x1 - 0.3*x2
```
### Example_Lorenz

Reconstructing the classical Lorenz attractor with a workflow shown earlier. The result of the code is
```
f1 = -10*x1 + 10*x2
f2 = 28*x1 - x2 - x1*x3
f3 = -2.6667*x3 + x1*x2
```

Actual phase portrait looks like

![LorenzPhasePortrait](https://github.com/worriedlemon/oderecon/blob/main/GITHUB_GRAPHICS/Lorenz.png)

### Example_Lorenz_1var

Reconstructing the Lorenz attractor given only one state variable $z$. For better results, `ode113` solver is used, and additional options for tolerances are added:
```matlab
opts = odeset('RelTol',1e-13,'AbsTol',1e-15);
[t,y] = ode113(@Lorenz2,[0:h:Tmax],[0.1,0,-0.1],opts); %solve ODE
```
Once we have only $z$ variable, we must reconstruct two other variables. Usually, a direct application of Takens's theorem is performed, and the system is reconstructed in delay variables. Alternatively, derivatives are also often used. Here, we reconstruct two missing variables via numerical integration of $z$ using Bool's rule:
```matlab
val = u(:,3); %3rd variable
[iv, ind] = integrate_bool(val,0,h);
[iv2, ind2] = integrate_bool(iv,0.1,4*h);
```
We can find analytically that the resulting formula should contain negative powers, so we generate a degree-lexicographic ordering starting from $-1$:
```matlab
sigma = deglexord(-1,3,3);
```
Then, we do not use `PolyRegression` but explicitly apply `delMinorTerms` to each equation to find $H$ and $T$. This is done just for illustration, and all that bulcky code can be rewritten with `PolyRegression`. The result of the code is
```
f1 =  x2
f2 =  x3
f3 = 720*x1 - 29.3333*x2 - 13.6667*x3 + 11*x1^-1*x2^2 + x1^-1*x2*x3 - 10*x1^3 - x1^2*x2
```
For more detail, we refer to the [original publication](https://doi.org/10.3390/math8020300).

### Example_Mem3var

Reconstructing B. Muthuswamy’s circuit equations from a three-dimensional data series, with one missing state variable. The result of the code is
```
f1 =  x2
f2 = 24519*x2 + 73530*x3 - 0.086055*x1^2*x2
f3 = 7353*x2 - 7353*x3 - 14700000*x4
f4 = 55.55*x3
```
For more detail, we also refer to the [original publication](https://doi.org/10.3390/math8020300).

### Example_RingTest

This example shows how ODERECON deals with a conservative system. The output of the code is
```
f1 = -9.9865e-05*x2
f2 = 1.0006*x1
```
The inaccurate coefficients are due to numerical differentiation `dX = diff(X,1,1)`. Replacing it with the following code:
```matlab
dX = diff4(X,tspan)
```
results in accurate equations:
```
f1 = -0.0001*x2
f2 =  x1
```

## Orthogonal polynomials

In this subsection all example code (files end with `_Orth` suffix) use **Orthogonal polynomials properties** for reconstruction.

### Example_1dim_Orth and Example_2dim_Orth

These examples are simple tests of the proposed method. They contain algortihm steps for reconstructing a simple algebraic equation

$$f(x) = 1.7x^3 - x^2 - 0.6x + 4$$

and a Himmelblau function

$$f(x,y) = (x^2 + y^2 - 11)^2 + (x + y^2 - 7)^2$$

respectively (`1dim` and `2dim`).

### Example_FHN_Orth

This example illustrates reconstructing FitzHugh-Nagumo system (see [Example_FHN](https://github.com/worriedlemon/oderecon/blob/main/examples/EXAMPLES.md#example_fhn)) from time-series using `Orthogonal polynomials`.

Firstly, the variable `delminor` indicates whether *delMinorTerms* algorithm should be applied to data.

```matlab
delminor = 1;     % delMinorTerms is being used
```

Then, the second important line is 

```matlab
F = orthpoly_t(sigma, t, x);     % Getting orthogonal polynomials matrix
```

Function `orthpoly_t` creates an orthogonal set of monomials $\theta_i(\mathbf{x})$ over time series `t`, where every single monomial is actually a polynomial of regular monomials with certain coefficients. Matrix `F` represents a linear map from regular polynomials $\mathbb{O}$ to orthogonal polynomials $\Theta$ and has a size of `mc`$\times$`mc` (or $\mathbb{|O|} \times \mathbb{|O|}$). `mc` is equal to `sizeof(sigma, 1)`. Matrix `F` is described in the [Algorithm](https://github.com/worriedlemon/ODERECON/tree/main#algorithm) section of file README.md.

Using matrix `F` we can simply obtain coefficients of each monomial by using Fourier series:

```matlab
E = EvalPoly(F', x, sigma);
for j = 1:mc
    for i = 1:eqc
        coefs(j, i) = trapz(x(:, i), E(:, j));
    end
end
```

A process is a little bit different for `delminor = 1` - there we use the following code

```matlab
y = diff4(x,t);     % derivatives
tol = 1e-3;         % tolerance
for j = 1:vc
    [coefs(:, j), ~, ~, coefs_reg(:, j)] = delMinorTerms_dy(t,x(:, j), x, y(:, j), F, sigma, tol, 0);
end
```

With the help of the matrix `F` we can also obtain regular polynomials' coefficients by simply multiplying the matrix transposed on the vector of orthogonal polynomials.

```matlab
coefs_reg =  F' * coefs;
```

### Example_Rossler_Orth

The process is similar to FitzHugh-Nagumo reconstruction, but the data series comes from the Rossler attractor.
