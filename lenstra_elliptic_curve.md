# Lenstra's Elliptic Curve Factorization Method

## History

Lenstra's Elliptic Curve Method was introduced in a paper written by Hendrik Lenstra, a Dutch mathematician, in 1987 [^1]. The method is inspired by Pollard's $p-1$ Method, but utilizes a group of points on elliptic curves defined by a finite field \cite{Stein_2009}.

## Mathematical Overview

We note the definition of an elliptic curve in $\mathbb{R}^2$ space:

$$y^2 = x^3 + ax + b$$

We then consider the discrete analogue for a finite field with $N$ elements, namely $\mathbb{Z}_N$. In this case, we say that coefficients $a,b \in \mathbb{Z}_N$. Let us call this group $E(N)$. We also set a restriction that $4a^3 + 27$ is invertible modulo $N$, ensuring that there are no singular points on the elliptic curve, which would result in the inability to draw a tangent line. Drawing tangent lines is crucial to the definition of addition under the elliptic curve group $E(N)$.

Along with all points on the elliptic curve, we also consider a new, general, point called the "point at infinity", denoted by $\mathcal{O}$. For our purposes, we do not utilize $\mathcal{O}$, but nonetheless, it is recognized in point addition.

We define a group operator $+$ on $E(N)$ as an algorithm:

### **Algorithm 1**: Group $+$ Operator over $E(N)$

>Given $P,Q \in E(N)$, we find $R=P+Q=Q+P$. Note that $(x_1.y_1)=P$ and $(x_2,y_2)=Q$.
>1. Is either $P=\mathcal{O}$ or $Q=\mathcal{O}$? Return $R=Q$ if the former and $R=P$ if the latter. Return >$R=\mathcal{O}$ if both.
>
>2. Is $x_1=x_2$ and $y_1=-y_2$? Return $R=\mathcal{O}$.
>
>3. Is $P=Q$? Then $\lambda=(3{x_1}^2+a)/(2y_1)$
>
>4. Is $P \neq Q$? Then $\lambda = (y_1-y_2)/(x_1-x_2)$
>
>5. Return $R=(\lambda^2-x_1-x_2, -\lambda(\lambda^2-x_1-x_2)-(y_1-\lambda x_1))$

Note that addition is commutative in $E(N)$. We examine the multiple cases of addition between points in $E(N)$, and how they are formulated.

* For any $P \in E(N)$, it follows that $P + \mathcal{O} = P$.

* For any $P,Q \in E(N)$, if $P \neq Q$, then $R$, the sum, is the point that is at the intersection between the elliptic curve and a straight line drawn between $P$ and $Q$ which is then flipped across the $x$-axis. If this point does not exist, the result is $\mathcal{O}$. This is why if they are mirrored on the x-axis, the line does not intersect a point on the curve, so the result is $\mathcal{O}$.

* For any $P,Q \in E(N)$, if $P = Q$, then $R$, the sum, is the point that is at the intersection between the elliptic curve and a straight tangent line from $P$ which is then flipped across the $x$-axis. If this point does not exist, the result is $\mathcal{O}$.

Of particular note is the fact that, when calculating $\lambda$, we require the utilization of the inverse of the denominator modulo $N$. If $N$ is a prime, this is always guaranteed so long as the denominator is non-zero. However, if $N$ is composite, we only have $\phi(N) \neq N-1$ invertible elements, so it becomes possible for a non-zero denominator to not be invertible. We check before calculating by ensuring that the greatest common divisor between our denominator and $N$ is 1. Say we find a denominator $\alpha$ such that $\text{gcd}(\alpha, N) = d \neq 1$. Then if $d \neq N$, we have a non-trivial factor of $N$. This is the basis of factorization with elliptic curves.

How do we find such points $P$ and $Q$ such that addition will allow us to encounter this special algorithm-breaking case? We can take an initial elliptic curve that is guaranteed to have a starting point $P=(0,1)$, and compute $B!P$:

$$B!P = P + P + \cdots + P$$

This would require that $b = 1$, which is why $b$ is never mentioned in the calculations in \textbf{Algorithm 3}. If we ever encounter an addition which causes the division by an invertible element $\alpha$, then we do as mentioned above and ensure that $1<d<N$ and return $d$. If $d=N$, we were unlucky, and can try with a new initial value $a$, essentially creating a new elliptic curve. We, therefore, have two bounds: $B!$ for the number of point additions tested and $K$ for the number of random elliptic curves tested.

### **Algorithm 2**: Lenstra's Elliptic Curve Factorization Method

>Given integer $N\geq2$ and bound $B\geq1$, we find either non-trivial factor of $N$ or fail.
>1. Select random $a \in \mathbb{Z}_N$ such that $\gcd(4a^3+27, N) = 1$. Our elliptic curve is therefore $y^2=x^3+ax+1$, and point $P=(0,1)$ exists on the curve.
>
>2. Attempt to compute $B!P$. If a denominator $\alpha$ is ever encountered in either steps 3 or 4 of [**Algorithm 1**](#algorithm-1-group--operator-over) such that $\gcd(\alpha, N) = d \neq 1$, then check if $d \neq N$, and return if true. Otherwise, return failure or test with new $a$. If all denominators are relatively prime to $N$, return failure.

## Implementation

We implement [**Algorithm 1**](#algorithm-1-group--operator-over) and [**Algorithm 2**](#algorithm-2-lenstras-elliptic-curve-factorization-method) in [```lenstra_elliptic_curve.py```](./lenstra_elliptic_curve.py). Note that we provide an additional parameter for **Algorithm 2**, namely a number of times a random elliptic curve is generated for testing.


[^1]: My reference, with further explanation and a [supporting link](https://website.com).