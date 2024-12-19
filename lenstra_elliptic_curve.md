# Lenstra's Elliptic Curve Method

## History

Lenstra's Elliptic Curve Method was introduced in a paper written by Hendrik Lenstra, a Dutch mathematician, in 1987 [^1]. The method is inspired by Pollard's $p-1$ Method, but utilizes a group of points on elliptic curves defined by a finite field \cite{Stein_2009}.

[^1]: My reference, with further explanation and a [supporting link](https://website.com).

## Mathematical Overview

We note the definition of an elliptic curve in $\mathbb{R}^2$ space:

$$y^2 = x^3 + ax + b$$

We then consider the discrete analogue for a finite field with $N$ elements, namely $\Z_N$. In this case, we say that coefficients $a,b \in \Z_N$. Let us call this group $E(N)$. We also set a restriction that $4a^3 + 27$ is invertible modulo $N$, ensuring that there are no singular points on the elliptic curve, which would result in the inability to draw a tangent line. Drawing tangent lines is crucial to the definition of addition under the elliptic curve group $E(N)$.

Along with all points on the elliptic curve, we also consider a new, general, point called the "point at infinity", denoted by $\mathcal{O}$. For our purposes, we do not utilize $\mathcal{O}$, but nonetheless, it is recognized in point addition.

We define a group operator $+$ on $E(N)$ as an algorithm:
\begin{algorithm}[H]
\caption{Group $+$ Operator over $E(N)$}\label{alg:cap}
\begin{algorithmic}
\Require $N \geq 2, P \in E(N), Q \in E(N)$
\If{$P = \mathcal{O}$ or $Q = \mathcal{O}$}
    \State $R \gets Q$ \textbf{if} $P = \mathcal{O}$ \textbf{else} $P$
    \State\Return R \Comment{$P + \mathcal{O} = P$}
\EndIf
\State $(x_1,y_1) \gets P$ \Comment{$x_1,y_1 \in \Z_N$}
\State $(x_2,y_2) \gets Q$ \Comment{$x_2,y_2 \in \Z_N$}
\If{$x_1 = x_2$ and $y_1 = -y_2$}
    \State $R \gets \mathcal{O}$
    \State\Return R
\EndIf
\If{$P = Q$}
    \State $\lambda \gets (3{x_1}^2 + a)/(2y_1)$ \Comment{$\lambda$ is slope of tangent line to $P$}
\Else
    \State $\lambda \gets (y_1-y_2)/(x_1-x_2)$ \Comment{$\lambda$ is slope of line intersecting $P$ and $Q$}
\EndIf
\State $x_3 \gets \lambda^2 - x_1 - x_2$
\State $\nu = y_1 - \lambda x_1$
\State $y_3 \gets -\lambda x_3 - \nu$
\State $R \gets (x_3, y_3)$
\State\Return $R$
\end{algorithmic}
\end{algorithm}

Note that addition is commutative in $E(N)$. We examine the multiple cases of addition between points in $E(N)$, and how they are formulated.

\begin{itemize}
    \item[-] For any $P \in E(N)$, it follows that $P + \mathcal{O} = P$.
    \item[-] For any $P,Q \in E(N)$, if $P \neq Q$, then $R$, the sum, is the point that is at the intersection between the elliptic curve and a straight line drawn between $P$ and $Q$ which is then flipped across the $x$-axis. If this point does not exist, the result is $\mathcal{O}$.
    \item[-] For any $P,Q \in E(N)$, if $P = Q$, then $R$, the sum, is the point that is at the intersection between the elliptic curve and a straight tangent line from $P$ which is then flipped across the $x$-axis. If this point does not exist, the result is $\mathcal{O}$.
\end{itemize}

Of particular note is the fact that, when calculating $\lambda$, we require the utilization of the inverse of the denominator modulo $N$. If $N$ is a prime, this is always guaranteed so long as the denominator is non-zero. However, if $N$ is composite, we only have $\phi(N) \neq N-1$ invertible elements, so it becomes possible for a non-zero denominator to not be invertible. We check before calculating by ensuring that the greatest common divisor between our denominator and $N$ is 1. Say we find a denominator $\alpha$ such that $\text{gcd}(\alpha, N) = d \neq 1$. Then if $d \neq N$, we have a non-trivial factor of $N$. This is the basis of factorization with elliptic curves.

How do we find such points $P$ and $Q$ such that addition will allow us to encounter this special algorithm-breaking case? We can take an initial elliptic curve that is guaranteed to have a starting point $P=(0,1)$, and compute $B!P$:

\[B!P = P + P + \cdots + P\]

This would require that $b = 1$, which is why $b$ is never mentioned in the calculations in \textbf{Algorithm 3}. If we ever encounter an addition which causes the division by an invertible element $\alpha$, then we do as mentioned above and ensure that $1<d<N$ and return $d$. If $d=N$, we were unlucky, and can try with a new initial value $a$, essentially creating a new elliptic curve. We, therefore, have two bounds: $B!$ for the number of point additions tested and $K$ for the number of random elliptic curves tested.

\begin{algorithm}[H]
\caption{Lenstra's Elliptic Curve Method}\label{alg:cap}
\begin{algorithmic}
\Require $N \geq 2, B \geq 1, K \geq 1$
\State $k \gets 0$
\While{$k < K$}
\State $a \gets \text{Unif}(2,N-1)$ \Comment{Random integer $1<a<N$}
\If{$\text{gcd}(4a^3+27,N) \neq 1$}
\State \textbf{continue}
\EndIf
\State Compute $B!P$
\If{denominator $\alpha$ is ever encountered such that $\text{gcd}(\alpha,N) \neq 1$}
\State $d \gets \text{gcd}(\alpha,N)$
\If{$1<d<N$}
\State\Return d
\EndIf
\EndIf
\State $k \gets k + 1$
\EndWhile
\end{algorithmic}
\end{algorithm}

Note that, similar to \textbf{Algorithm 1} and \textbf{Algorithm 2}, $a$ is drawn from a discrete uniform distribution.

\subsection{Implementation}

We implement \textbf{Algorithm 4} in Sage.

\begin{code}
\caption{Lenstra's Elliptic Curve Method in Sage (Python)}
\begin{minted}[linenos,frame=lines]{python}
import random

def elliptic_curve_add(P, Q, a, N):
    (x1, y1) = P
    (x2, y2) = Q
    
    if P == Q:
        d = gcd(2*int(y1), N)
        if d != 1:
            return (False, None, d)

        l = (3*x1^2 + a)/(2*y1)
        
    else:
        d = gcd(int(x1)-int(x2), N)
        if d != 1:
            return (False, None, d)
        
        l = (y1-y2)/(x1-x2)
       
    v = y1-l*x1
    x_ = l^2-x1-x2
    R = (x_, -l*x_^3-v)
    return (True, R, None)
            
    

def lenstra_elliptic_curve(N, B, K):
    for _ in range(K):
        P = (Mod(0,N), Mod(1,N))
        while True:
            a = Mod(random.randint(2,N), N)
            if gcd(4*a^3 + 27, N) == 1:
                break
        
        P_ = P
           
        for b in range(B):
            if P == None:
                    break
            for _ in range(b):
                if P == None:
                    break
                (work, P, d) = elliptic_curve_add(P, P_, a, N)
                if not work:
                    if 1 < d and d < N:
                        return d
    
    return None
\end{minted}
\end{code}

Lines (3) through (24) define the addition of points in $E(N)$ as described in \textbf{Algorithm 3}. We initialize a point $P$ in line (30). We then find a suitable coefficient $a$ in lines (31) to (35). We then keep track of the initial point in line (37). We then perform nested additions of increasing count until we achieve $B!P$, employing the function defined for point addition. If this function ever returns an output stating that computation was not possible due to a non-invertible denominator, the greatest common divisor between the denominator and $N$ is returned. The main algorithm then recognizes the failure on line (46), checks if the divisor is non-trivial on line (47), and returns it on line (48). If not, it continues computing $B!P$. If the computation is ultimately successful, the algorithm cycles $K$ times maximum, utilizing different coefficients $a$.

We test the algorithm with large number factorizations, utilizing the same five (50 digit) numbers as in Listing 2 and Listing 5.

\begin{code}
\caption{Lenstra's Elliptic Curve Method Testing Input in Sage (Python)}
\begin{minted}[linenos,frame=lines]{python}
nums = [20450795081324708633669281467532419890072792846123,
        61609551856222489159014922280716157711388198115973, 
        28507764660819497990602573735852708834963899025723, 
        22687481133480947572271781795641028016535047861045, 
        35682687494674502027123657314702931270785057320590]

for n in nums:
    print(lenstra_elliptic_curve(n, 100, 10))
\end{minted}
\end{code}

Our output is the following:

\begin{code}
\caption{Lenstra's Elliptic Curve Method Testing Output in Sage (Python)}
\begin{minted}[linenos,frame=lines]{python}
521
863
7
5
2
\end{minted}
\end{code}