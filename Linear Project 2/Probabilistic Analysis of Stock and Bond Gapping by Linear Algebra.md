---
editor_options:
  markdown:
    wrap: 72
output: pdf_document
---

# Probabilistic Analysis of Stock and Bond Overnight Movement by Linear Algebra

**By Sam Lu and Sebastian Reyes**

# General Introduction

This project analyzes stock and bond behavior via data from the SPDR S&P
500 ETF Trust (SPY) and the iShares 20+ Year Treasury Bond ETF (TLT),
respectively. We utilized data ranging from December 2004 to April 2025.

There were two main questions this project hoped to investigate.
Firstly, we aimed to learn whether Markov chains may be used to model
short term and/or long-term stock market behavior. Secondly, we aimed to
learn whether fitting a polynomial via least-squares approximation to
stock market data could provide insights into behavior not observable
from raw data.

At its core, our project researched overnight movements in the stock
market, in which a given stock's price when the market opens is higher
or lower than its previous closing price. We chose to investigate the
relationship between overnight movement in the SPY index and the TLT
index.

# Markov Chains: Introduction

A ***Markov chain*** is defined as an evolving system comprising of a
sequence of ***stages*** where each stage is in one of a finite number
of ***states***, and where the state of a given stage is only dependent
on the state of the stage immediately before it. To phrase it
differently, the Markov property--one of the defining features of a
Markov chain--states that $P(X_{t+1} = s_j \mid X_t = s_i)$, or that the
probability of a given state $X$ at time $t+1$ is only dependent on the
given state X at time $t$. This probability is known as a ***transition
probability***.

***Stochastic matrices***, also commonly known as ***probability
matrices***, are mathematical tools which may be used to represent
Markov chains. They are square matrices of dimension $n$, where $n$ is
the number of possible states in the system. Each column represents a
current state $j$, and each row represents a subsequent state $i$. The
entries of a stochastic matrix $P_{ij}$ are the transition probabilities
corresponding to each state combination. An example stochastic matrix is
displayed in **Figure A**, in which the entry $P_{12}$ indicates that
given a current state of 2, there is a probability of 0.4 (or 40%) that
the next state will be 1. This gives rise to a key property of
stochastic matrices: $\sum_{i} P_{ij} = 1 \; \text{for all } j$. In
other words, the sum of each column must be equal to 1, as there must be
a 100% chance that stage $t+1$ is in one of the $n$ possible finite
states, regardless of the state at stage $t$.

**Figure A:** $$
P = 
\begin{bmatrix}
0.7 & 0.4 \\
0.3 & 0.6 \\
\end{bmatrix}
$$

Any column of a probability matrix may also be isolated to form a
***state vector*** $\overrightarrow{s}$ of dimension $n$. The
probabilities of each state at a stage $m$ can be found by raising each
component of the vector to the power of $m$. $$
\vec{s}^{(m)} = 
\begin{bmatrix}
s_1^{m} \\
\vdots \\
s_n^{m}
\end{bmatrix}
$$ This new vector is known as the ***m***$^{th}$ state vector.
Accordingly, raising $\overrightarrow{s}$ to the power of 0 will result
in the ***initial state vector***. State vectors may also be used to
investigate the long-term behavior of a system. This can be done by
finding the vector $\overrightarrow{s}^{(m)}$ as $m$ approaches
infinity. $$
\lim_{m \to \infty} \vec{s}^{(m)} = \vec{s}
$$ This results in what is known as a ***steady-state vector***.

# Least-Squares Approximation: Introduction

***Least-Squares Approximation*** is defined as a method to fit a line
through a set of data points by minimizing the sum of the squared
distances between the fitted line and each data point. In linear algebra
terms, least-squares approximation is seeking to find a line on
$\mathbb{R}^2$, $\begin{bmatrix} a \\ b \end{bmatrix}$ , that best fits
a series of points $(x_1, y_1)...(x_n, y_n)$. The values $a$ and $b$ are
the slope and intercept of the line, respectively, and the line of "best
fit" is defined as the line that minimizes
$s = (ax_1 + b - y_1)^2 ... + (ax_n + b - y_n)^2$, where $ax_n+b$ is the
value predicted by the line, and $y_n$ is the value of the data point.
In terms more closely aligned with linear algebra, $s$ may also be
defined by the quantity below.

$$
s =
\left\| 
\begin{bmatrix}
ax_1 + b \\
\vdots \\
ax_n + b
\end{bmatrix}
-
\begin{bmatrix}
y_1 \\
\vdots \\
y_n
\end{bmatrix}
\right\|^2
$$

This implies that the below value is to be minimized.

$$
\left\|
\begin{bmatrix}
x_1 & 1 \\
\vdots & \vdots \\
x_n & 1
\end{bmatrix}
\begin{bmatrix}
a \\
b
\end{bmatrix}
-
\begin{bmatrix}
y_1 \\
\vdots \\
y_n
\end{bmatrix}
\right\|^2
$$

Defining
$A = \begin{bmatrix} x_1 & 1 \\ \vdots & \vdots \\ x_n & 1 \end{bmatrix}$,
$\vec{v} = \begin{bmatrix} a \\ b \end{bmatrix}$, and
$\vec{b} = \begin{bmatrix} y_1 \\ \vdots \\ y_n \end{bmatrix}$, this
leaves us with the goal of minimizing the quantity
$\left\| A\vec{v} - \vec{b} \right\|^2$. Knowing this, the establishment
of three main facts must be established prior to continuing:

1.  If $\overrightarrow{v}$ minimizes
    $\left\| A\vec{v} - \vec{b} \right\|^2$,
    $\overrightarrow{b}-A\overrightarrow{v}$ is normal to every vector
    in $im(A)$
    $\Rightarrow (A \overrightarrow{w}) \cdot (\overrightarrow{b} - A \overrightarrow{v}) = 0$.
2.  Given two vectors $\overrightarrow{w_1}$ and $\overrightarrow{w_2}$,
    $\overrightarrow{w_1}\cdot\overrightarrow{w_2}=\overrightarrow{w_1^T}\;\overrightarrow{w_2} \Rightarrow (A \overrightarrow{w})^T(\overrightarrow{b}-A \overrightarrow{v}) = 0 \Rightarrow \overrightarrow{w^T} A^T(\overrightarrow{b}-A \overrightarrow{v}) = 0$.
3.  Given a $\overrightarrow{w}$ such that
    $\overrightarrow{w} \cdot \overrightarrow{w'} = 0$ for **all**
    $\overrightarrow{w'}$,
    $\Rightarrow \overrightarrow{w'} = \overrightarrow{w} = 0 \Rightarrow \overrightarrow{w} \cdot (A^T(\overrightarrow{b} - A \overrightarrow{v}) = 0$
    **UNSURE, CHECK MATH THEORY HERE**

Thus, the ***normal equation*** to be used is defined as
$(A^TA) \overrightarrow{v} = A^T \overrightarrow{b}$.

# Markov Chains: Application

We utilized Python to create a dataframe (an object used to store
tabular data) which displayed the share price of SPY and the share price
of TLT at the time the markets opened and at the time the markets closed
each day. **Figure 1** contains an example of this.

**Figure 1:** $$
\begin{array}{c|cccc}
\text{Date} & \text{SPY Open} & \text{SPY Close} & \text{TLT Open} & \text{TLT Close} \\
\hline
2025\text{-}04\text{-}01 & 557.45 & 560.97 & 107.83 & 108.20 \\
2025\text{-}04\text{-}02 & 555.05 & 564.52 & 108.01 & 109.00 \\
2025\text{-}04\text{-}03 & 545.11 & 537.70 & 109.25 & 108.73 \\
\end{array}
$$

We then used this data to calculate the share price fluctuation that
took place the previous night. This was calculated as
$(Open_t - Close_{t-1})/Close_{t-1}$. For example, the SPY overnight
movement corresponding to April 2nd, 2025 according to **Figure 1**
would be $(555.05-560.97)/560.97$, resulting in a movement of
approximately $-0.011$, or $-1.1\%$. We then normalized these movements
using Z-scores. In this manner, 4 possible states were established,
depending on overnight behavior for that date. These states are
explained in **Figure 2**.

**Figure 2:** $$
\begin{array}{|c|c|c|c|}
\hline
\textbf{State} & \textbf{Stock Movement (Z)} & \textbf{Bond Movement (Z)} & \textbf{Interpretation} \\
\hline
1 & \geq 0 & \geq 0 & \text{Stock up, Bond up} \\
\hline
2 & \geq 0 & < 0 & \text{Stock up, Bond down} \\
\hline
3 & < 0 & \geq 0 & \text{Stock down, Bond up} \\
\hline
4 & < 0 & < 0 & \text{Stock down, Bond down} \\
\hline
\end{array}
$$

Following the establishment of these 4 states, we once again utilized
Python to create a table displaying how many times the different states
"fed into" one another. This table is displayed in **Figure 3**, where
each column represents the current state, and each row represents the
subsequent state. For instance, entry (1, 2) of **Figure 3** indicates
that there were 305 instances where a day categorized as State 2 was
immediately followed by a day categorized as State 1 (i.e. Monday =
State 2, Tuesday = State 1).

**Figure 3:** $$
\begin{array}{c|cccc}
 & \text{State 1} & \text{State 2} & \text{State 3} & \text{State 4} \\
\hline
\text{State 1} & 219 & 305 & 318 & 207 \\
\text{State 2} & 292 & 482 & 560 & 279 \\
\text{State 3} & 317 & 539 & 419 & 250 \\
\text{State 4} & 221 & 287 & 227 & 195 \\
\end{array}
$$

Following the use of **Figure 3**, we opted to scale each column such
that each individual entry would be a proportion. That is, to scale
entry (1, 1), we took the sum of column 1 (219, 292, 317, and 221),
indicating that there were a total of 1049 days in our data which were
categorized as being in State 1. Dividing entry (1, 1) by 1049, we found
that a proportion of 0.2088, or 20.88% of days categorized as State 1
were immediately followed by *another* day categorized as State 1. In
this manner, we reduced the entire table, as displayed in **Figure 4**.

**Figure 4:** $$
\begin{array}{c|cccc}
 & \text{State 1} & \text{State 2} & \text{State 3} & \text{State 4} \\
\hline
\text{State 1} & 0.2088 & 0.1891 & 0.2087 & 0.2223 \\
\text{State 2} & 0.2784 & 0.2988 & 0.3675 & 0.2997 \\
\text{State 3} & 0.3022 & 0.3342 & 0.2749 & 0.2685 \\
\text{State 4} & 0.2107 & 0.1779 & 0.1490 & 0.2095 \\
\end{array}
$$

Slightly changing this table to take the form of a matrix, we found
ourselves with **Figure 5**. This was our final ***probability
matrix***, or ***stochastic matrix***.

**Figure 5:** $$
\begin{bmatrix}  
0.2088 & 0.1891 & 0.2087 & 0.2223 \\  
0.2784 & 0.2988 & 0.3675 & 0.2997 \\  
0.3022 & 0.3342 & 0.2749 & 0.2685 \\  
0.2107 & 0.1779 & 0.1490 & 0.2095  
\end{bmatrix}  
$$

We then utilized this stochastic matrix, in combination with Python, to
find the ***steady-state vector*** for the system, as displayed in
**Figure 6.**

**Figure 6:**

$$
\overrightarrow{s} = \begin{bmatrix}0.2050 \\0.3152 \\0.2980 \\0.1817\end{bmatrix}
$$

This displayed multiple key findings.

1.  $P(1+2) > P(3+4)$. That is, it is more likely that overnight stock
    movement will be positive rather than negative. This is intuitive,
    as stocks tend to rise.\
2.  $P(1+3) \approx P(2+4)$. Once more, this is to be expected, as
    20-year bonds are usually constant/steady.\
3.  $P(2+3) > P(1+4)$. This is logical due to stocks and bonds
    historically being negatively correlated. As such, it is natural
    that divergence is more likely than convergence.
