# Rule-Cycle-Length-Verifier
A tool to verify the maximum cycle lengths of additive cellular automata.

**Statement 1.** In additive cellular automata, any state transition from state $b$ of length $n$ to state $b'$ can be written in matrix form: $b' = Ab$.
Here $A$ is $n \times n$ state transition matrix.

**Statement 2.** Steven Wolfram in [his 1984 paper](https://content.wolfram.com/sw-publications/2020/07/algebraic-properties-cellular-automata.pdf) has shown that any state in an additive cellular automaton eventually enters a cycle. 

In matrix transition terms, it means that for any state $b$ and sufficiently large number $q$, there is a value $p$ such that:
$A^(p + q) * b = (A^q) * b$.

The maximum cycle length corresponds to the largest such value of $p$ for all states $b$.

In practice, choosing $q = 2*p$ is sufficient.

See also [the article](https://arxiv.org/pdf/cond-mat/9305012) by Shin-ichi Tadaki.

**Statement 3.**

If $p$ is a cycle length, then $k*p$ for any natural $k > 0$ is also a cycle length.  
Likewise, if $p$ and $k$ are cycle lengths and $k$ < $p$, then $k$ must be a divisor of $p$.

This is easy to see from the formula in the statement 2.  

---

