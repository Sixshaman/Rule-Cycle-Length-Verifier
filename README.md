# Rule-Cycle-Length-Verifier
A tool to verify the maximum cycle lengths of additive cellular automata.

**Statement 1.** In additive cellular automata, any state transition from state $b$ of length $n$ to state $b'$ can be written in matrix form: $b' = Ab$.
Here $A$ is $n \times n$ state transition matrix.

**Statement 2.** Steven Wolfram in [his 1984 paper](https://content.wolfram.com/sw-publications/2020/07/algebraic-properties-cellular-automata.pdf) has shown that any state in an additive cellular automaton eventually enters a cycle. 

In matrix transition terms, it means that for any state $b$ and sufficiently large number $q$, there is a value $p$ such that:  
$A^{p + q} * b = A^q * b$. (eq. 1)

The maximum cycle length corresponds to the largest such value of $p$ for all states $b$.

In practice, choosing $q = 2*p$ is sufficient.

See also [the article](https://arxiv.org/pdf/cond-mat/9305012) by Shin-ichi Tadaki.

**Statement 3.**

If $p$ is a cycle length for state $b$, then $k*p$ for any natural $k > 0$ is also a cycle length for state $b$.  
Likewise, if $p$ and $k$ are cycle lengths and $k$ < $p$, then $k$ must be a divisor of $p$.

This is easy to see from the equation 1.  

**Statement 4.** From the same [Wolfram's paper](https://content.wolfram.com/sw-publications/2020/07/algebraic-properties-cellular-automata.pdf): for cyclic universes, it's sufficient to check a single state $b$ with a single cell "on" .  

I.e. if for such state $b$ the equation $A^{p + q} * b = A^q * b$ holds true, it holds true for all other states.  
Thus, $p$ is the maximum cycle length for the transition matrix $A$.  
In other words, instead of calculating $p$ for all states $b$, we can calculate it for just a single state with a single cell "on".

This follows from Lemma 3.4 in the paper, as well as the statement at the end of the section 2:  
_"By virtue of the cyclic symmetry between the sites it suffices to consider the case_ $j=0$_"._

---

### Rule 90

Cycle lengths for rule 90 were completely described in [Wolfram's paper](https://content.wolfram.com/sw-publications/2020/07/algebraic-properties-cellular-automata.pdf): 
- For $N = 2^k$, the cycle length is always $1$;
- For even $N$ not in the form of $2^k$, the cycle length $p(N)$ is $2 * p(N/2)$;
- For odd $N$, the cycle length divides $2^{sord_2(N)}$, where $sord_2(N)$ is the least power $k$ such that $2^k = \pm1\mod N$.

The $sord_2(N)$ function for odd $N$ corresponds to OEIS sequence [A003558](https://oeis.org/A003558).  
_Note:_ because $sord_2(N)$ is only defined for odd $N$, the value `A003558(N)` is equal to $sord_2(2N + 1)$.

### Rule 150

Rule 150 (evidently) has almost the same behavior as rule 90:
- For $N = 2^k$, the cycle length is always $2^{k-1}$;
- For even $N$ not in the form of $2^k$, the cycle length $p(N)$ is $2 * p(N/2)$;
- For odd $N$, the cycle length divides $2^{sord_2(N)}$.

### Calculating cycle lengths for rule 90 and rule 150 in circle 1D grids

Note that the behavior for rule 90 and rule 150 differs only at powers of 2.  
For other even numbers, the cycle length for $N$ can be obtained by doubling the cycle length for $N/2$, both for rule 90 and rule 150.  
For odd numbers, the cycle length for $N$ divides $2^{sord_2(N)}$, again both for rule 90 and rule 150.

Combining the conditions for all  $N$, we get the following rule: the cycle length divides $2^m * 2^{sord_2(\frac{N}{2^{m}})}$,  
where $m$ is the highest power of $2$ dividing $N$.  

Note that the number $2^m * 2^{sord_2(\frac{N}{2^{m}})}$ divides the cycle length for _all possible_ $N$, both for rule 90 and rule 150.  

This suggests the following algorithm to calculate the cycle length:  
1. (Optional step) Test if $p = 2^m * 2^{sord_2(\frac{N}{2^{m}})}$ really divides the maximum cycle length. This can be done using the statements 2, 3, and 4.  
   Take $q = 2*p$ and check if $A^{p + q} * b = A^p * b$ for a single-cell-on state $b$.  
   If this statement is true, the cycle length for $N$ must be one of the divisors of $p$.
2. (Optional step) Find the precise value of $q$ such that $A^{p + q} * b = A^p * b$ for a single-cell-on state $b$.  
   This can be done using a binary search from $q = 0$ to $q = 2 * p$.  
   If this step is skipped, $q$ is kept at $2 * p$.
3. For all prime divisors $d$ of $p$, check if $p/d$ still divides the maximum cycle length.  
   This is done by checking if $A^{\frac{p}{d} + q} * b = A^{\frac{p}{d}} * b$ for all prime divisors $d$. If the equation still holds, set $p = p / d$.  
   The divisors of $p$ are $m$ copies of $2$, plus all divisors for $sord_2(\frac{N}{2^{m}})$-th Mersenne number.  
   The Mersenne number divisors can be obtained from [A001265](https://oeis.org/A001265) and [A046051](https://oeis.org/A046051).
4. The final number $p$ is the maximum cycle length for the rule given by the transition matrix $A$.

# Implementation  

The application consists of two parts: a C++ part for quick $A^{p + q} * b = A^p * b$ validation and a Python script which calls C++ part with different $N$.  
The C++ part requires Boost (for long integer library `cpp_int` and `dynamic_bitset`).

There are two Python scripts right now:
- `gather_cycle_length_rule_90.py` to calculate maximum cycle lengths for rule 90 from $N = 1$ to $N = 1000$;
- `gather_cycle_length_rule_90.py` to calculate maximum cycle lengths for rule 150 from $N = 1$ to $N = 1000$.
