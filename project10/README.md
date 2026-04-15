# Introduction
This project was completed as a group and tasked us with implementing the Baum-Welch algorithm for parameter estimation in Hidden Markov Models (HMMs). The Baum-Welch algorithm is an Expectation-Maximization (EM) approach that allows an HMM to learn its own parameters — initial state probabilities, transition probabilities, and emission probabilities — directly from unlabeled observation sequences, without needing any ground-truth state labels.

Our implementation takes in four data structures (observation sequences, initial probabilities, transition probabilities, and emission probabilities) and uses them to iteratively refine the model until it converges on a locally optimal set of parameters. We made no assumptions about the number of hidden states, the size of the observation alphabet, or the domain of the data, keeping our solution fully general.

To structure our work, we used Object-Oriented Programming (OOP). Having already implemented Viterbi, Forward, Backward, and Posterior decoding in our existing HMM class, the Baum-Welch algorithm was a natural extension — the E-step draws directly on our Forward and Backward methods to compute expected counts, and the M-step uses those counts to re-estimate parameters, completing the full unsupervised learning pipeline.
# Pseudocode

```
## Pseudocode

### `_safe_log(x)`
- Input: probability value `x`
- If `x <= 0`
  - Return `-∞`
- Else
  - Return `log(x)`

---

### `_synchronize_emission()`
- For each state in the HMM:
  - For each emission in that state:
    - Add emission to global emission set

- For each emission in global emission set:
  - For each state:
    - If state does not contain emission:
      - Add emission to state
      - Set emission probability to `0`

---

### `_expectation_values(observations)`
- Convert observations to list
- Let `n_states = number of hidden states`
- Let `T = length of observations`

- Compute `forward_mat`
- Compute `backward_mat`

- Compute sequence log-likelihood:
  - `log_prob = log-sum-exp(last column of forward_mat)`

- Compute gamma:
  - `gamma = exp(forward_mat + backward_mat - log_prob)`

- Initialize `xi` as zeros of shape `(T-1, n_states, n_states)`

- For each time step `t` from `0` to `T-2`:
  - For each state `i`:
    - For each state `j`:
      - Compute:
        - forward probability at `(i, t)`
        - transition probability `i → j`
        - emission probability of next symbol from state `j`
        - backward probability at `(j, t+1)`
      - Combine into log-xi
      - Convert to normal probability
      - Store in `xi[t, i, j]`

- Return:
  - forward matrix
  - backward matrix
  - log-likelihood
  - gamma
  - xi

---

### `sequence_log_likelihood(observations)`
- Call `_expectation_values(observations)`
- Extract `log_prob`
- Return `log_prob`

---

### `gamma_matrix(observations)`
- Call `_expectation_values(observations)`
- Extract `gamma`
- Return `gamma`

---

### `xi_tensor(observations)`
- Call `_expectation_values(observations)`
- Extract `xi`
- Return `xi`

---

### `posterior_decoding(observations)`
- Compute `gamma = gamma_matrix(observations)`

- Convert gamma to log-space:
  - If `gamma > 0` → `log(gamma)`
  - Else → `-∞`

- For each position:
  - Select state with highest posterior probability

- Convert indices to state objects

- Return:
  - posterior matrix
  - state path

---

### `baum_welch(sequences, max_iter, tol, pseudocount)`

- Let `n_states = number of hidden states`
- Create sorted list of emission symbols
- Create dictionary mapping symbol → index

- Initialize empty list `history`
- Set `prev_loglik = -∞`

- For each iteration up to `max_iter`:

  - Initialize:
    - `beta_counts = 0`
    - `trans_counts = 0`
    - `trans_denom = 0`
    - `emit_counts = 0`
    - `emit_denom = 0`
    - `total_loglik = 0`

  - For each sequence:
    - Convert sequence to list
    - Compute `(log_prob, gamma, xi)`

    - Add `log_prob` to total log-likelihood

    - Update initial state counts:
      - `beta_counts += gamma[:, 0]`

    - Update transition counts:
      - `trans_counts += sum(xi over time)`
      - `trans_denom += sum(gamma except last column)`

    - For each position `t`:
      - Identify observed symbol
      - Find its index
      - Add `gamma[:, t]` to emission counts

    - Update emission denominator:
      - `emit_denom += sum(gamma)`

  - Update initial probabilities:
    - Normalize `beta_counts` with pseudocount

  - Update transition probabilities:
    - Normalize `trans_counts` with pseudocount

  - Update emission probabilities:
    - Normalize `emit_counts` with pseudocount

  - Append `total_loglik` to history

  - If improvement < `tol`:
    - Stop training

  - Update `prev_loglik`

- Return `history`
```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Tien Nguyen
I found the Baum–Welch algorithm challenging to grasp at first, particularly the idea of learning hidden state behavior without directly observing it. While the overall algorithm felt abstract, the gamma and xi computations were more concrete and easier to implement. As I connected these computations back to the E-step, I began to understand how the model uses probabilities instead of fixed paths to learn from the data. This process helped me develop a deeper understanding of how expectation-based methods work in probabilistic models.

# Generative AI Appendix
CHATGPT was used to explain the Baum–Welch algorithm, assist with pseudocode planning, and support debugging throughout the implementation process.