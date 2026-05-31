# Kalman Filter and RTS Smoother (TASC)

Implements the Kalman filter (forward pass) and Rauch-Tung-Striebel
smoother (backward pass) for the state-space model in Rho et al. (2026):

## Usage

``` r
kalman_smoother_cpp(Y, W, A, C, Q, R, z0, P0)
```

## Arguments

- Y:

  Observed data matrix (N x T). Use NA for unobserved entries.

- W:

  Observation / loading matrix (N x r)

- A:

  State transition matrix (r x r). Pass diag(r) for random-walk
  dynamics.

- C:

  State drift vector (r x 1)

- Q:

  State noise covariance (r x r)

- R:

  Observation noise covariance (N x N, diagonal in practice)

- z0:

  Initial state mean (r x 1)

- P0:

  Initial state covariance (r x r)

## Value

A list with z_smooth, P_smooth, P_cross, z_pred, z_upd. P_cross is an r
x r x (T-1) cube. Slice t (C++ 0-indexed, t=0,...,T-2) stores P(t+1, t
\| T) (0-indexed), i.e. P(t+2, t+1 \| T) in 1-indexed Shumway-Stoffer
notation. Formula: P(t+1\|T) \* J_t^T (eq. 6.68-6.69).

## Details

State: z(t+1) = A z(t) + C + eta(t), eta(t) ~ N(0, Q) Observation: y_t =
W z_t + eps_t, eps_t ~ N(0, R)

Observation rows with NA (treated post-intervention) are automatically
dropped at each time step so only control-unit rows update the filter.

The P update uses the numerically stable Joseph form: P(t\|t) = (I - K
W_obs) P(t\|t-1) (I - K W_obs)^T + K R_obs K^T
