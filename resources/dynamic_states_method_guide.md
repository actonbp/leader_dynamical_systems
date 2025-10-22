Dynamic States: Method Selection Guide
=====================================

Purpose
-------
Brief guidance on when to use derivative‐based operationalizations (velocity/acceleration) versus model‐based parameters (ω, ζ) for leadership identity dynamics.

Derivative Proxies (Velocity/Acceleration)
------------------------------------------
- What they are: Local, timescale‐specific features derived from first and second derivatives of within‐person residuals (via GOLD).
- Use when:
  - Short or sparse time series (e.g., 7 daily observations per week).
  - Questions focus on short‑term “ebb/flow” and “speed of return”.
  - You want a transparent, low‐assumption pipeline and descriptive features.
- Inference:
  - Aggregate per ID × context (e.g., mean |velocity|, mean acceleration).
  - Use multilevel models for population‐level estimates and moderators.
- Caveats:
  - Proxies, not structural parameters; sensitive to sampling, embedding (tau), and imputation.
  - Acceleration mixes restoring and damping effects unless modeled with x and x′.
  - Report sensitivity (tau 3–5; linear vs spline; optional normalization like mean|v|/sd(x)).

Model‑Based Parameters (ω, ζ)
----------------------------
- What they are: Natural frequency (ω) and damping ratio (ζ) from a second‐order ODE: x'' = −2ζω x' − ω² x + β u + ε.
- Use when:
  - You need interpretable, design‑invariant dynamics for macro‑level theory or cross‐study comparability.
  - You have longer/denser series or can estimate across multiple windows.
  - You need explicit input modeling (events) and error structure (e.g., AR(1)).
- Inference:
  - Estimate per person/context, then pool in hierarchical models to obtain population means/variances and moderators.
  - Report units (ω depends on time scale), and boundary cases (overdamped, unstable fits).

Timescale and Design Dependence
-------------------------------
- Derivatives and ω depend on the time unit (per day vs per minute). Be explicit about sampling and embedding.
- If identity dynamics vary across time (“change in dynamics”), allow time‑varying windows (derivatives) or time‑varying parameters (state‑space/ctsem).

Practical Recommendation
------------------------
- Use the derivative approach as the default for shortitudinal identity data; it is valid and defensible for micro‑dynamics.
- Add an ODE parameterization as a robustness check when you need interpretable dynamics.
- Always document timescale, tau, imputation choice, and sensitivity outcomes.

