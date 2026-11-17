# Solving the Damped Wave Equation with the FFT

In the previous assignment, we used the FFT to solve the diffusion equation by allowing the coefficients of a trigonometric interpolant to change with time. Each Fourier coefficient satisfied a first-order ordinary differential equation that we could solve analytically.

We will now extend this approach to the one-dimensional **damped wave equation**. We will keep the same periodic domain and spatial grid. However, the equation now involves a second time derivative, so we must specify both an initial displacement and an initial velocity. We will also allow an external forcing term. For a general forcing term, we may not have a convenient analytical solution to the coefficient ODEs, so we will solve them numerically.

As in the previous assignments, you are encouraged to use AI tools to assist with programming. However, you are responsible for ensuring that your implementation is mathematically correct and for understanding and explaining your results.

Consider

$$
\frac{\partial^2 u}{\partial t^2}+\alpha\frac{\partial u}{\partial t}
=c^2\frac{\partial^2 u}{\partial x^2}+h(x,t),
$$

on $0\leq x<2\pi$, with periodic boundary conditions and initial conditions

$$
u(x,0)=f(x),\qquad u_t(x,0)=g(x).
$$

Here $c>0$ is the wave speed, $\alpha\geq0$ is the damping coefficient, and $h(x,t)$ is the forcing. Unlike the diffusion assignment, $\alpha$ now measures damping rather than diffusion. Assume the initial data and forcing are periodic.

## Part 1. Implementing and testing RK4

Before solving the PDE, we will test our time-stepping method on a single second-order ODE:

$$
q''(t)+\alpha q'(t)+a q(t)=b(t).
$$

1. Introduce $v(t)=q'(t)$ and rewrite this equation as a first-order system for $Y(t)=(q(t),v(t))^T$.

2. Implement the classical fourth-order Runge-Kutta method for a general system $Y'=F(t,Y)$:

$$
\begin{aligned}
   K_1&=F(t_n,Y_n),\\
   K_2&=F(t_n+\Delta t/2,Y_n+(\Delta t/2)K_1),\\
   K_3&=F(t_n+\Delta t/2,Y_n+(\Delta t/2)K_2),\\
   K_4&=F(t_n+\Delta t,Y_n+\Delta t K_3),\\
   Y_{n+1}&=Y_n+\frac{\Delta t}{6}(K_1+2K_2+2K_3+K_4).
\end{aligned}
$$

   Your implementation should work for a vector of unknowns so that you can reuse it for the PDE.

3. Test your code with $\alpha=0.5$, $a=4$, $b(t)=-\sin(2t)$, $q(0)=1$, and $q'(0)=0$. Verify by substitution that $q(t)=\cos(2t)$ is the exact solution.

4. Compute the error in $q(1)$ using $\Delta t=1/5,\ 1/10,\ 1/20,\ 1/40,$ and $1/80$. Plot the error against $\Delta t$ on a log-log scale. Estimate the observed order using

$$
p=\log_2\!\left(\frac{E_{\Delta t}}{E_{\Delta t/2}}\right).
$$

   Does your implementation exhibit fourth-order convergence?

## Part 2. Constructing the Fourier method

As in the diffusion assignment, use $N=2m$ equally spaced points

$$
x_j=\frac{2\pi j}{N},\qquad j=0,\ldots,N-1,
$$

and seek an approximation of the form

$$
p_N(x,t)=\sum_{k=-(m-1)}^{m-1}a_k(t)e^{ikx}+a_m(t)\cos(mx).
$$

As before, the cosine represents the single Nyquist entry in the FFT. On the grid, $e^{imx_j}=e^{-imx_j}=\cos(mx_j)=(-1)^j$.

1. Compute $\partial p_N/\partial t$, $\partial^2p_N/\partial t^2$, and $\partial^2p_N/\partial x^2$. Represent the sampled forcing by a trigonometric interpolant in the same basis, with coefficients $H_k(t)$.

2. Substitute these expressions into the PDE and derive an ODE for each $a_k(t)$, including the Nyquist coefficient. Introduce $v_k(t)=a_k'(t)$ and write the resulting first-order system. Which two sets of initial Fourier coefficients are needed?

3. Explain how to compute those initial coefficients using the FFT of $f(x_j)$ and $g(x_j)$. State your normalization convention and frequency ordering. For MATLAB's FFT convention, normalized coefficients are `fft(values)/N`, and reconstruction is `ifft(N*coefficients)`. Alternatively, you may evolve the unnormalized FFT arrays throughout, provided the forcing and inverse transform use the same convention.

4. Arrange all displacement coefficients and all velocity coefficients into a single vector, retaining FFT order within each group. Implement its right-hand side and use your RK4 function to advance the state. Multiplication by the second-derivative multipliers can be done entry by entry; you do not need to form a large matrix.

5. Write a solver that takes $f$, $g$, $h$, $c$, $\alpha$, $N$, a final time $T$, and a time step, and returns the approximate displacement and velocity on the grid at requested output times. Note that the forcing must be evaluated at the correct time of **each RK4 stage**, including both midpoint stages. Reconstruct displacement and velocity with inverse FFTs. Check the size of the imaginary parts before reporting real-valued results.

**Note:** The spatial interval still has length $2\pi$, so the second-derivative multiplier is $-k^2$. No extension of the data is needed. Some FFT implementations label the Nyquist entry with $-m$; its second derivative uses the same multiplier $-m^2$ with either sign.

## Part 3. Verifying the PDE solver with an exact solution

We can test a forced PDE solver by choosing a solution first and then determining the data that make it satisfy the equation. This is called the **method of manufactured solutions**.

Let $c=1$, $\alpha=0.5$, and choose

$$
u(x,t)=\sin(3x)s(t),\qquad s(t)=e^{-t}\sin(2t).
$$

1. Determine $f(x)$, $g(x)$, and $h(x,t)$ so that this is the exact solution. Show your derivatives and verify the PDE and initial conditions.

2. Use $N=32$ and compute the numerical solution at $T=0.5$ with $\Delta t=1/100,\ 1/200,\ 1/400,$ and $1/800$. Compute

$$
E_{N,\Delta t}(T)=\max_{0\leq j \lt N}|p_N(x_j,T)-u(x_j,T)|.
$$

   Plot the error against $\Delta t$ and estimate the order as in Part 1. Plot one numerical solution together with the exact solution.

3. There are two potential sources of approximation error in this method: spatial discretization, from representing the solution with finitely many Fourier modes, and time stepping, from solving the coefficient ODEs numerically. Explain why $N=32$ can represent the spatial dependence of this exact solution without approximation error. Which source of error should therefore dominate before roundoff becomes significant?

4. Fix $\Delta t=1/200$ and repeat the calculation with $N=16,\ 32,\ 64,$ and $128$. Does increasing $N$ substantially improve the answer? Explain why a spatial refinement study alone cannot verify the accuracy of this solver.

## Part 4. Exploring resonance and damping

We will use numerical experiments to investigate two questions: how does damping change a wave's motion, and what happens when a periodic force repeatedly pushes the wave?

Set $c=1$ and use $N=64$ and $\Delta t=0.005$ throughout this part. For a single spatial mode $\sin(kx)$, the solution has the form

$$
u(x,t)=A_k(t)\sin(kx).
$$

The signed amplitude $A_k(t)$ describes the motion: a negative value means that the displacement has reversed sign. For $k=1$, record the solution at $x=\pi/2$ to obtain $A_1(t)$ directly.

You may use the following facts without deriving them. An unforced, undamped mode has **natural angular frequency** $c|k|$. If that mode is initially displaced and released from rest, its **critical damping coefficient** is $\alpha_{\mathrm{crit}}=2c|k|$. Below the critical value, positive damping produces oscillations that decrease in amplitude (**underdamped**). At the critical value, the displacement returns to zero as quickly as possible without oscillating (**critically damped**). Above it, the return is slower and still does not oscillate (**overdamped**).

1. **Displace the wave and let it go.** Set $f(x)=\sin(x)$, $g(x)=0$, and $h(x,t)=0$. Compute the solution for $0\leq t\leq8$ with $\alpha=0,\ 1,\ 2,$ and $4$. Plot the four amplitude histories on the same axes, with a horizontal line at zero. Identify the undamped, underdamped, critically damped, and overdamped cases. Which curves cross zero?

2. **Push the wave at different frequencies.** Now start from rest, with $f(x)=0$ and $g(x)=0$, and apply the periodic forcing

$$
h(x,t)=\sin(x)\cos(\omega t).
$$

   The factor $\sin(x)$ determines the spatial pattern of the force and excites only mode $k=1$. With $c=1$, this mode has natural angular frequency $c|k|=1$. The factor $\cos(\omega t)$ determines how rapidly the force oscillates in time. We choose $\omega=0.5,\ 1,$ and $1.5$ to push this same spatial mode below, at, and above its natural frequency.

   Set $\alpha=0$ and compute the solution for $0\leq t\leq40$ for each of these forcing frequencies. Plot the amplitude histories in separate panels with the same vertical scale. Which forcing frequency produces a response that keeps growing? This buildup when repeated forcing matches a natural frequency is called **resonance**. Describe what you observe in the other two cases as well.

3. **Add damping to the resonant experiment.** Keep the same forcing with $\omega=1$ and repeat the calculation with $\alpha=0.2,\ 0.5,$ and $2$. Compare these amplitude histories with the undamped result from the previous experiment. Does the response continue growing, or does it approach an oscillation of roughly constant amplitude?

4. **Start with many frequencies, but force only one.** Use

$$
f(x)=\frac{\sin(x)}{\frac54-\cos(x)},\qquad g(x)=0,
\qquad h(x,t)=\sin(4x)\cos(4t).
$$

   This smooth periodic initial displacement has infinitely many nonzero coefficients, while the force acts only on the $k=4$ mode at its natural frequency. Set $\alpha=0$ and compute the solution for $0\leq t\leq40$. Plot the initial displacement and several later snapshots. Which part of the solution grows? Does forcing one spatial frequency cause the other frequencies to grow as well?

   In particular, for your analysis, use spatial frequencies $1$, $4$, and $8$. Use the normalized Fourier coefficients from your solver to track these modes separately. For these sine-wave data, the signed amplitude of mode $k$ is

$$
A_k(t)=-2\operatorname{Im}(a_k(t)).
$$

   Here $a_k$ means the coefficient for signed frequency $k$, not array entry $k$. Plot $A_1(t)$, $A_4(t)$, and $A_8(t)$ in separate panels, each with its own vertical scale so that the smaller amplitudes remain visible.

5. **Choose enough damping to meet a displacement limit.** Keep the initial data and forcing from experiment 4. Suppose a vibrating structure must stay within an allowable displacement: no point may move more than $2$ units from equilibrium during $0\leq t\leq60$. We will use this as a simplified failure criterion. Unlike measuring mode 4 alone, this requirement includes all the spatial frequencies and how they combine. Find the minimum damping coefficient that meets this limit over the simulation.

   For each trial damping coefficient $\alpha$, run your PDE solver to $T=60$ and compute

$$
M(\alpha)=\max_{\substack{0\leq t_n\leq60\\0\leq j \lt N}}|p_N(x_j,t_n;\alpha)|,
\qquad F(\alpha)=M(\alpha)-2,
$$

   using every spatial grid point and every time step, including $t=0$. A positive value of $F$ means the displacement limit is exceeded; a nonpositive value means the limit is met on the sampled grid. Note that we can now frame the damping threshold as a root finding problem, $F(\alpha)=0$.

   Use the **bisection method**, either automatically or manually, by finding an initial bracket for the damping value. Verify that the lower endpoint exceeds the displacement limit and the upper endpoint meets it. At each midpoint, run the PDE solver to evaluate $F$ and update the bracket. Stop when the bracket width is less than $0.01$. Report the final bracket, with a failing lower endpoint and a passing upper endpoint, and choose its upper endpoint as your estimate of the minimum damping that meets the limit. Include a table of the trial damping coefficients and their values of $F$.

   Plot $\max_j|p_N(x_j,t;\alpha)|$ against time for your chosen damping coefficient, together with a horizontal line at $2$. As a check, extend the simulation to $T=80$ and evaluate the maximum over the entire interval $0\leq t_n\leq80$. Does your choice still meet the limit? What does this tell you about using a finite simulation to make a prediction about later motion? Compare your result with the critical damping coefficient for mode 4. Is critical damping necessary to meet the displacement requirement through $T=60$?

## Submission

Submit your code and a short report containing the derivations, requested plots, error studies, and explanations of your results. State the grid size, time step, and FFT normalization used in each experiment.

Your report should emphasize how Fourier differentiation and numerical time stepping work together, how periodic forcing excites individual modes, and how damping changes the response near resonance. Use the exact-solution and convergence studies to support the correctness of your implementation.
