#import "ieee.typ": *
#import "@preview/physica:0.8.1":bra, ket, braket, mel
#import "@preview/algorithmic:0.1.0"
#import algorithmic: algorithm
#import "@preview/lovelace:0.1.0": *

#let Ry = $op("R")_"y"$

#show: setup-lovelace

#show: ieee.with(
  title: [Variational Quantum Circuit Optimization via Interior Point Method],
  abstract: [Variational quantum algorithms (VQAs) have emerged as a promising approach to leverage near-term quantum devices and demonstrate quantum supremacy. However, highly expressive variational quantum circuits, known as ansatzes, often encounter trainability issues, such as the barren plateau problem, where gradients of trainable parameters vanish exponentially with increasing qubit numbers. To address these challenges, this project introduces specially designed circuit templates, focusing on the State Efficient Ansatz (SEA) that utilizes Schmidt decomposition to reduce problem size and amplify gradients. In particular, a sub-case of SEA is studied, where a layer of $Ry$ rotations serves as the Schmidt coefficient layer. The expressibility of this ansatz is investigated by casting it as an optimization problem. Then the problem is solved using the primal-dual interior point method, alongside a comparison with the conventional stochastic gradient descent method. The results demonstrate the superior convergence of the interior point method. Additionally, the observed solution structures offer potential refinements to the problem formulation, opening avenues for future generalization. This study contributes to the advancement of variational quantum algorithms, offering insights into optimizing circuit design and overcoming trainability limitations.
  ],
  authors: (
    (
      name: "Jiaxin Huang",
      department: [Dept. of Computer Science],
      organization: [University of Hong Kong],
      email: "jiaxin.huang@connect.hku.hk"
    ),
  ),
  bibliography-file: "refs.bib",
)

= Introduction

Variational quantum algorithms hold the promise of demonstrating quantum supremacy on noisy intermediate-scale quantum (NISQ) devices in the near future @Preskill2012QuantumCA. However, it has been observed that highly expressive variational quantum circuits, known as ansatzes, can suffer from trainability issues @zhang2021toward. One such issue is the barren plateau phenomenon, which occurs when the gradients of the trainable parameters exponentially vanish as the number of qubits increases @mcclean2018barren.

To address these challenges, researchers have proposed specially designed circuit templates that aim to strike a balance between expressibility and trainability @PhysRevX.11.041011. Among these templates is the State Efficient Ansatz (SEA) @liu2022mitigating, which leverages Schmidt decomposition @nielsen2010quantum to reduce the problem size by a factor of two, thereby quadratically amplifying the gradients. The structure of SEA is illustrated in Figure [insert figure number]. Unlike general ansatzes that require unitary universality, SEA only assumes wavefunction universality for the Schmidt coefficient layer, denoted as $U_1$. However, an efficient method for constructing $U_1$ that satisfies these properties remains unclear.

In this project, we focus on studying a specific sub-case of SEA, where $U_1$ is instantiated as a layer of $Ry$ rotations, given by:
$ U_1(bold(theta))=times.circle.big_(i=1)^N Ry (theta_i), $
where $theta_i$ represents the rotation angles. Our objectives are twofold. First, we investigate the expressibility of this ansatz, exploring whether it can approximate any arbitrary wavefunctions. Second, given the critical importance of trainability in VQA research, we employ the primal-dual interior-point method to optimize this ansatz, comparing its efficiency with the conventional gradient descent method commonly used in the literature.

= Literature Review

The optimization of variational quantum circuits currently relies on classical optimizers, as there is no established theory for designing quantum optimizers. In particular, machine learning optimizers have gained significant popularity due to their success in the field and their similarity to variational quantum circuits. Various machine learning optimization techniques have been widely adopted, ranging from the vanilla stochastic gradient descent method to more advanced gradient-based methods such as RMSprop and Adam @graves2013generating @kingma2014adam.

However, it has recently been discovered that universal quantum ansatzes may encounter the barren plateau problem, which refers to the exponential vanishing of gradients for trainable parameters as the number of qubits, denoted as $N$, increases. This phenomenon poses a significant limitation on the scalability of gradient-based optimization methods @mcclean2018barren. Consequently, alternative approaches that do not rely on gradients, known as gradient-free methods, have gained attention and consideration in the optimization of variational quantum circuits.

In a recent study by Monroig et al. @bonet2023performance, four gradient-free optimization methods (SLSQP, COBYLA, CMA-ES, and SPSA) were compared across multiple systems of varying sizes. This work aimed to assess the effectiveness and performance of these methods in addressing the barren plateau problem. By evaluating their performance on different systems, the researchers sought to provide insights into the applicability and potential limitations of gradient-free optimization techniques in the context of variational quantum circuits.

In @diez2021quantum, Díez-Valle et al. also formulated the placement of quantum gates in a circuit template as a quadratic unconstrained binary optimization (QUBO) problem. Traditionally, the placement of quantum gates within a circuit has been approached using heuristic algorithms and expert knowledge.

Traditionally, the discrete placement of quantum gates within a circuit is approached using heuristic algorithms and expert knowledge. Only continuous parameters within the gates are optimized. Inspired by classical automated machine learning (AutoML, @he2021automl), Díez-Valle et al. proposed a novel perspective by casting the gate placement problem as a quadratic unconstrained binary optimization (QUBO) problem, leveraging the power of optimization techniques commonly used in combinatorial optimization @diez2021quantum.


= Problem formulation <sec:formulation>

Let $ket(phi)$ be a $2N$ qubit quantum state to be synthesized. It admits Schmidt decomposition
$ ket(phi) &= sum_(j=1)^(2^N) omega_j ket(a_j)ket(b_j)\
          &= (U times.circle V)sum_(j=1)^(2^N) omega_j ket(j)ket(j), $
for some unitaries $U, V$ and $omega_j >= 0$.

Let $ket(psi)$ be the output of the SEA circuit, with SC layer $U_1=times.circle.big_(i=1)^N Ry$ and some LBC layer $U_2, U_3$. Below we analyze the effect of SEA on $2N$ basis states $ket(0)^(times.circle 2N)$

We start with the SC layer. Let $theta in RR$ be a rotation angle. The $Ry(theta)$ rotation gate applied to the basis state $ket(0)$ produces the effect 
$ ket(0) -> Ry (theta) ket(0) 
  & = sin theta ket(0) + cos theta ket(1)\
  // & = vec(delim: "[", sin theta_i, cos theta_i)
 $
// In vector notation, $ ket(bold(theta)) = times.circle.big_(i=1)^N  $

// Let $gamma_i$ be the $i$-th entry of the following vector of length $2^N$:

Therefore, given a vector $bold(theta) in RR^N$ of $N$ rotation angles, the SC layer $U_1(bold(theta))$ applied to $N$ basis states $ket(0)^(times.circle N)$ gives
$ ket(0)^(times.circle N) -> ket(bold(theta)) 
  & := U_1(bold(theta))ket(0)^(times.circle N)\
  & = (times.circle.big_(i=1)^N Ry (theta_i)) ket(0)^(times.circle N)\
  & = times.circle.big_(i=1)^N (Ry (theta_i) ket(0))\
  & = times.circle.big_(i=1)^N (sin theta_i ket(0)+cos theta_i ket(1))\
  & = sum_(i=1)^(2^N) mu_i ket(i),
 $
where $ bold(mu) := times.circle.big_(i=1)^N vec(delim: "[", sin theta_i, cos theta_i) in RR^(2^N). $
As an example, when $N=2$, $bold(mu)=vec(delim: "[", sin theta_1 sin theta_2, sin theta_1 cos theta_2, cos theta_1 sin theta_2, cos theta_1 cos theta_2).$

Then, the remaining layers give
$ ket(0)^(times.circle 2N) -> ket(psi) 
    & := op("SEA") ket(0)^(times.circle 2N)\
    // &= (U' times.circle V') op("CNOT")(times.circle.big_(i=1)^N op("R")_"y" (theta_i)ket(0),ket(0)^(times.circle N))\
    &= (U_2 times.circle U_3) op("CNOT")(ket(bold(theta)), ket(0)^(times.circle N))\
    &= (U_2 times.circle U_3) sum_(i=1)^(2^N) mu_i op("CNOT")( ket(i), ket(0)^(times.circle N))\
    &= (U_2 times.circle U_3) sum_(i=1)^(2^N) mu_i ket(i)ket(i).
 $

Therefore, the synthesis loss can be characterized by the fidelity 
$ F & = braket(psi, phi)\ 
  &= (sum_(i=1)^(2^N) mu_i bra(i)bra(i)) (U_2^dagger times.circle U_3^dagger)\
  & quad quad quad (U times.circle V)(sum_(j=1)^(2^N) omega_j ket(j)ket(j))\
  // & = 1 - 1/2 (sum_(i=1)^(2^N) mu_i bra(i)bra(i))(sum_(j=1)^(2^N) omega_j ket(j)ket(j))\
  & = sum_(i=1)^(2^N) sum_(j=1)^(2^N) mu_i omega_j mel(i, U_2^dagger U, j)mel(i, U_3^dagger V, j).\
  // & = 1 - 1/2 sum_(i=1)^(2^N) sum_(j=1)^(2^N) mu_i omega_j delta_(i,j)delta_(i,j)\
  // & = 1 - 1/2 sum_(i=1)^(2^N) mu_i omega_i.
 $
// Since $mel(i, U'U, j)mel(i, V'V, j) in [-1, 1]$, 
// $L$ is

This objective function involes optimization over unitaries $U_1$ and $U_2$, which is computationally difficult. We reduce the problem by assuming $U_2^dagger U=I$ and $U_3^dagger V=P$ for some permutation matrix $P$.

The rearrangement inequality states that, for any fixed pair of $bold(mu)$ and $bold(omega)$, $F$ is maximized when $P$ is a permutation that rearranges $bold(mu)$ to the same order as $bold(omega)$, i.e.,
$ omega_k < omega_l <=> (P bold(mu))_k < (P bold(mu))_l. $

Therefore, we could write
$ F = (P bold(mu))^T bold(omega). $


We further assume that $bold(mu)$ is homogenious across all of its $2^N$ dimensions. That is, if $bold(mu)^*$ is attainable, so is $P bold(mu)^*$ for all permutation $P$. (In fact this is not true due to the tensor product structure of $bold(mu)$.) Under this assumption, we can drop the matrix $P$ by assuming $bold(omega)$ is given in ascending order, i.e.,
$ F = bold(mu)^T bold(omega). $


// For any fixed $ket(phi)$, the performance of SEA can be characterized by the following optimization problem:
// $ limits("minimize")_(bold(theta)) & - F  \
//   "subject to" & 0 <= theta_i <= pi/2, forall i\
//   "where" & F = max_sigma sum_(i=1)^(2^N) mu_(i) omega_(sigma_i)\
//   & bold(mu)=times.circle.big_(i=1)^N vec(delim: "[", sin theta_i, cos theta_i)\
//     // & omega = op("sch-coef")(phi)
//  $

// An heuristic approximation of the above problem is

Moreover, we note that Schmidt decomposition guarantees $bold(omega) succ.eq 0$. Therefore, $F$ is maximized when $bold(u) succ.eq 0,$ which is equivalent to $0 prec.eq bold(theta) prec.eq pi/2 bold(I).$

In summary, synthesizing a given state $ket(phi)$ is equivalent to the following optimization problem:
$ limits("minmize")_(bold(theta)) & -bold(mu)^T bold(omega)\ "subject to" & 0 prec.eq bold(theta) prec.eq pi/2 bold(I)\ "where" & bold(mu)=times.circle.big_(i=1)^N vec(delim: "[", sin theta_i, cos theta_i)\ & omega_i <= omega_(i+1), forall i. $

// The performance of the SEA circuit can be characterized by average loss 
// $ limits(EE)_(omega) min_(bold(theta)) - F $
// and worse-case loss
// $ limits("maximize")_(omega) min_(bold(theta)) - F $

// $ dif^2/(dif theta_k dif theta_l) F = - sum_(i=1)^(2^N) mu''_(i)omega_i,
//  $
// where 
// $ bold(mu)'':=times.circle.big_(i=1)^(k-1) vec(delim: "[", sin theta_i, cos theta_i) times.circle vec(delim: "[", sin theta_k, cos theta_k)
//  $
// $ nabla^2 f =  $

// $ bold(omega)  $

// We now demonstrate that the above problem is concave. Firstly, $sin$ and $cos$ functions are positive and concave in $(0,pi/2)$. Therefore, each $mu_i$, as a product of $sin(theta_i)$ and $cos(theta_j)$, is also concave. Hence, $sum_(i=1)^(2^N) mu_(i) omega_i$, as a non-negative combination of $mu_i$, is also concave. 

// We now determine the convexity of this optimization problem. We first note that $sum_(i=1)^(2^N) mu_(i)omega_i$ is a non-negative combination of $mu_i$, and $omega_i$ can be arbitrary. Therefore, the problem is convex (concave) if and only if $mu_i$ is convex (concave) for all $i$.


= Algorithm


The interior point method was chosen to solve the optimization problem. There is no need to implement phase 1 method because a feasible starting point can be easily found: currently, it's chosen as $bold(theta)^((0))=pi/4 bold(1)$, which is the center of the feasible set. To accelerate convergence, the primal-dual form of the algorithm is used. 

Several important functions are explained below. The objective function is an inner product of two vectors with length $2^N$:
$ f_0(bold(theta)) := - (times.circle.big_(i=1)^N vec(delim: "[", sin theta_i, cos theta_i))^T bold(omega). $
The constraints are 
$ bold(g)(bold(theta)) := op("concat")(-bold(theta),bold(theta)- pi / 2 bold(1)) prec.eq 0, $
which corresponds to $0 <= theta_i <= pi/2$ exactly.

The residual is
$ bold(r)(bold(theta), bold(lambda)) = vec(delim: "[", bold(r)_"dual", bold(r)_"cent") := vec(delim: "[", nabla f_0(bold(theta)) + D bold(g)(bold(theta))^T bold(lambda),-op("diag")(bold(lambda))bold(g)(bold(theta))-1/t bold(1)). $
Note that there is no primal residual as the optimization problem has no equality constraints.

In the primal-dual form, solving the optimization method is equivalent to find the root
$ bold(r)(bold(theta)^*, bold(lambda)^*)=0. $

By Newton method, this can be iteratively solved by solving the linear approximation
$ bold(r)(bold(theta)^((k)), bold(lambda)^((k))) + [D bold(r)(bold(theta)^((k)), bold(lambda)^((k)))](Delta bold(theta), Delta bold(lambda))^T=0, $
where $D bold(r)$ is the Jacobian
$ D bold(r)(bold(theta), bold(lambda)) = mat(delim: "[", nabla^2 f_0(bold(theta))+sum_(i=1)^(2N)lambda_i nabla^2 g_i(bold(theta)), D bold(g)(bold(theta))^T; -op("diag")(bold(lambda))D bold(g)(bold(theta)), -op("diag")(bold(g)(bold(theta)))). $

After finding the direction $(Delta bold(theta), Delta bold(lambda))$, the backtracking line search is used to determine the step size. Backtracking line search finds the larest step size $s$ which:
- $0<s<1$;
- maintains the primal constraint $bold(g)(bold(theta)) prec.eq 0$;
- maintains the dual constraint $bold(lambda) succ.eq 0$;
- reduces $||bold(r)||_2$ by at least a factor of $alpha s$ (currently, $alpha$ is set as 0.5).

Lastly, the algorithm terminates when dual residual is  smaller than a predetermined threshold $epsilon_"feas"$, and the surrogate duality gap 
$ eta := -bold(g)(bold(theta))^T bold(lambda) $
is smaller than another threshold $epsilon$.

The complete pseudocode is shown in @algo.

To evaluate the performance of the interior point method, a stochastic gradient descend (SGC) method is also implemented. The gradient descend method simply computates the gradient of the objective function
$ nabla f_0(bold(theta))=vec(delim: "[", - (vec(delim: "[", cos theta_1, -sin theta_1)times.circle(times.circle.big_(i=2)^N vec(delim: "[", sin theta_i, cos theta_i)))^T bold(omega), dots.v, - ((times.circle.big_(i=1)^(N-1) vec(delim: "[", sin theta_i, cos theta_i))  times.circle vec(delim: "[", cos theta_N, -sin theta_N))^T bold(omega)) $
and update parameters
$ bold(theta) <- bold(theta) + s nabla f_0(bold(theta)) $
by some step size $s$.


#figure(supplement: "Algorithm", caption: [Primal-dual interior point method])[
#set align(left)
#pseudocode(
  no-number,
  [*input:* $bold(omega) succ.eq 0, epsilon_"feas" >0, epsilon > 0, nu > 1,$],
  no-number,
  align(right)[$beta in (0,1), alpha in (0,1)$],
  no-number,
  [*output:* $bold(theta), f_0(bold(theta))$],
  $bold(lambda) <- bold(1)_(2N)$,
  [*repeat*], ind,
    $t <- 2 nu N \/ eta$,
    
    
    [*solve* $A vec(delim: "[", Delta bold(theta), Delta bold(lambda))+bold(r)=0$ for $vec(delim: "[", Delta bold(theta), Delta bold(lambda))$],

    $s <- min(0.99, -lambda_i\/Delta lambda_i : Delta lambda_i < 0)$,
    [*repeat*], ind,
      $s <-  beta s$, ded,
    [*until* $bold(g)(bold(theta) + s Delta bold(theta)) prec.eq 0, $ and $||bold(r)(bold(theta) + s Delta bold(theta), bold(lambda) + s Delta bold(lambda))||$],
    no-number,ded,
    align(right)[$<= (1-alpha s)||bold(r)(bold(theta), bold(lambda))||$],ind,
  $bold(theta) <- bold(theta) + s Delta bold(theta)$,
  $bold(lambda) <- bold(lambda) + s Delta bold(lambda)$,
  $bold(r)_"dual" <- bold(r)(bold(theta), bold(lambda))[:N]$,
  $eta <- -bold(g)(bold(theta))^T bold(lambda)$,ded,
  [*until* $||bold(r)_"dual"|| <= epsilon_"feas", eta <= epsilon$ ],
  [*return* $bold(theta), f_0(bold(theta))$]
)]<algo>

= Evaluation

#let subfigure = figure.with(caption: "", kind: "subfigure", supplement: none, numbering: "(a) ", outlined: false)

#show figure: it => {
  if it.kind == "subfigure" {
    align(bottom, it)
  } else {
    locate(loc => {
      let q = query(figure.where(kind: "subfigure").after(loc), loc)
      q.first().counter.update(0)
    })
    it
  }
}

#show ref: it => {
  if it.element != none and it.element.func() == figure and it.element.kind == "subfigure" {
    locate(loc => {
      let q = query(figure.where(outlined: true).before(it.target), loc).last()
      ref(q.label)
    })
  }
  it
}

We evaluate the optimization methods by using randomly sampled $2N$-qubit quantum state $ket(phi)$. The coefficient $bold(omega)$ can be obtained by singular value decomposition
$ M_(ket(phi)) = U op("diag")(bold(omega)) V, $
where $M_(ket(phi))$ is the inverse vectorziation of $ket(phi)$ by rewriting the $2^(2N)$ dimentional vector as a $2N times 2N$ matrix.

We tested two scenarios where $N=5$ and $N=10$. For the interior point optimizer, the stopping criteria is when $bold(r)_"dual" <= epsilon_"feas"$ and $eta <= epsilon$. We used $epsilon_"feas"=epsilon=10^(-2)$ for both $N=5$ and $N=10$ cases. For the stochastic gradient descend optimizer, there is no explicit stopping criteria. So we terminat it after the same iteration number of the interior point optimizer on the same problem. In experiments, we observe that the typical iteration number is 51 for $N=5$ and 97 for $N=10$.

The optimization results are shown in @fig:wave and @fig:f. In @fig:wave. we observe that the wavefunction generated by the interior point method (orange) exhibits a closer approximation to $bold(omega)$ (green), which agrees with the fact that the objective function decreases faster for interior point method in @fig:f. Moreover, in both @fig:f5 and @fig:f10, the curve for interior point method exhibits larger fluctuation than the curve for stochastic gradient descend method. This can be explained by the fact that the interior point method is not a descend method. It steadily descreases the residual $bold(r)$, but not necessarily the objective function $-F$.

Finally, we note that for large quantum systems, for example $N=10$ in @fig:wave10, the optimal value exhibits a very clear pattern: the $bold(omega)$ is almost a straight line from 0 to $sqrt(3/(2^N))$, and $bold(mu)$ is roughly four disjoint  line sengments. Such clear pattern may allow us to drop some assumptions in @sec:formulation, hence further improving synthesis accuracy.
 

#figure(
    stack(spacing: 1em,
      [#subfigure(image("../figs/5 wave.png", width:50%), caption: [$2N=10$ quibits])<fig:wave5>],
      [#subfigure(image("../figs/10 wave.png", width:50%), caption: [$2N=20$ quibits])<fig:wave10>],
    ),
    caption: [The wavefunction $mu$ and $omega$ after optimization.]
  ) <fig:wave>

#figure(
    stack(spacing: 1em,
      [#subfigure(image("../figs/5 f.png", width:50%), caption: [$2N=10$ quibits])<fig:f5>],
      [#subfigure(image("../figs/10 f.png", width:50%), caption: [$2N=20$ quibits])<fig:f10>],
    ),
    caption: "Fidelity over optimization iterations."
  ) <fig:f>
// #figure(grid(columns: 2, row-gutter: 2mm, column-gutter: 1mm,
//   image("../figs/5 f.png", width:50%),
//   image("../figs/5 f.png", width:50%), 

//   "a) test 1",
//   "b) test2"
//   ),
//   caption: "Caption"
// )<label>

// #figure(#box())

= Discussion
In this project, the formulation of a variational quantum algorithm as an optimization problem has been investigated. The problem was tackled using two optimization methods, namely the primal-dual interior point method and the vanilla stochastic gradient descent method. Through comprehensive experiments, it has been demonstrated that the interior point method exhibits faster convergence compared to the stochastic gradient descent method. This result sheds light on the potential for developing new optimizers that do not solely rely on the gradients of the objective function. By exploring alternative optimization methodologies, researchers may overcome the limitations imposed by the barren plateau problem in large variational circuits and pave the way for more efficient and effective variational quantum algorithms.

However, it is important to consider the computational efficiency implications associated with the interior point method. This approach necessitates the calculation of Hessian matrices for the objective and constraint functions, as well as solving a linear system with a dimension of $3N$ at each iteration. As a consequence, when optimizing large quantum circuits with a substantial number of tunable parameters, careful consideration must be given to strike a balance between computational efficiency and optimization effectiveness.

Moreover, the solution structures observed throughout this project offer valuable insights into potential refinements of the assumptions employed in the problem formulation. By refining these assumptions, it may be possible to improve computational efficiency or generalize the problem formulation, thereby expanding its applicability to a wider range of scenarios.

In summary, this project has successfully addressed the optimization of variational quantum algorithms by formulating them as optimization problems and exploring different solution methodologies. The findings highlight the superiority of the interior point method in terms of convergence speed and provide inspiration for the development of novel optimizers that can circumvent the barren plateau problem. Additionally, the computational efficiency considerations and the potential for refining the problem formulation assumptions underscore the need for continued research in this field to unlock the full potential of variational quantum algorithms.