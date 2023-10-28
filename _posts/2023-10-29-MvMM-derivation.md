---
layout: post
title:  "MvMM: Multi-variate Mixture Model"
date:   2023-10-29 10:20:00
blurb: "MvMM: Multi-variate Mixture Model"
# og_image: /assets/img/content/StudentT/gamma_distribution.svg
---


# MvMM derivation

Yuanye Liu 2023.10


## Introduction

This is a simple but detailed derivation for MvMM [[Link]](https://ieeexplore.ieee.org/document/8458220). The prerequisite is GMM-EM which you can refer to [[GMM]](https://henrylau7.github.io/2023/10/28/GMM). The PDF version is in my Github [Notes Page](https://github.com/HenryLau7/Notes/blob/main/2023-10-29-MvMM-derivation.pdf), where you can give me a precious star if you think it's useful for you and I'll appreciate it :)

## Contents
1. Likelihood function [[Link]](#likelihood-function)
2. EM Solution[[Link]](#em-solution)
3. Spatial Constraints. [[Link]](#spatial-regularization--constraint)
4. Registration in MvMM. [[Link]](#registration-in-mvmm)
5. HCMMI [[Link]](#hetero-coverage-multi-modality-images-hcmmi)


## Notation
* Common space: $\Omega$
    * $x\in\Omega$: pixel
    * $I(x) / I_i(x)$: intensity
* Image / Sequence / Modality:
    * $i\in N_I$
    * $N_I=3$ in this task: LGE, T2, bSSFP
* Type: defined on common space
    * $k\in K$
    * Latent variable: $s(x)=k$
    * Prior: $\pi_{kx}$ or $\pi_k$
* Subtype: defined on a specific image and its type
    * $c\in C_{ik}$ (NOTE: related to image and type)
    * Latent variable: $z_i(x)=c$ or $c_{ik}$
    * Prior: $\tau_{ikc}$
    * Gaussian parameters: $\mu_{ikc}, \sigma^2_{ikc}$


## Likelihood function
First we derive the incomplete likelihood function for $\theta=\{\pi,\tau,\mu,\sigma\}$

$$
\begin{aligned}
L(\theta\mid I)&= \prod_{x\in\Omega}\sum_{k\in K}\pi_{kx}\prod_{i\in N_I} \sum_{c\in C_{ik}} \tau_{ikc} \Phi(I_i(x);\mu_{ikc},\sigma^2_{ikc} )\\
:&=P(I\mid\theta)
\\
&=\prod_{x\in\Omega}P(I(x)\mid \theta)
\\
&=\prod_{x\in\Omega}\sum_{k\in K}P(I(x)\mid s(x)=k,\theta)P(s(x)=k\mid\theta)
\\
&=\prod_{x\in\Omega}\sum_{k\in K}\pi_{kx}\prod_{i\in N_I}P(I_i(x)\mid s(x)=k, \theta)
\\
&=\prod_{x\in\Omega}\sum_{k\in K}\pi_{kx}\prod_{i\in N_I} \sum_{c\in C_{ik}}P(z_i(x)=c\mid s(x)=k,\theta)P(I_i(x)\mid z_i(x)=c,\theta)
\\
&= \prod_{x\in\Omega}\sum_{k\in K}\pi_{kx}\prod_{i\in N_I} \sum_{c\in C_{ik}} \tau_{ikc} \Phi(I_i(x);\mu_{ikc},\sigma^2_{ikc} )
\end{aligned}
$$


## EM solution
Recall the EM solution, we need the complete log-likelihood function in E-step, so we now derive the complete likelihood function.

### Complete likelihood
Here we use the trick, introduce the indicator function $1\{(s(x)=k)\}$ and $1\{(z_i(x)=c)\}$.

$$
\begin{aligned}
L(\theta\mid I,S,Z)
:
&=P(I,S,Z\mid\theta)
\\
&=\prod_{x\in\Omega}\prod_{k\in K}\{\pi_{kx}\prod_{i\in N_I} \prod_{c\in C_{ik}} \{ \tau_{ikc} \Phi(I_i(x);\mu_{ikc},\sigma^2_{ikc} )\}^{1\{z_i(x)=c\}}\}^{1\{s(x)=k\}}\\
\end{aligned}
$$

### Complete log-likelihood
Take log on both sides.

$$
\begin{aligned}
l(\theta\mid I,S,Z)

&=\log L(\theta\mid I,S,Z)
\\
&=\sum_{x\in\Omega}\sum_{k\in K}1_{(s(x)=k)}\log\pi_{kx}\sum_{i\in N_I} \sum_{c\in C_{ik}} 1_{(z_i(x)=c)}\{\log\tau_{ikc}+ \log\Phi(I_i(x);\mu_{ikc},\sigma^2_{ikc} )\}
\\
&=\sum_{x\in\Omega}\sum_{k\in K}1_{(s(x)=k)}\log\pi_{kx} + \sum_{x\in\Omega}\sum_{k\in K}\sum_{i\in N_I} \sum_{c\in C_{ik}}1_{(s(x)=k)}1_{(z_i(x)=c)}\{\log\tau_{ikc}+ \log\Phi(I_i(x);\mu_{ikc},\sigma^2_{ikc} )\}

\end{aligned}
$$


### E-step: $Q$ function
In E-step, we are going to take conditioanl expectation for complete log-likelihood function we have derived above, with respect to latent varaible given obervation $I(x)$ and current pamameter $\theta^{[m]}$, i.e. $S,Z\mid I(x),\theta^{[m]}$. 

So the Q-function is,

$$
\begin{aligned}
Q(\theta\mid \theta^{[m]})
&=\mathbb{E}_{S,Z|\theta^{[m]}}l(\theta\mid I,S,Z)
\\
&=\sum_{x\in\Omega}\sum_{k\in K}\mathbb{E}_{S,Z|\theta^{[m]}}\{1_{(s(x)=k)}\}\log\pi_{kx} 
\\
&+ \sum_{x\in\Omega}\sum_{k\in K}\sum_{i\in N_I} \sum_{c\in C_{ik}}\mathbb{E}_{S,Z|\theta^{[m]}}\{1_{(s(x)=k)}1_{(z_i(x)=c)}\}\{\log\tau_{ikc}+ \log\Phi(I_i(x);\mu_{ikc},\sigma^2_{ikc} )\}
\end{aligned}
$$


### E-step: Compute posterior function
We can find that it is the indicator functions that are related to the expected variable, which are the posterior,

$$
\begin{aligned}
\mathbb{E}_{S,Z\mid I,\theta^{[m]}}\{1_{(s(x)=k)}\}
&=P(s(x)=k\mid I,\theta^{[m]} )
\\
&=\frac{P(I(x)\mid s(x)=k,\theta^{[m]})P(s(x)=k\mid\theta^{[m]})}{\sum_{l\in K}P(I(x)\mid s(x)=l,\theta^{[m]})P(s(x)=l\mid\theta^{[m]})}
\\
&=\frac{P(I(x)\mid s(x)=k,\theta^{[m]})\pi_{kx}^{[m]}}{\sum_{l\in K}P(I(x)\mid s(x)=l,\theta^{[m]})\pi_{lx}^{[m]}}
\\
&:=P_{kx}^{[m+1]}

\end{aligned}
$$

And then

$$
\begin{aligned}
&\mathbb{E}_{S,Z\mid I,\theta^{[m]}}\{1_{(s(x)=k)}1_{(z_i(x)=c)}\}
\\
&=P(s(x)=k,z_{i}(x)=c_{ik}\mid I,\theta^{[m]} )
\\
&=P(z_i(x)=c_{ik}\mid s(x)=k,I,\theta^{[m]})P(s(x)=k\mid I,\theta^{[m]})
\\
&=P(z_i(x)=c_{ik}\mid s(x)=k,I,\theta^{[m]})P_{kx}^{[m+1]}
\\
&=\frac{P(z_i(x)=c_{ik}\mid s(x)=k,\theta^{[m]})P(I(x)\mid z_i(x)=c_{ik},s(x)=k,\theta^{[m]})}{P(I(x)\mid s(x)=k,\theta^{[m]})} P_{kx}^{[m+1]}
\\
&=\frac{\tau_{ikc}\Phi(I_i(x)\mid\mu_{ikc}^{[m]},\sigma_{ikc}^{2{[m]}})}{P(I_i(x)\mid s(x)=k,\theta^{[m]})} P_{kx}^{[m+1]}
\\
&:=P_{ikcx}^{[m+1]}

\end{aligned}
$$

So the $Q$ function is,

$$
\begin{aligned}
Q(\theta\mid \theta^{[m]})

&=\sum_{x\in\Omega}\sum_{k\in K}\mathbb{E}_{S,Z|I,\theta^{[m]}}\{1_{(s(x)=k)}\}\log\pi_{kx} 
\\
&+ \sum_{x\in\Omega}\sum_{k\in K}\sum_{i\in N_I} \sum_{c\in C_{ik}}\mathbb{E}_{S,Z|I,\theta^{[m]}}\{1_{(s(x)=k)}1_{(z_i(x)=c)}\}\{\log\tau_{ikc}+ \log\Phi(I_i(x);\mu_{ikc},\sigma^2_{ikc} )\}
\\
&=\sum_{x\in\Omega}\sum_{k\in K}P_{kx}^{[m+1]}\log\pi_{kx} 
\\
&+ \sum_{x\in\Omega}\sum_{k\in K}\sum_{i\in N_I} \sum_{c\in C_{ik}}P_{ikcx}^{[m+1]}\{\log\tau_{ikc}+ \log\Phi(I_i(x);\mu_{ikc},\sigma^2_{ikc} )\}
\end{aligned}
$$

### M-Step: Maximize the $Q$ function 
In the M-step, we want to maximize $Q(\theta\mid \theta^{[m]})$, with respect to $\theta=(\pi_{kx},\tau_{ikc},\mu_{ikc},\sigma^2_{ikc})$. We can obtain the closed-form solution by letting the derivative equal to zero.

#### for $\tau$

$$
\tau_{ikc}^{[m+1]}=\arg \max_{\tau_{ikc}}Q(\theta\mid\theta^{[m]})
\\
s.t. \sum_{d\in C_{ik}} \tau_{ikd}=1
$$

Construct Lagrange multiplier

$$
L(\theta,\lambda)=Q(\theta\mid \theta^{[m]})+\lambda(1-\sum_{d\in C_{ik}}\tau_{ikd})
\\
\begin{aligned}
\frac{\partial L}{\partial \tau_{ikc}}

&=\sum_{x\in\Omega}\frac{P_{ikcx}^{[m+1]}}{\tau_{ikc}}+\lambda=0\\

\sum_{d\in C_{ik}}\tau_{ikd}&=1

\end{aligned}
$$

Solution:

$$
\tau_{ikc}=\frac{\sum_{x\in \Omega}P_{ikcx}^{[m+1]}}{\sum_{d\in C_{ik}}\sum_{x\in \Omega}P_{ikdx}^{[m+1]}}
$$

#### For $\pi_{kx}$
when no spatial constraint is applied, $\pi_{kx}=\pi_{k}$.

$$
\pi_{k}^{[m+1]}=\arg \max_{\pi_{kx}} Q(\theta\mid \theta^{[m]})
\\
s.t. \sum_{k\in K}  \pi_{k} =1
$$

Solution

$$
\pi_{k}^{[m+1]}=\frac{\sum_{x\in \Omega}P_{kx}^{[m+1]}}{\sum_{l\in K}\sum_{x\in\Omega}P_{kx}^{[m+1]}}
$$


#### For $\mu_{ikc}$

$$
\mu_{ikc}^{[m+1]}=\arg\max_{\mu_{ikc}} Q(\theta\mid\theta^{[m]})
$$

Solution

$$
\frac{\partial Q(\theta\mid\theta^{[m]})}{\partial \mu_{ikc}}=\sum_{x\in\Omega}P_{ikcx}^{[m+1]}(I_i(x)-\mu_{ikc})=0\\

\mu_{ikc}^{[m+1]}=\frac{\sum_{x\in \Omega }P_{ikcx}^{[m+1]}I_i(x)}{\sum_{x\in \Omega}P_{ikcx}^{[m+1]}}
$$


#### For $\sigma^{2[m+1]}_{ikc}$


$$
\sigma^{2[m+1]}_{ikc}=\arg \max_{\sigma^{2}_{ikc}} Q(\theta\mid\theta^{[m]})
$$


$$
\frac{\partial Q(\theta\mid\theta^{[m]})}{\partial \sigma^{2}_{ikc}}=\sum_{x\in\Omega}P_{ikcx}^{[m+1]}(-\frac{1}{2}\frac{1}{\sigma^2_{ikc}}+\frac{(I_i(x)-\mu_{ikc}^{[m+1]})^2}{2 (\sigma^{2}_{ikc})^2})=0

\\

\sigma^{2[m+1]}_{ikc}=\frac{\sum_{x\in \Omega }P_{ikcx}^{[m+1]}(I_i(x)-\mu_{ikc}^{[m+1]})^2}{\sum_{x\in \Omega}P_{ikcx}^{[m+1]}}
$$

We can note that, here we use the newest $\mu^{[m+1]}$ for the update of $\sigma^{2[m+1]}$

## Spatial regularization / constraint
The motivation of spatial constrain is, pixels with the same intensity distribution in medical images can come from different structures.

### Probabilistic Atlas
Probabilistic atlas is a common used method for spatial constraint [[Link]](https://www.sciencedirect.com/science/article/pii/S1361841504000271).

After the introducing of atlas, the prior now is,

$$
\pi_{kx}=P(s(x)=k\mid\theta)=\frac{\pi_kP_A(s(x)=k)}{NF}=\frac{\pi_kP_A(s(x)=k)}{\sum_{l\in K}\pi_l P_A(s(x)=l)}
$$

where $NF$ is normalization factor. $NF=\sum_{l\in K}\pi_l P_A(s(x)=l)$.

Similarly, we want to maximize the Q-function, we denote $Q_{\pi}$ referred to the $\pi$ part.


$$
\begin{aligned}

Q_{\pi}
&=\sum_{x\in\Omega}\sum_{k\in K} P_{kx}^{[m+1]}\log \pi_{kx}=\sum_{x\in\Omega}\sum_{k\in K} P_{kx}^{[m+1]}\log \frac{\pi_kP_A(s(x)=k)}{\sum_{l\in K}\pi_l P_A(s(x)=l)}
\\
&=\sum_{x\in\Omega}\sum_{k\in K} P_{kx}^{[m+1]}\log{\pi_k}

-\sum_{x\in\Omega}\sum_{k\in K} P_{kx}^{[m+1]}\log \frac{\sum_{l\in K}\pi_l P_A(s(x)=l)}{P_A(s(x)=k)}
\\
&=\sum_{x\in\Omega}\sum_{k\in K} P_{kx}^{[m+1]}\log{\pi_k}

-\sum_{x\in\Omega}\sum_{k\in K} P_{kx}^{[m+1]}
\{\log {\sum_{l\in K}\pi_l P_A(s(x)=l)}-\log\sum_{k\in K}{P_A(s(x)=k)}\}
\end{aligned}
$$

Then,

$$
\begin{aligned}
\frac{\partial Q_{\pi}}{\partial \pi_{k}}

&=\sum_{x\in \Omega}\frac{P_{kx}^{[m+1]}}{\pi_k}-\sum_{x\in\Omega}\sum_{j\in K}\frac{P_{jx}^{[m+1]}P_A(s(x)=k)}{\sum_{l\in K}\pi_l P_A(s(x)=l)}

\\

&=\sum_{x\in \Omega}\frac{P_{kx}^{[m+1]}}{\pi_{k}}-\sum_{x\in \Omega }\frac{P_A(s(x)=k)\sum_{j\in K}P_{jx}^{[m+1]}}{\sum_{l\in K}\pi_l P_A(s(x)=l)}
\\
&=\sum_{x\in \Omega}\frac{P_{kx}^{[m+1]}}{\pi_{k}}-\sum_{x\in \Omega }\frac{P_A(s(x)=k)}{\sum_{l\in K}\pi_l P_A(s(x)=l)}

=0
\end{aligned}
$$

Since there is no closed-form solution. So we assume on the left hand side, $\pi_k=\pi_{k}^{[m+1]}$, the right hand side, $\pi_l=\pi_l^{[m]}$. And we denote the constant ${C_x^{[m]}=\sum_{l\in K}\pi_l ^{[m]}P_A(s(x)=l)}$. So then we have a iterative ${\pi_k^{[m+1]}}$,

$$
\pi_{k}^{[m+1]}=\frac{\sum_{x\in \Omega}P_{kx}^{[m+1]}}{\sum_{x\in\Omega}(P_A(s(x)=k)/C_x^{[m]})}
$$

We need to prove that $\pi_k^{[m+1]}$ can increase the likelihood,

$$
Q_{\pi}^{[m]}=\sum_x\sum_k P_{kx}^{[m+1]}\log \frac{\pi_k^{[m]}P(A_{kx})}{C_x^{[m]}}
\\
Q_{\pi}^{[m+1]}=\sum_x\sum_k P_{kx}^{[m+1]}\log \frac{\pi_k^{[m+1]}P(A_{kx})}{C_x^{[m+1]}}
\\
Q_{\pi}^{[m+1]} - Q_{\pi}^{[m]}=\sum_x\sum_k P_{kx}^{[m+1]}\log \frac{\pi_k^{[m+1]}}{\pi_k^{[m]}}\frac{C_x^{[m]}}{C_x^{[m+1]}}
$$

We see

$$
\begin{aligned}

\frac{\pi_k^{[m+1]}}{\pi_k^{[m]}}\frac{C_x^{[m]}}{C_x^{[m+1]}}
&=\frac{\sum_{y\in\Omega}P_{ky}^{[m+1]}}{\sum_{y\in\Omega}(P_A(s(y)=k)/C_y^{[m]})}\frac{1}{\pi_k^{[m]}}\frac{C_x^{[m]}}{C_x^{[m+1]}}
\\
&=\frac{\sum_{y\in\Omega}P_{ky}^{[m+1]}}
{
{\sum_y \pi_k^{[m]}{P_A(s(y)=k)} / {C_y^{[m]}}}

}

\frac{C_x^{[m]}}{C_x^{[m+1]}}

\\
&=\frac{\sum_{y\in\Omega}P_{ky}^{[m+1]}}

{\sum_y \pi_{ky}^{[m]} }
\frac{C_x^{[m]}}{C_x^{[m+1]}}

=\frac{\sum_{y\in\Omega}P_{ky}^{[m+1]} / {C_x^{[m+1]}}}

{\sum_y \pi_{ky}^{[m]} /{C_x^{[m]}}}
\end{aligned}
$$

So we have,


$$
\begin{aligned}
Q_{\pi}^{[m+1]} - Q_{\pi}^{[m]}

&=\sum_x\sum_k P_{kx}^{[m+1]}\log \frac{\sum_{y\in\Omega}P_{ky}^{[m+1]} / {C_x^{[m+1]}}}

{\sum_y \pi_{ky}^{[m]} /{C_x^{[m]}}}

\\
&=  {C_x^{[m+1]}} \sum_k \sum_x P_{kx}^{[m+1]}/ {C_x^{[m+1]}} \log \frac{\sum_{y\in\Omega}P_{ky}^{[m+1]} / {C_x^{[m+1]}}}

{\sum_y \pi_{ky}^{[m]} /{C_x^{[m]}}}

\\

&=  {C_x^{[m+1]}} |\Omega|\sum_k \sum_x P_{kx}^{[m+1]}/ {(C_x^{[m+1]}|\Omega|)} \log \frac{\sum_{y\in\Omega}P_{ky}^{[m+1]} / {(C_x^{[m+1]}|\Omega|)}}

{\sum_y \pi_{ky}^{[m]} /{(C_x^{[m]}|\Omega|)}}

\\

&= {C_x^{[m+1]}} |\Omega|\cdot KL(P\Vert \pi) \ge 0

\end{aligned}
$$

So, we proved that $Q_{\pi}$ increases with the update of $\pi_{k}$. It is the GEM (General EM) algorithm.

### Initialization
We start from $\theta^{[0]}$, and then compute the posterior $P_{ik}^{[1]}$ and $P_{ickx}^{[1]}$, and then iteratively.

For $\pi_k^{[0]}$:

$$
\pi_k^{[0]}=\frac{\sum_{x}P(A_{kx})}{\sum_l\sum_x P(A_{lx})}
$$

For $\tau_{ikc}^{[0]}$:

$$
\tau_{ikc}^{[0]}=\frac{1}{|C_{ik}|}
$$

For $\mu_{ikc}^{[0]},\sigma^{2[0]}_{ikc}$:

$$
\mu_{i k c}^{[0]}=\left\{\begin{array}{ll}
\mu_{i k}^{[0]}+a \sigma_{i k}^{[0]}, & \left|C_{i k} \geq 2\right| \\
\mu_{i k}^{[0]}, & \left|C_{i k}=1\right|
\end{array},\left(\sigma_{i k c}^{[0]}\right)^{2}=\frac{\left(\sigma_{i k}^{[0]}\right)^{2}}{\left|C_{i k}\right|}\right.
$$

where

$$
\mu_{i k}^{[0]}=\frac{\sum_{x} I_{i}(x) p\left(A_{k x}\right)}{\sum_{x} p\left(A_{k x}\right)}, \text { and }\left(\sigma_{i k}^{[0]}\right)^{2}=\frac{\sum_{x}\left(I_{i}(x)-\mu_{i k}^{[0]}\right)^{2} p\left(A_{k x}\right)}{\sum_{x} p\left(A_{k x}\right)}
$$

## Registration in MvMM
Two types of misalignment in Myocardial Segmentation:
* Inter-slice: motion shift.
* Misalignment between atlas and common space.

For the slice reg, we model it as a rigid transformation

$$
P(I_i(x)\mid c_{ik},\theta,G_{i,s})=\Phi_{ikc}(I_i(G_{i,s}(x)))
$$

where $G_{i,s}$ the transformations for correcting.

The atlas deformation, modeled as a FFD, denoted as $D$,

$$
P_A(S(x)=k\mid D)=P_A(s(D(x))=k)=A_k(D(x))
$$

So the prior is $\pi_{kx\mid D}=P(s(x)=k\mid D)$.

The likelihood:

$$
\begin{aligned}
L(\theta,D,\{G_{i,s}\})
&=\prod_{x\in \Omega}\sum_{k\in K} \pi_{kx\mid D}P(I(x)\mid s(x)=k,\theta,\{G_{i,s}\})
\\
&=\prod_{x\in \Omega}\sum_{k\in K} \pi_{kx\mid D}\prod_{i\in N_I}P(I_i(x)\mid s(x)=k,\theta,G_{i,s})
\\
&=\prod_{x\in \Omega}\sum_{k\in K} \pi_{kx\mid D}\prod_{i\in N_I}\sum_{c_{ik}}\tau_{ikc}\Phi_{ikc}(I_i(G_{i,s}(x)))
\end{aligned}
$$

So the log-likelihood:

$$
LL(\theta,D,\{G_{i,s}\})=\log L(\theta,D,\{G_{i,s}\})=\sum_{x\in \Omega}\log\{\sum_{k\in K} \pi_{kx\mid D}\prod_{i\in N_I}\sum_{c_{ik}}\tau_{ikc}\Phi_{ikc}(I_i(G_{i,s}(x)))\}
$$

We use ICM (Iterative Conditional Mode) Optimization, optimize some while fixing others.

Here we fix segmentation parameters $\theta$, and optimize the transformation $\{G\},D$ .

$$
\begin{aligned}
\frac{\partial LL}{\partial G_{i,s}}
&=\sum_{x\in\Omega}\frac{1}{LH(x)}\sum_{k\in K}\prod_{j\ne i}P(I_j(x)\mid s(x)=k,\theta,G_{j,s})
\\
&\cdot\sum_c \tau_{ikc}\frac{\partial \Phi_{ikc}(I_i(G_{i,s}(x)))}{\partial (I_i(G_{i,s}(x)))}
\cdot
\frac{\partial (I_i(G_{i,s}(x)))}{\partial G_{i,s}(x)}
\cdot
\frac{\partial G_{i,s}(x)}{\partial  G_{i,s}}

\\
&=\sum_{x\in\Omega}\frac{1}{LH(x)}\sum_{k\in K}\prod_{j\ne i}P(I_j(x)\mid s(x)=k,\theta,G_{j,s})
\\
&\cdot\sum_c \tau_{ikc}\Phi_{ikc}'\nabla I_i(G_{i,s}(x)) \nabla G_{i,s}(x)

\end{aligned}
$$

For the deformation $D$

$$
\begin{aligned}
\frac{\partial LL}{\partial D}
& = \sum_{x\in\Omega} \frac{1}{LH(x)}\sum_{k\in K} \frac{\partial \pi_{kx\mid D}}{\partial D}P(I\mid s(x),\theta,\{G_{i,s}\})

\end{aligned}
$$


Where $\frac{\partial \pi_{kx\mid D}}{\partial D}=\frac{\partial A_{k}(D(x))}{\partial D}$,

$$
\frac{\partial A_{k}(D(x))}{\partial D}=\left.\nabla A_{k}\right|_{y=\left[y_{1}, y_{2}, y_{3}\right]} \times\left[\frac{\partial y_{1}}{\phi_{d}}, \frac{\partial y_{2}}{\phi_{d}}, \frac{\partial y_{3}}{\phi_{d}}\right]^{\mathrm{T}}
$$

Where $y=D(x),\{\phi_d\}$ are parameters of FFD.

So, the update is,

$$
\begin{aligned}
G_{i, s}^{[k+1]} & =G_{i, s}^{[k]}+l_{G} \cdot \frac{\partial L L}{\partial G_{i, s}} \\
D^{[k+1]} & =D^{[k]}+l_{D} \cdot\left(\frac{\partial L L}{\partial D}+\lambda \cdot \frac{\partial C_{\text {smooth }}}{\partial D}\right)
\end{aligned}
$$

where $\dfrac{\partial C_{\text {smooth }}}{\partial D}$ is a regularization term.

## Hetero-Coverage Multi-Modality Images (HCMMI)

Introduction of HCMMI: owing to each modality may have different coverage of the subject, so we emplopy multiple MvMMs to address this problem.

First, The whole common space is divided into several spaces:

$$
\Omega = \cup_{v=1}^{N_{sr}}\Omega_{v}
$$

where $\Omega_v$ expresses the $v_{th}$ non-overlapping region.

So our likelihood is the summation of each region $\Omega_v$,

$$
LL_{HC}=\sum_{v=1}^{N_{sr}}LL_{\Omega_v}=\sum_{v=1}^{N_{sr}}\sum_{x\in \Omega_v}\log\sum_{k\in K}\pi_{kx}\prod_{i\in N_v}(\sum_{c\in C_{ik}}\Phi_{ikc}(I_i(x)))
$$

The parameter update is simliar above, but $x\in\Omega_v$ instead of $x\in\Omega$ and $\prod_{i\in N_v}$ instead of $\prod_{i\in N_I}$. Specially, if $\|\Omega_v\|=1$, it's a GMM. If $\|\Omega_v\|\ge 2$, it's MvMM.
