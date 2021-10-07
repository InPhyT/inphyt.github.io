# Neuronal Model Calibration Framework

## Introduction

Hello,

We would like to introduce a neuronal model calibration framework, with the idea of being very user friendly for the beginners, but retaining high customization capabilities for the domain expert. 

This post is organized as follows: 

* [Section 1](#1. Overview) gives an overview of the issues addressed by this framework;
* [Section 2](#2. Worked Example) goes through the current interface to showcase some capabilities (we tried to keep it as light and readable as possible);
* [Section 3](#3. Possible Applications from Literature) provides some literature example that fits naturally within our framework;
* [Section 4](#4. Near Future Developments) concerns the developments we are going to make in the near future;
* [Section 5](#5. Long Term Developments / Needed Contributions) discusses the areas where we'd need contribution from people more experienced than us;
* [Section 6](#6. Example Calibration Results) illustrates some aspects of the current framework by looking at an user case with images.  

We have developed a model calibration pipeline, where the model is currently limited to be a [DifferentiaEquations.jl](https://github.com/SciML/DifferentialEquations.jl)-compatible system, and though our intentions are very ambitious (see [section 4](4. Near Future Developments) and [section 5](#5. Long Term Developments / Needed Contributions)), the full support is now only guaranteed for Hodgkin-Huxley (HH) type models. In this context, the pipeline focuses on calibrating HH-type models w.r.t. voltage traces coming from electrophysiological recordings.

The pipeline assists the modeler in the tasks of electrophysiological data collection, model specification and model fitting. Many calibration features have been implemented (the main focus has been around the methodology first devised in [Druckmann2007](https://www.frontiersin.org/articles/10.3389/neuro.01.1.1.001.2007/full), and then extended by numerous papers, e.g. Allen Brain Atlas's [technical whitepapers]()), and pretty much every main optimization technique that to our knowledge is available in Julia has been integrated (e.g. [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl), [Turing.jl](https://github.com/TuringLang/Turing.jl), [Evolutionary.jl](https://github.com/wildart/Evolutionary.jl), [Metaheuristics.jl](https://github.com/jmejia8/Metaheuristics.jl), [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl), etc. ) and they are all accessible and fully customizable through a unified interface.

## 1. Overview

Before having a look at the currently developed interface, let me be more specific about the purpose and the setting of the pipeline. We will also briefly present the ideas behind and cite the main references.

Compartmental neuron models are inherently hard to calibrate w.r.t. electrophysiological recordings, because of:

1. Their high degree of nonlinearity, which may induce:
   - **simulation issues**:  solver selection and optimization, computational time, etc. ;
   - **loss function complexity** ;
   - **train set identification**: difficulty in identifying a train set "that covers the whole dynamical range" of the neuron, which is equivalent to asking which is the best minimal set of stimuli to use as a train set, so that the calibrated model performs best on the test set, see [Druckmann2011](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002133) . 
2. **Stochasticity of the underlying physical system**: this implies that given a neuron, the same input (injected current) does not produce the same output (voltage trace) every time it is applied. This requires :
   - **feature-based parameter estimation**: since the voltage trace varies from sweep to sweep, its raw points cannot be used as a reliable way to measure a parameters set performance (and thus to calibrate a model). One must resort to the the extraction of features (defined by a domain expert) from that trace.
   - **calibration refinements**: it could be necessary to perform successive calibration techniques to obtain a better fit.

These issues have been addressed in our work as follows:

- **simulation issues**: we employed [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) efficient solving capabilities, and implemented an (almost) automated HH-type model conversion from the *gating variable rates formalism* (GVRF) to the *boltzmann terms formalism* (BTF), which makes out-of-dynamic-range derivatives more controllable during calibration ([Chen2010](https://ieeexplore.ieee.org/document/5491188) and [Thot2011](https://link.springer.com/article/10.1007%2Fs00422-011-0459-1)). This is achieved thanks to a set of structs representing membrane `Current`s and `Gating_Variable`s that are aware of the formalism they are written in ;
- **train set identification**: the pipeline allows for simultaneous multi-sweep calibration, where sweep grouping by stimulus type and amplitude is supported (and in some cases automatized), and so is the extraction of electrophysiological features means and standard deviations from such groups. These groups of sweeps are represented using a struct `Type_Amplitude`, which contains all sweeps obtained injecting the same stimulus (defined by a "type" which can be a *ramp*, a *square* and so on, and an "amplitude" expressed in pA). Thus the ideas behind [Druckmann2011](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002133) can be easily tested, so as every other experiment ;
- **feature-based parameter estimation**: the pipeline integrates electrophysiological features through a `Feature` struct, whose use is highlighted below. `Feature`s also support *combination*, as in [Hay2011](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002107), even though more general `combination`s have been made possible. `Feature`s values extracted from the data and from the simulated model are regularly compared during calibration, and the comparison is done through a loss function: this (together with many other details) can be independently specified for each `Feature`, and many common features have been already implemented. If the user wishes to adopt a custom feature, much support for [Allen SDK](https://allensdk.readthedocs.io/en/latest/) and [eFel](https://efel.readthedocs.io/en/latest/) is already there ;
- **complex loss functions**: this probably determined the poor performance of derivative-based optimization algorithm such as the BFGS from [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) (even though it's fair to say that we have little experience about hyperparameter tuning of optimization algorithms, see [section 5](#5. Long term developments / Needed Contributions)), in favor of evolutionary algorithms, especially the ones supporting multi-objective optimization, that could be found in [Metaheuristics.jl](https://github.com/jmejia8/Metaheuristics.jl) and [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl) ;
- **calibration refinements** : this pipeline has been built using a modeling package that we hope to publish in the near future. As such, it supports multi-step calibration, so that multiple calibration techniques may be stacked together, each one improving on the results of the previous one.

## 2. Worked Example

Being our work integrated with [Allen SDK](https://allensdk.readthedocs.io/en/latest/), let us now showcase how one may calibrate an HH-type model taken from [Pospischil2008](https://link.springer.com/article/10.1007%2Fs00422-008-0263-8) w.r.t. two `Type_Amplitude`s measured from cell number 488697163 in the [Allen Cell Types](https://celltypes.brain-map.org/) mouse database, one corresponding to a square pulse of amplitude 110pA, and the other to a ramp pulse of amplitude 800pA.  

First of all we're going to import all the relevant custom modules:

```julia
# main.jl

using ModelingFramework # This is a meta-package we developed that eases the production of domain-specific modeling packages like this one.

using Computational_Neuroscience # This module exports all computational neuroscience specific functionalities, such as the Sweep, Feature and Type_Amplitude structs and loss functions.

include(joinpath(@__DIR__,"src/Models/HH__Na_Kd_M_L_mf.jl" ))  # Load the model (see below).
include(joinpath(@__DIR__,"src/Models/Neuron_Parameters.jl" )) # A script containing all default parameter sets obtained from the literature for each HH-type model defined in Neuron_Models. Besides, it also contains all the membrane `Current`s derived from literature.
```

As stated in the last comment, the script `Neuron_Parameters.jl` contains all `Current`s. They are responsible for the *formalism switching* as described in the **simulation issues**. But how is a `Current` defined, and how does it allow for formalism switching ?<br>We will explore it in the next paragraphs.

### 2.1 Membrane currents definition

First we should remember that a  gating variable $x$ is written in the GVRF if it satisfies a differential equation of the form:

$$
\begin{equation} 
\frac{dx}{dt}(V) = \alpha(V)(1-x)-\beta(V)x
\end{equation} 
$$
In a specific module `Currents.jl`, we represented the $\alpha$s and the $\beta$s as `Gating_Variable_Rate`s:

```julia
# Currents.jl

# struct that represents a gating variable rate
mutable struct Gating_Variable_Rate{F<:Function}
    evaluation::F # the actual formula for the α or the β. 
	symbols::Vector{Symbol} # symbols representing the parameters contained in the `evaluation` field
end
```

A gating variable $x$ is instead written in the BTF if it satisfies:

$$
\begin{equation} 
\frac{dx}{dt}(V) = \frac{x_{\infin}(V) - x(V) }{\tau_x(V)}
\end{equation} 
$$
So we represented the $x_{\infin}$s and the $\tau_x$s as `Gating_Variable_Boltzmann_Term`s:

```julia
# Currents.jl

# struct that represents a gating variable boltzmann term
mutable struct Gating_Variable_Boltzmann_Term{F<:Function}
    evaluation::F # the actual formula for the x_∞ or the τₓ. 
	symbols # symbols representing the parameters contained in the `evaluation` field
end
```

Since simulation is more efficient in the BTF ([Thot2011](https://link.springer.com/article/10.1007%2Fs00422-011-0459-1)), we opted to represent gating variables always in that setting:

```julia
# Currents.jl

# struct that represents a gating variable
mutable struct Gating_Variable
    symbol # a symbol identifing this gating variable, e.g. :m
    x₀::Gating_Variable_Boltzmann_Term
    τ::Gating_Variable_Boltzmann_Term 
	# other fields that aren't relevant to our purposes
end
```

Anyway, the user is free to use a constructor that accepts `Gating_Variable_Rate`s, that internally calculates the `Gating_Variable_Boltzmann_Term`s and outputs the `Gating_Variable` in BTF:

```julia
# Currents.jl

# You may skip the implementation details, here reported just for completeness.
function Gating_Variable(symbol::Symbol; α::Gating_Variable_Rate, β::Gating_Variable_Rate) 
    symbols::Vector{Symbol} = vcat(α.symbols,  β.symbols)

    x₀_evaluation = function(V::Float64,p::Vector{Float64})
        return α.evaluation(V, p[1:length(α.symbols)])/(α.evaluation(V, p[1:length(α.symbols)]) + β.evaluation(V, p[length(α.symbols)+1:end]))
    end

    x₀::Gating_Variable_Boltzmann_Term = Gating_Variable_Boltzmann_Term(symbols,x₀_evaluation)

    τ_evaluation = function(V::Float64,p::Vector{Float64})
        return 1/( α.evaluation(V,p[1:length(α.symbols)]) + β.evaluation(V,p[length(α.symbols)+1:end]) )
    end

    τ::Gating_Variable_Boltzmann_Term = Gating_Variable_Boltzmann_Term(symbols,τ_evaluation)

    return Gating_Variable(symbol = symbol,  x₀ = x₀, τ = τ)
end
```

The link between the two formalisms is given by:

$$
\begin{equation}
x_{\infin}(V) = \frac{\alpha(V)}{\alpha(V)+\beta(V)}
\end{equation}
$$

$$
\begin{equation}
\tau_x(V) = \frac{1}{\alpha(V)+\beta(V)}
\end{equation}
$$

Of course, each membrane current $I$ is given by at most two gating variables (here $x$ and $y$) in the form:
$$
\begin{equation}
I(V) = g_I x^a y^b (V-V_I)
\end{equation}
$$
So a `Current` is represented as:

```julia
# Currents.jl

# struct representing a membrane current
mutable struct Current
    activating::Union{Gating_Variable,Nothing}   # the activating gating variable
    deactivating::Union{Gating_Variable,Nothing} # the inactivating gating variable
   
    activating_exponent::Float64
    deactivating_exponent::Float64
		# other fields that would just clutter the presentation, hidden from the user by the outer constructors.
		# ...
end
```

Thus, in the `Neuron_Parameters.jl`, you would find the sodium current `I_Na` (the *pospischil* interfix reminds of the article this current has been taken from, [Pospischil2008](https://link.springer.com/article/10.1007%2Fs00422-008-0263-8)) implemented as :

```julia
# Neuron_Parameters.jl

## I_Na
### m
#### αₘ
#### (V,p) -> ( m_α_A * ( V - V_T - 13 ) / ( exp( - ( V - V_T - 13 ) / m_α_C) - m_α_D ) )
const αₘ_pospischil_symbols = [:m_α_A ,:V_T,:m_α_C]
const αₘ_pospischil_evaluation(V,p) = ( - p[1] * ( V - p[2] - 13.0 ) / ( exp( - ( V - p[2] - 13.0 ) / p[3]) - 1.0 ) )
const αₘ_pospischil = Gating_Variable_Rate(αₘ_pospischil_symbols,αₘ_pospischil_evaluation )
#### βₘ
#### ( m_β_A * ( V_rest - V_T - 40 ) / ( exp( ( V_rest - V_T - 40 ) / m_β_C) - m_β_D ) )
const βₘ_pospischil_symbols = [:m_β_A,:V_T, :m_β_C ]
const βₘ_pospischil_evaluation(V,p) =  ( p[1] * ( V - p[2] - 40 ) / ( exp( ( V - p[2] - 40 ) / p[3]) - 1.0 ) )
const βₘ_pospischil = Gating_Variable_Rate(βₘ_pospischil_symbols, βₘ_pospischil_evaluation )

const m_pospischil = Gating_Variable(:m, α = αₘ_pospischil, β =  βₘ_pospischil)

### h
#### αₕ
#### (V,p) ->  ( h_α_A * exp( - ( V_rest - V_T - 17 ) / h_α_C ) )
const αₕ_pospischil_symbols = [:h_α_A ,:V_T , :h_α_C]
const αₕ_pospischil_evaluation(V,p) =  ( p[1] * exp( - ( V - p[2] - 17.0 ) / p[3] ) )
const αₕ_pospischil = Gating_Variable_Rate(αₕ_pospischil_symbols, αₕ_pospischil_evaluation  )

#### βₕ
#### ( m_β_A * ( V_rest - V_T - 40 ) / ( exp( ( V_rest - V_T - 40 ) / m_β_C) - m_β_D ) )
const βₕ_pospischil_symbols = [:h_β_A,:V_T ,:h_β_C ]
const βₕ_pospischil_evaluation(V,p) =  ( p[1] / ( 1.0 + exp( - ( V - p[2] - 40.0  ) / p[3]) ) )
const βₕ_pospischil = Gating_Variable_Rate(βₕ_pospischil_symbols, βₕ_pospischil_evaluation )

const h_pospischil = Gating_Variable(:h, α = αₕ_pospischil,  β =   βₕ_pospischil)

const I_Na_pospischil = Current(m_pospischil, h_pospischil)
```

This last snippet shows how the user may define his own currents. In the next paragraph we shall see how currents are used to switch from the GVRF to the BTF.

### 2.2 HH-type Model definition and formalism switching

Let us now turn to model definition. In the future, we plan to make it way easier (see section [HH-type model definition](#HH-type model definition)), but for now the user should just adopt the standard DifferentialEquations.jl interface. The name of the following model is to be interpreted as "HH-type model containing the `Na`, `Kd`, `M` and `L` current (taken from [Pospischil2008](https://link.springer.com/article/10.1007%2Fs00422-008-0263-8))".

```julia
# HH__Na_Kd_M_L.jl

## BT deterministic_dynamic
function HH__Na_Kd_M_L_BT_deterministic_dynamic!(I_series::F) where {F <:Function}
    function parametrized_HH__Na_Kd_M_L_BT_deterministic_dynamic!(du, u, p, t)

        C,   g_Na, V_Na, g_K,     g_M,   V_K,   g_L, V_Ca, g_l, V_l,
        V_m, κ_m,  τ₀_m, τ_max_m, V_τ_m, σ_m, 
        V_h, κ_h,  τ₀_h, τ_max_h, V_τ_h, σ_h,
        V_n, κ_n,  τ₀_n, τ_max_n, V_τ_n, σ_n,
        V_p, κ_p,  τ₀_p, τ_max_p, V_τ_p, σ_p,
        V_q, κ_q,  τ₀_q, τ_max_q, V_τ_q, σ_q,
        V_r, κ_r,  τ₀_r, τ_max_r, V_τ_r, σ_r        = p
        
        # State variables
        V = @view u[1]
        m = @view u[2]
        h = @view u[3] 
        n = @view u[4] 
        p = @view u[5]
        q = @view u[6]
        r = @view u[7]

        # State variables
        dV = @view du[1]
        dm = @view du[2]
        dh = @view du[3] 
        dn = @view du[4] 
        dp = @view du[5]
        dq = @view du[6]
        dr = @view du[7]
        
        dV .= (.-1.0 ./ C ) .* ( g_Na .* (m .^ 3) .* h .* ( V .- V_Na ) .+ g_K .* (n .^ 4 ) .* ( V .- V_K ) .+ g_M .* p .* (V .- V_K) .+ g_L .* (q.^2)*r*(V .- V_Ca ) .+ g_l .* ( V .- V_l ) .- I_series(t)   )
        @. dm = ( (1/2) * (1 + tanh( (V - V_m) / κ_m ) ) - m ) / ( τ₀_m + τ_max_m * ( 1 - tanh( ( V - V_τ_m  ) / σ_m ) ^ 2 ) )  
        @. dh = ( (1/2) * (1 + tanh( (V - V_h) / κ_h ) ) - h ) / ( τ₀_h + τ_max_h * ( 1 - tanh( ( V - V_τ_h  ) / σ_h ) ^ 2 ) )
        @. dn = ( (1/2) * (1 + tanh( (V - V_n) / κ_n ) ) - n ) / ( τ₀_n + τ_max_n * ( 1 - tanh( ( V - V_τ_n  ) / σ_n ) ^ 2 ) ) 
        @. dp = ( (1/2) * (1 + tanh( (V - V_p) / κ_p ) ) - p ) / ( τ₀_p + τ_max_p * ( 1 - tanh( ( V - V_τ_p  ) / σ_p ) ^ 2 ) ) 
        @. dq = ( (1/2) * (1 + tanh( (V - V_q) / κ_q ) ) - q ) / ( τ₀_q + τ_max_q * ( 1 - tanh( ( V - V_τ_q  ) / σ_q ) ^ 2 ) ) 
        @. dr = ( (1/2) * (1 + tanh( (V - V_r) / κ_r ) ) - r ) / ( τ₀_r + τ_max_r * ( 1 - tanh( ( V - V_τ_r  ) / σ_r ) ^ 2 ) ) 

    end
end
   
## noise
function HH__Na_Kd_M_L_noise!()
    
    function parameterized_HH__Na_Kd_M_L__noise!(du, u, p, t)

        du .= repeat([0.0], 7)

    end
    
end;

```

Notice that we didn't report here the GVRF dynamics: it is not needed.

Usually one finds the HH models and parameters values in the GVRF in the literature, but may want to switch to the BTF for better computational performance. This means that the reasearcher usually knows the GVRF parameters values (Cfr. the code above):

```julia
C,     g_Na,  V_Na,        g_K,   g_M,   V_K,   g_L,   V_Ca,  g_l, V_l, V_T,
m_α_A, m_α_C,              m_β_A, m_β_C,
h_α_A, h_α_C,              h_β_A, h_β_C,
n_α_A, n_α_C,       	   n_β_A, n_β_C,
p_0_B, p_0_C, τ_p_max_exp, τ_p_A, τ_p_B, τ_p_C,
q_α_A, q_α_B, q_α_C,       q_β_A, q_β_B, q_β_C,
r_α_A, r_α_B, r_α_C,       r_β_A, r_β_B, r_β_C
```

And wants to calculate the corresponding BTF parameters values so that each current's dynamical role is preserved (Cfr. the code above):

```julia
 C,   g_Na, V_Na, g_K,     g_M,   V_K,   g_L, V_Ca, g_l, V_l,
V_m, κ_m,  τ₀_m, τ_max_m, V_τ_m, σ_m, 
V_h, κ_h,  τ₀_h, τ_max_h, V_τ_h, σ_h,
V_n, κ_n,  τ₀_n, τ_max_n, V_τ_n, σ_n,
V_p, κ_p,  τ₀_p, τ_max_p, V_τ_p, σ_p,
V_q, κ_q,  τ₀_q, τ_max_q, V_τ_q, σ_q,
V_r, κ_r,  τ₀_r, τ_max_r, V_τ_r, σ_r
```

The researcher also wants to get the BTF equations from the GVRF equations, but this is straightforward since all BTF equations usually have the exact same form (more about it in the [BTF functional forms](#BTF functional forms) future development section), so one just needs to care about the number of them. In the following, we show how the pipeline eases this process.

From [Pospischil2008](https://link.springer.com/article/10.1007%2Fs00422-008-0263-8)  we extracted al the GVRF parameters values and coded them in a dictionary `GVRF_parameters_values` and placed it in the `Neuron_Parameters.jl` then we wrote the BTF model in the `parametrized_HH__Na_Kd_M_L_BTF_deterministic_dynamic` enclosed function (from a coding point of view it is a straightforward procedure, since BTF equations are all the same) . The last `HH__Na_Kd_M_L_noise!` closure is there for compatibility with Bayesian parameter estimation methods, but we won't go into the details of this here, so we will just leave it null.

A dedicated set of functionalities allows for automatically getting the parameters values for the BTF from the GVRF parameters values. Infact, in the `main.jl` we have:

```julia
# main.jl

BTF_parameters_values,BTF_currents = convert_GVRF_parameters_to_BTF_parameters([I_Na_pospischil, I_Kd_pospischil, I_M_pospischil, I_L_pospischil ], GVRF_parameters_values, plots = true  )
```

Where `I_Na_pospischil, I_Kd_pospischil, I_M_pospischil, I_L_pospischil` are all `Current`s written in the GVRF, while `GVRF_parameters_values` are the default GVRF parameters values derived from [Pospischil2008](https://link.springer.com/article/10.1007%2Fs00422-008-0263-8),

```julia
julia> GVRF_parameters_values

OrderedDict{Symbol,Float64}( :C           => 1.0   , :g_Na    => 50.0, :V_Na  => 50.0, :g_K => 5.0, :g_M => 0.004 , :V_K => -90.0, :g_L => 0.1, :V_Ca => 120.0,  :g_l => 0.01, :V_l => -70.61,  :V_T => -50.0 , 
                             :m_α_A       => 0.32    , :m_α_C => 4.0 ,                                  # αₘ
                             :m_β_A       => 0.28    , :m_β_C => 5.0 ,                                  # βₘ
                             :h_α_A       => 0.128   , :h_α_C => 18.0,                                  # αₕ
                             :h_β_A       => 4.0     , :h_β_C => 5.0 ,                                  # βₕ
                             :n_α_A       => 0.032   , :n_α_C => 5.0 ,                                  # αₙ
                             :n_β_A       => 0.5     , :n_β_C => 40.0,                                  # βₙ
                             :p_0_B       => 35.0    , :p_0_C => 10.0,                                  # p_x₀
                             :τ_p_max_exp => 4000.0  , :τ_p_A => 3.3 , :τ_p_B => 35.0, :τ_p_C => 20.0,  # p_τ
                             :q_α_A       => 0.055   , :q_α_B => 27.0, :q_α_C => 3.8 ,                  # α_q
                             :q_β_A       => 0.94    , :q_β_B => 75.0, :q_β_C => 17.0,                  # β_q
                             :r_α_A       => 0.000457, :r_α_B => 13.0, :r_α_C => 50.0,                  # α_r
                             :r_β_A       => 0.0065  , :r_β_B => 15.0, :r_β_C => 28.0                   # β_r 

)


```

The `plots` argument allows for a visual check that the fitting of the gating variables `m`,`n` etc according to the BTF functional forms was ok. It proved to be a good way of debugging copy-pasted BT and GVRF models. See [BTF functional forms](#BTF functional forms) to read how we may improve this process. The `BTF_parameters_values` is a dictionary containing the parameters values in the BTF:

```julia
julia> BTF_parameters_values

OrderedDict{Symbol,Float64}( :C     => 1.0   , :g_Na  => 50.0, :V_Na  => 50.0, :g_K => 5.0, :g_M => 0.004 , :V_K => -90.0, :g_L => 0.1, :V_Ca => 120.0,  :g_l => 0.01, :V_l => -70.61, 
                            :V_m => -24.289, :κ_m => 15.923,   :τ₀_m => 0.0287668,:τ_max_m => 0.0902789, :V_τ_m => -28.367, :σ_m => 42.521, 
                            :V_h => -28.37,  :κ_h => -7.98066, :τ₀_h => 0.029524, :τ_max_h => 5.93525,   :V_τ_h => -33.409, :σ_h => 14.874, 
                            :V_n => -25.782, :κ_n => 20.8588,  :τ₀_n => 0.169517, :τ_max_n => 1.51954,   :V_τ_n => -38.718, :σ_n => 34.855, 
                            :V_p => -35.0,   :κ_p => 20.0,     :τ₀_p => 14.836,   :τ_max_p => 1086.13,   :V_τ_p => -46.939, :σ_p => 30.256,
                            :V_q => -33.337, :κ_q => 8.87387,  :τ₀_q => 0.025755, :τ_max_q => 7.01596,   :V_τ_q => -38.257, :σ_q => 15.475,
                            :V_r => -58.488, :κ_r => -40.4226, :τ₀_r => 99.961,   :τ_max_r => 349.793,   :V_τ_r => -67.161, :σ_r => 65.92   

) 
```

While `BTF_currents` is an array of `Current`s written in the BTF. For calibration purposes, some algorithms need priors for those parameters. Since, to the best of our knowledge, such priors are often specified as ranges, we implemented a function that returns, for each parameter, a `Uniform` prior centered on its value and spanning a user-defined interval measured in orders of magnitude (`range`):

```julia
# main.jl

const parameters_priors_HH__Na_Kd_M_L_BTF = get_uniform_priors(BTF_parameters_values, range = 1 )
```

Initial conditions of (compartmental) models depend on parameters values (and resting potential, see [indico.ictp](http://indico.ictp.it/event/a13235/session/99/contribution/371/material/0/1.pdf)), so now that we have initial parameters values we may continue our HH-type model instantiation . This order is not compulsory and further arguments may be needed depending on the application, but we won't go into these details.

### 2.3 HH-type model instantiation

Our framework represents a model with a dedicated struct:

```julia
# HH_Model.jl

mutable struct HH_model{D <: Function, N <: Union{Nothing, Function}, E <: Function } <: AbstractDiffEqSystem
    deterministic_dynamic::D 
    noise::N
    calibration_history::Calibration_History    
    IC_evaluation_function::E
    model_kwargs::NamedTuple{<:Any, <:Tuple{Vararg{Any}}}
end

```

Where:

- `deterministic_dynamic`: in our case, this would correspond to the deterministic dynamic closure `HH__Na_Kd_M_L_deterministic_dynamic!` ;
- `noise`: the pipeline supports Bayesian parameter optimization (here not covered), featuring integration with [Turing.jl](https://github.com/TuringLang/Turing.jl). In this report it will be set to a null system ;
- `calibration_history`: this field keeps track of the calibrated parameters, initial conditions and (multivariate) posteriors during calibration, here not covered. This can be initialized with initial parameters values and initial conditions ;
- `model_kwargs`: allows for passing arguments to `deterministic_dynamic`and `noise` during calibration ;
- `IC_evaluation_function`: the field representing the function used to calculate initial conditions from resting potential and parameters values. The user should implement this function at least for the BTF, following the boilerplate:

```julia
# HH__Na_Kd_M_L.jl

## see indico.ictp : http://indico.ictp.it/event/a13235/session/99/contribution/371/material/0/1.pdf
function get_HH__Na_Kd_M_L_BTF_IC(V_rest::Float64, all_parameters_values::Vector{Float64})
    C,   g_Na, V_Na, g_K,     g_M,   V_K,   g_L, V_Ca, g_l, V_l,
    V_m, κ_m,  τ₀_m, τ_max_m, V_τ_m, σ_m, 
    V_h, κ_h,  τ₀_h, τ_max_h, V_τ_h, σ_h,
    V_n, κ_n,  τ₀_n, τ_max_n, V_τ_n, σ_n,
    V_p, κ_p,  τ₀_p, τ_max_p, V_τ_p, σ_p,
    V_q, κ_q,  τ₀_q, τ_max_q, V_τ_q, σ_q,
    V_r, κ_r,  τ₀_r, τ_max_r, V_τ_r, σ_r     = all_parameters_values

    m₀::Float64 = (1/2) * (1 + tanh( (V_rest - V_m) / κ_m ) )
    h₀::Float64 = (1/2) * (1 + tanh( (V_rest - V_h) / κ_h ) )
    n₀::Float64 = (1/2) * (1 + tanh( (V_rest - V_n) / κ_n ) )
    p₀::Float64 = (1/2) * (1 + tanh( (V_rest - V_p) / κ_p ) )
    q₀::Float64 = (1/2) * (1 + tanh( (V_rest - V_q) / κ_q ) )
    r₀::Float64 = (1/2) * (1 + tanh( (V_rest - V_r) / κ_r ) )

    return vcat(V_rest, m₀, h₀, n₀, p₀, q₀, r₀ )
end
```

Thus one may instantiate the model via a dedicated constructor:

```julia
# main.jl

HH__Na_Kd_M_L_BTF= Compartmental_Model(HH__Na_Kd_M_L_BTF_deterministic_dynamic!, parameters_symbols_HH__Na_Kd_M_L_BTF, parameters_priors_HH__Na_Kd_M_L_BTF, get_HH__Na_Kd_M_L_BTF_IC; initial_parameters_values = BTF_parameters_values);
```

Let us now turn to experimental dataset loading.

### Dataset loading

As already stated, the pipeline is interfaced to [Allen SDK](https://allensdk.readthedocs.io/en/latest/) via [PyCall.jl](https://github.com/JuliaPy/PyCall.jl). This means that we may load all information stored in the Allen Brain Atlas regarding cell 488697163 of the mouse cortex this way:

```julia
# main.jl

const manifest_path = raw"path\to\manifest.json" # path to previously-created allensdk manifest, or to the location where we want to save the manifest at.
const cell_dataset  = get_cell_dataset(manifest_path,488697163) 
```

See [Dataset loading developments](#Dataset loading developments) for future improvements of this section.

The `cell_dataset` variable is an `OrderedDict` containing:

```julia
julia> cell_dataset
Dict{String, Any} with 3 entries:
"sweeps" => DefaultDict{String, DefaultOrderedDict{String, Vector{Sweep}, F} where F, Vector{Sweep}}("ramp"=>DefaultOrderedDict("800.0 pA"=>[Sweep{Float64}([-94.2813, -94.25, -94.2813, -94.3125, -94.3125, -94.1875, -94.2813, -94.2813, -94.2813, -94.3125  … …  
"cell_ephys_features"  => DataFrameRow…
"cell_sweeps_features" => Dict{Any, Any}()
```

For this brief overview, we are interested in the `sweeps` key, which points to a `DefaultDict` that groups the sweeps by input type:

```julia
julia> cell_dataset["sweeps"]
DefaultDict{String, DefaultOrderedDict{String, Vector{Sweep}, F} where F, Vector{Sweep}} with 8 entries:
"ramp" => DefaultOrderedDict("800.0 pA"=>[Sweep{Float64}([-94.2813, -94.25, -94.2813, -94.3125, -94.3125, -94.1875, -94.2813, -94.2813, -94.2813, -94.3125  …  -14.0, -14.0, -14.0, -14.0, -14.0, -14.0, -14.0, -14.0, -14.0, -14.0], [0.0, 0.0, 0.0, 0.…  
"square - 2s suprathreshold" => DefaultOrderedDict("110.0 pA"=>[Sweep{Float64}([-94.75, -94.75, -94.6875, -94.7188, -94.75, -94.7188, -94.75, -94.75, -94.6875, -94.7188  …  -94.7813, -94.75, -94.75, -94.6875, -94.7188, -94.7813, -94.75, -94.75, -94.7188, -94.7813], [0.0…  
"test" => DefaultOrderedDict("0.0 pA"=>[Sweep{Float64}([-13.1875, -13.1875, -13.125, -13.1875, -13.1875, -13.1875, -13.1562, -13.1562, -13.1875, -13.1875  …  -13.2187, -13.2812, -13.1875, -13.125, -13.1875, -13.2187, -13.2187, -13.2187, -13.1562, -…  
"long square" => DefaultOrderedDict("-110.0 pA"=>[Sweep{Float64}([-93.8438, -93.8125, -93.7813, -93.8125, -93.8125, -93.875, -93.8125, -93.8438, -93.7813, -93.8125  …  -95.125, -95.1563, -95.1563, -95.1563, -95.1875, -95.125, -95.1563, -95.125, -95.0938, …  
"square - 0.5ms subthreshold" => DefaultOrderedDict("-200.0 pA"=>[Sweep{Float64}([-93.875, -93.875, -93.875, -93.8125, -93.8438, -93.7813, -93.8438, -93.8125, -93.875, -93.8438  …  -94.0625, -94.0313, -94.0313, -94.0625, -94.0625, -94.0313, -94.0313, -94.0625, -94.0938, …  
"noise 2" => DefaultOrderedDict("259.375 pA"=>[Sweep{Float64}([-95.0, -95.0, -95.0, -95.0, -95.0, -94.9375, -95.0, -94.9688, -94.9688, -95.0  …  -93.6875, -93.7188, -93.7813, -93.7188, -93.7813, -93.7188, -93.7188, -93.7813, -93.6875, -93.75], [0.0, 0…  
"short square" => DefaultOrderedDict("400.0 pA"=>[Sweep{Float64}([-94.0313, -94.0313, -94.0625, -94.0, -94.0313, -94.0313, -94.0, -94.0313, -94.0313, -94.0625  …  -93.9063, -93.875, -93.9063, -93.9688, -93.875, -93.9063, -93.9063, -94.0313, -93.875, -93.90…  
"noise 1" => DefaultOrderedDict("255.0 pA"=>[Sweep{Float64}([-95.1563, -95.1563, -95.125, -95.1563, -95.1563, -95.125, -95.125, -95.125, -95.1563, -95.125  …  -95.875, -95.8125, -95.8438, -95.8438, -95.8438, -95.8438, -95.8438, -95.8125, -95.8438, -95…
```

And then by input amplitude:

```julia
cell_dataset["sweeps"]["square - 2s suprathreshold"]
DefaultOrderedDict{String, Vector{Sweep}, DataType} with 3 entries:
"110.0 pA" => [Sweep{Float64}([-94.75, -94.75, -94.6875, -94.7188, -94.75, -94.7188, -94.75, -94.75, -94.6875, -94.7188  …  -94.7813, -94.75, -94.75, -94.6875, -94.7188, -94.7813, -94.75, -94.75, -94.7188, -94.7813], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  …  
"150.0 pA" => [Sweep{Float64}([-94.3438, -94.375, -94.375, -94.3125, -94.3125, -94.3438, -94.3438, -94.3438, -94.375, -94.4063  …  -94.0, -94.0313, -93.9375, -94.0, -93.9375, -93.9688, -93.9375, -94.0, -93.9688, -93.9375], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.…  
"190.0 pA" => [Sweep{Float64}([-94.9688, -94.9688, -95.0313, -94.9688, -94.9375, -94.9688, -94.9688, -94.9688, -95.0, -94.9688  …  -94.7813, -94.7813, -94.8438, -94.8125, -94.7813, -94.9063, -94.7813, -94.8438, -94.8125, -94.8125], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0…
```

Once one has selected the type and the amplitude, all the recorded experimental sweeps corresponding to that stimulus configuration are presented:

```julia
julia> cell_dataset["sweeps"]["square - 2s suprathreshold"]["110.0 pA"]
4-element Vector{Sweep}:
 Sweep{Float64}([-94.7500034570694, -94.7500034570694, -94.6875005364418, -94.7187557220459, -94.7500034570694, -94.7187557220459, -94.7500034570694, -94.7500034570694, -94.6875005364418, -94.7187557220459  …  -94.7812511920929, -94.7500034570694, -94.7500034570694, -94.6875005364418, -94.7187557220459, -94.7812511920929, -94.7500034570694, -94.7500034570694, -94.7187557220459, -94.7812511920929], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.034999999999999996, 0.04, 0.045000000000000005, 0.05  …  10019.96, 10019.964999999998, 10019.970000000001, 10019.975, 10019.98, 10019.985, 10019.99, 10019.994999999999, 10020.0, 10020.005], (150000, 2004000), (v = 0.001, i = 1.0e-6, t = 0.001), (sweep_number = 61, sampling_rate_Hz = 200000.0, number_of_spikes = 12, metadata = Dict{Any, Any}("num_spikes" => 12, "stimulus_relative_amplitude" => 1.10000002384186, "stimulus_name" => "Square - 2s Suprathreshold", "pre_vm_mv" => -80.482666015625, "bridge_balance_mohm" => 23.0769233703613, "slow_vm_mv" => -80.482666015625, "post_vm_mv" => -80.8379364013672, "id" => 488697462, "stimulus_units" => "Amps", "pre_noise_rms_mv" => 0.0377454534173012…)))
 Sweep{Float64}([-95.18750154972076, -95.15625381469727, -95.15625381469727, -95.12500607967377, -95.15625381469727, -95.15625381469727, -95.18750154972076, -95.28125220537186, -95.18750154972076, -95.18750154972076  …  -94.50000667572021, -94.50000667572021, -94.46875149011612, -94.43750375509262, -94.43750375509262, -94.43750375509262, -94.53125441074371, -94.46875149011612, -94.50000667572021, -94.43750375509262], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.034999999999999996, 0.04, 0.045000000000000005, 0.05  …  10019.96, 10019.964999999998, 10019.970000000001, 10019.975, 10019.98, 10019.985, 10019.99, 10019.994999999999, 10020.0, 10020.005], (150000, 2004000), (v = 0.001, i = 1.0e-6, t = 0.001), (sweep_number = 62, sampling_rate_Hz = 200000.0, number_of_spikes = 13, metadata = Dict{Any, Any}("num_spikes" => 13, "stimulus_relative_amplitude" => 1.10000002384186, "stimulus_name" => "Square - 2s Suprathreshold", "pre_vm_mv" => -80.9035186767578, "bridge_balance_mohm" => 23.0769233703613, "slow_vm_mv" => -80.9035186767578, "post_vm_mv" => -80.4853668212891, "id" => 488697464, "stimulus_units" => "Amps", "pre_noise_rms_mv" => 0.0276526212692261…))
...
```

Each `Sweep` is a struct representing a sweep. It has the following definition:

```julia
# _Sweep.jl

struct Sweep{T}
    v::Vector{T}
    i::Vector{Float64}
    t::Vector{Float64}
    index_range::Tuple{Int64,Int64}
    UOM::NamedTuple{(:v, :i, :t), Tuple{Float64, Float64, Float64}}
    extra_data::NamedTuple{<:Any, <:Tuple{Vararg{Any}}}
end
```

Where:

- `v`: the membrane potential time series ;
- `i`: the stimulus time series ;
- `t`: the time steps at which `v` and `i` are simultaneously recorded. Equal time stepping [is required for Allen SDK to work properly](https://github.com/AllenInstitute/AllenSDK/issues/2210) ;
- `index_range`: index range (subset of `1:length(t)`) that includes the stimulus but excludes test pulses. See [Allen Elecrophysiology Overview Whitepaper](http://help.brain-map.org/display/celltypes/Documentation) ;
- `UOM`: units of measure with respect to SI of the three time series `v`,`i` and `t`. So if their units are , respectively, mV, pA and sec, `UOM` will be `(v = 1e-3, i = 1e-12, t = 1e0)`. This field is required so that the modeler may use the units she/he prefers (this may influence solver's performance), while the pipeline is constantly aware of current units and may convert them when needed ;
-  `extra_data`: extra information about the sweep or the cell it comes from.

Now that we have the model and the dataset, we only need one more ingredient, `Feature`s, and then we'll be ready for calibration setup.

### 2.4 Features

As stated in the **feature-based parameter estimation** bullet point, it makes little sense to use the raw membrane potential time series as reference data to calibrate against. A more optimal solution would be the identification of a set of features of these time series that is both minimalistic and fully encodes electrophysiological dynamics. Our framework supports the use of features via the `Feature` struct:

```julia
# _Feature.jl

mutable struct Feature{E <: Function,L <: Function,M <: Function,S <: Union{Function,Nothing}, G <: Function}
    name::String
    default_standard_deviation::Union{Float64,Nothing}
    UOM::Float64
    evaluation::E
    loss::L
    mean_function::M
    std_function::S
    kwargs::NamedTuple{<:Any, <:Tuple{Vararg{Any}}}
    loss_kwargs_evaluation::G
end
```

Some relevant fields:

- `evaluation`: the actual "formula" that intuitively takes a sweep an input and evaluated the corresponding feature ;
-  `loss`: the loss function associated to this feature, i.e. a function that takes two values of this `Feature` (one coming from the experimental trace, and the other from the set of parameters that is currently evaluated during calibation), and returns a `Float64`; 
- other fields: necessary fields to allow for the most general features to be implemented .

Now that we have features, we may move on to setup the calibration.

### 2.5 Calibration setup

The idea behind the following setup is the utilization of a framework that can be generalized to as many other fields than neuroscience as possible. The information needed to setup the calibration is a set of associations between experimental data and features that we want to calibrate our model against simultaneously.

In our case, the "unit of data" is a group of sweeps coming from the same cell when injected with a particular stimulus (defined by its type and amplitude), that we encoded in a `Type_Amplitude` struct:

```julia
# _Sweep.jl

## a struct representing a type_amplitude
struct Type_Amplitude <: AbstractMeasurement
    type::String
    amplitude::String
    sweeps::Tuple{Vararg{Sweep}}
end
```

We may get some `Type_Amplitude` instances by using an utility function on the `cell_dataset`:

```julia
# main.jl

const square_110pA = get_type_amplitude(cell_dataset, "square - 2s suprathreshold", "110.0 pA")
const ramp_800pA   = get_type_amplitude(cell_dataset, "ramp", "800.0 pA")
```

Also, the pipeline exports some standard electrophysiological features that we group as suggested in [Druckmann2011](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002133):

```julia
# main.jl

const square_features =  (allen_firing_frequency_Hz, allen_latency_to_first_spike_msec, allen_isi_cv, allen_spike_width_at_hh_msec,   allen_action_potential_peak_mV, allen_fast_trough_depth_mV)
const ramp_features =  (allen_firing_frequency_Hz, allen_latency_to_first_spike_msec, allen_isi_cv, allen_spike_width_at_hh,   allen_action_potential_peak_mV, allen_fast_trough_depth_mV, spike_height_slope)
```

Here I just used 7, but there are about 24 features implemented out of the box, coming from the various [Allen technical whitepapers]() and the cited relevant literature, that could be roughly divided in two groups: Druckmann-like (features compatible with a loss that takes experimental standard deviations into account), and L2-like (features that are best suited with a L2 loss, like [LeMasson and Maex, 2001](https://www.researchgate.net/publication/329683440_Introduction_to_equation_solving_and_parameter_fitting). and Maex's phase plane trajectory density). Some features that do not fall into these paradigms, like [Xu2019](https://ieeexplore.ieee.org/document/8682825)'s coincidence factor, can anyway be integrated with the framework.

Now we want to set up associations, to tell the pipeline which features we would like to calibrate w.r.t. which data, here we choose to calibrate both `Type_Amplitude`s w.r.t. the same features set, so we write:

```julia
# main.jl

const square_110_association = Association("square_110_association",square_110pA,square_features)
const ramp_800_association = Association("ramp_800_association",ramp_800pA,ramp_features)
```

To understand the remainder, we need a brief excursus on multi-phase calibration. It is a calibration technique that is useful when:

- there are external forces driving the physical system dynamics that are not included in the model ;
- these external forces only apply at discrete times, the number of which and/or their location in time being known in advance.

Say that we are in the case where we know when those unaccounted external forces happen. The time intervals between each external intervention will be called `Phase`, and within each `Phase` the modeler could be interested in developing a calibration strategy composed of different passages, each of which will be represented by a `Step` struct.

When it comes to neuronal models, multi-phase calibration is not needed (we just need one `Phase`, since we calibrate w.r.t. to the entire train set at once), but one may be interested in using different algorithms and setups in a particular order, or better, one may need `Step`s:

```julia
# ModelingFramework.jl

mutable struct Step{F <: Union{Nothing, <: Function} } <: AbstractStep
    name::String       # The name of the step. Used for later storage in `caibration_history`
    method::F          # It must have signature (data::D, Model::M; kwargs...) where {D <: Tuple{Vararg{<:AbstractData}}, M<: AbstractModel}
    kwargs::NamedTuple # kwargs to be passed to `method`
en`

```

Thus, we may instantiate a Step as follows:

```julia
# main.jl

const save_idxs = [1] # corresponds to voltage trace
bbo_borgmoea_step = Step("bbo_borgmoea_step",bbo_borgmoea_method, (variable_parameters_symbols = parameters_symbols_HH__Na_Kd_M_L_BTF, save_idxs = save_idxs, modelingtoolkit = true, cell_dataset_optimization = Cell_Dataset_Optimization( true,"block_average",  excluded_features = ("custom_voltage_base_mV_druck","custom_voltage_base_mV_l2","custom_steady_state_voltage_stimend_mV_druck", "custom_steady_state_voltage_stimend_mV_l2", "custom_voltage_after_stim_mV_druck", "custom_voltage_after_stim_mV_l2", "spike_times_msec")), maxtime = 30*60 ) )
```

The `bbo_borgmoea_method` is a method exported by the `ComputationalNeuroscience` that implements calibration using the `:borg_moea` algorithm from [BlackBoxOptim.jl](https://github.com/robertfeldt/BlackBoxOptim.jl). There are analogous methods tha tgive access to the algorithms from [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl), [Turing.jl](https://github.com/TuringLang/Turing.jl), [Evolutionary.jl](https://github.com/wildart/Evolutionary.jl) and [Metaheuristics.jl](https://github.com/jmejia8/Metaheuristics.jl). All the arguments that follow are given as `kwargs` to `bbo_borgmoea_method`. The `parameters_symbols_HH__Na_Kd_M_L_BTF` argument allows for quickly including/excluding individual parameters from calibration (here we choose to calibrate all of them).   The `modelingtoolkit` argument tells the pipeline to `modelingtoolkitize` the model, for incresed performance. To allow for it, stimuli are internally represented as functions (right now only square and ramp stimuli support this feature, we plan to add others in the future). The `Cell_Dataset_Optimization(true, "block_average")` specifies that we wish to apply an heuristic to trim the sweeps so that only the relevant portion and some equal margin is considered, while the `"block_average"` argument specifies the technique to remove points in the sweep in order to determine lower bound on solver stepping and perform other secondary optimization. The optimization procedure is automatically stopped when it detects that a feature of the optimized dataset differs from the same feature of the unoptimized one more than an experimental standard deviation.

We will perform a single-step calibration.

```julia
single_step_calibration!((square_110_association, ramp_800_association), HH__Na_Kd_M_L_BTF, bbo_borgmoea_step)
```

At the end of the calibration, resulting parameters sets (together with optimized initial conditions and posterior if required by custom calibration techniques), may be found in the `calibration_history` field of the model:

```julia
julia> get_calibratables_values(get_last_configuration(HH__Na_Kd_M_L_BTF.calibration_history).parameters)
46-element Vector{Float64}:
   8.329481760571682
  70.5016621571624
  62.09685457172623
   9.858797725026903
   0.009504285040026086
 -80.1864762458009
   ⋮
 -64.14007514627994
 -80.64025855432278
  15.919456056391812
 164.9340181305733
 -25.45128186775263
  80.67340483048163
```



## 3. Possible Applications from Literature

Although this work is far from over, we believe there is some value in it, the interested reader may confirm that applications such as  [Buhry2012](https://www.sciencedirect.com/science/article/abs/pii/S0925231211006618?via%3Dihub) easily fit in our framework, since the necessary algorithms are already integrated, and the possibility of defining custom `Feature`s  is already there. On the other hand, the integration of techniques like the one found in [Thot2011](https://link.springer.com/article/10.1007%2Fs00422-011-0459-1), may still require some effort. Other techniques, like the one based on the *coincidence factor* in [Xu2019](https://ieeexplore.ieee.org/document/8682825) and the *phase plane trajectory density* in [Achard2006](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0020094) and [LeMasson and Maex, 2001](https://www.researchgate.net/publication/329683440_Introduction_to_equation_solving_and_parameter_fitting), have already been integrated.

## 4. Near Future Developments

### Recipe

We aim to give a quick working recipe to effectively calibrate many HH-type (maybe even compartmental) models.

### HH-type model definition

We will implement a system that allows the user to define an HH model as follows:

```julia
HH_Na_Kd_L = HH_model(I_Na, I_Kd, I_L)
```

## 5. Long Term Developments / Needed Contributions

Developing a model calibration framework requires a large and heterogeneous amount of expertise that we are currently missing. 

Although we haven't given up yet, despite following the advices in [Druckmann2011](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002133) the optimization algorithms always output a "best parameters set" that is indifferent to the applied stimulus (within a certain range). See [section 6](#6. Example Calibration Results) for a better visualization of this issue. We probably ruled out the possibility that it is a bug, and it may be correlated to the following other areas where we lack expertise:
- **Optimization algorithm tuning**. Here is a [Twitter thread]() with references and the question in more details when it comes to evolutionary algorithms. Also, [Turing.jl](https://github.com/TuringLang/Turing.jl) is fully supported and usable with the same interface above (just set `method = :turing_neuro`), and although it works well on another project we are developing on epidemiological models, it achieves very poor performance on HH-type models.  
- **Modeling choices**. It may be that the HH-type models we tried so far (the most complete of them being the one shown above) may be inappropriate, or we could be missing some relations between the cell we select the train sweeps from and the model itself. We don't know whether multi-compartmental models are unavoidable, as the literature seems to suggest, or if we could dramatically improve calibration performance by just more fine tuning of the single-compartment models we are currently using. Also, this framework may need modifications in orderer to be compatible with GLIFs and multi-compartmental models.
- **Optimization techniques**. We have not been able to implement every possible calibration technique. A particularly relevant and missing one are Kalman Filters. Some possible references on the subject may be [Lankarany2014](https://www.sciencedirect.com/science/article/abs/pii/S0925231214001155?via%3Dihub), [Campbell2019](https://www.mdpi.com/2076-3417/10/2/550) (for the inverse problem), [Moye2018](https://mathematical-neuroscience.springeropen.com/articles/10.1186/s13408-018-0066-8) and [Lillacci2010](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000696), but there may be better ones. Any help would be greatly appreciated. Other techniques we did not implement are the ones presented in [Vavoulis2012](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002401#pcbi.1002401-Kitagawa1) , [Haufler2007](https://www.sciencedirect.com/science/article/abs/pii/S092523120600347X) and [Meliza2014](https://link.springer.com/article/10.1007%2Fs00422-014-0615-5).

We would be very grateful if anybody could help with these issues. 

## 6. Example Calibration Results

Here is the output of an 8-hour run, where we calibrated the above `HH__Na_Kd_M_L` model w.r.t. the `square_110pA` sweeps set (which counted 4 sweeps), using the following features:

- **time_of_slow_trough_fraction**: time from the spike peak to the slow through, divided by the `isi` length, adimensional ;
- **AP_width_hh_msec**:  width at half-height of the action potential, expressed in msec ;
-  **adaptation_index**: adaptation index, adimensional ;
-  **firing_frequency_Hz**: firing frequency, expressed in Hz ;
-  **isi_cv**: interspike interval (durations) coefficient of variation, adimensional;
-  **latency_to_first_spike_msec**: time from start of recording to first peak, expressed in msec ;
-  **duration_of_first_isi_msec**: duration of first interspike interval, expressed in msec ;
-  **action_potential_peak_mV**: action potentials' peaks height, expressed in mV ;
-  **fast_trough_depth_mV**: fast through depth, expressed in mV ;
-  **average_isi_msec**: average interspike interval duration, expressed in msec ;
-  **resting_potential_mV**: resting potential, expressed in mV ;
-  **slow_trough_depth_mV**: slow through depth, expressed in mV .

More details about these features are provided in the Allen technical whitepaper [Electrophysiology Overview](http://help.brain-map.org/display/celltypes/Documentation).

The achieved performance, in units of experimental standard deviation (see [Druckmann2007](https://www.frontiersin.org/articles/10.3389/neuro.01.1.1.001.2007/full)) is:

- **time_of_slow_trough_fraction**: 2.23066 ;
- **AP_width_hh_msec**:  5.36808 ;
-  **adaptation_index**: 1.02698 ;
-  **firing_frequency_Hz**: 1.96846 ;
-  **isi_cv**: 1.32066 ;
-  **latency_to_first_spike_msec**: 3.13225 ;
-  **duration_of_first_isi_msec**: 0.12058 ;
-  **action_potential_peak_mV**: 7.20770 ;
-  **fast_trough_depth_mV**: 1.49870 ;
-  **average_isi_msec**: 1.76509 ;
-  **resting_potential_mV**: 9.74667 ;
-  **slow_trough_depth_mV**: 3.67717 

which amounts to an aggregated loss of 39.063, that while being not too far from the 3 standard deviations per features goal (`3*12 = 36`), it is not very much satisfying either. Infact, if one looks at a simulated trace (blue) plotted together with one of the sweeps used for calibration:

![calibration_results.PNG](https://github.com/InPhyT/Twitter/blob/main/Neuroscience/reports/images/calibration_results.PNG?raw=true)

She/he may notice that:

- Spike coincidence is not a visual feature to take into consideration when qualitatively judging the results of a calibration run, since there is great variance in spike times (see plot below) ;
- The parameters set found does not take the injected current into consideration, as it keeps spiking well beyond the end of the current stimulus.

To address the second issue, we tried the following strategies without success:

- incorporate a `latency_to_last_spike` feature ;
- calibrate w.r.t. 4 `Type_Amplitude`s (recall that each  `Type_Amplitude` is a sweeps set corresponding to one stimulus waveform ), 3 of which did not produce any spike in the real neuron ;

None of these secceded.

We wrap up by showing more plots relative to the above run:

**Model currents, gating variables dynamics (with zoom in on first spike)**:

![time_series.PNG](https://github.com/InPhyT/Twitter/blob/main/Neuroscience/reports/images/time_series.PNG?raw=true)

**Same with voltage on the x-axis**:

![phase_space.PNG](https://github.com/InPhyT/Twitter/blob/main/Neuroscience/reports/images/phase_space.PNG?raw=true)

**All the experimental sweeps used**:

![type_amplitude.PNG](https://github.com/InPhyT/Twitter/blob/main/Neuroscience/reports/images/type_amplitude.PNG?raw=true)

In our opinion, this example shows both the possibilities of this framework, and the issues that still persist.

We hope to have even slightly caught your attention, don't hesitate to reach out if you think you may give some much needed help and suggestions!

[@johnmyleswhite](https://twitter.com/johnmyleswhite), [@KristofferC89](https://twitter.com/KristofferC89), [@briochemc](https://twitter.com/briochemc), [@benskuhn](https://twitter.com/benskuhn), [@drfeldt](https://twitter.com/drfeldt), [@ranjan_ananth](https://twitter.com/ranjan_ananth), [@baggepinnen](https://twitter.com/baggepinnen), [@ChrisRackauckas](https://twitter.com/ChrisRackauckas) [@vchuravy](https://twitter.com/vchuravy), [@rogerluorl18](https://twitter.com/rogerluorl18)  [@TuringLang](https://twitter.com/TuringLang), [@xukai92](https://twitter.com/xukai92), [@Hong_Ge2](https://twitter.com/Hong_Ge2), [@cameron_pfiffer](https://twitter.com/cameron_pfiffer), [@martin_trapp](https://twitter.com/martin_trapp), [@MathieuEmile](https://twitter.com/MathieuEmile), [@torfjelde](https://twitter.com/torfjelde), [@sethaxen](https://twitter.com/sethaxen), [@ANoackJensen](https://twitter.com/ANoackJensen), [@HarrisonDWilde](https://twitter.com/HarrisonDWilde)

@AlexGuetMcC , @AmitChaudhariMD, @Anna__K__K , @bradleyvoytek   , @misicbata       , @ChristianPehle  , @DamienQuerlioz  , @cdmeliza   , @MelissaRWarden  , @DrValiante      , @gmiotto         , @Segev_Lab   , @neurodidact     , @julie_grollier  , @lauracrucianel1 , @ogrisel         , @TheBrunoCortex  , @wake_sleep  , @RichCompNeuro   , @sergeydoestweet , @ShaulDr     , @neuronJoy       , @srikipedia      , @SylvainSaighi

@bcm_neurosci   ,@ChenLabNeuro   ,@ChicagoSfN     ,@C_R_A_N_I_A    ,@EITN_Paris     ,@ELSCbrain      ,@EPFL_en        ,@BlueBrainPjt   ,@INCForg        ,@HBPBrainSim    ,@HBP_Education  ,@KCNHub         ,@NeuralCompLab  ,@NeurIPSConf    ,@NIN_knaw       ,@NeuroPSI_saclay,@RadioSpin_EU   ,@SurLabMIT      















