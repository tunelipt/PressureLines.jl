# Correção da pressão nos tubos
import SpecialFunctions: besselj, besselj0
import FFTW: rfft, irfft, rfftfreq
export TubeSection, tubecoefs, pressratio, phaseangle, amplphase
export PressureLine, PressCorrect

#using Makie
#export plot_pressline!, plot_pressline


struct TubeSection{T<:AbstractFloat}
    "Length of tube"
    L::T
    "Radius of the tube"
    R::T
    "Volume at the end of the tube section"
    V::T
    "Effective volume ratio"
    Vr::T
    "Isentropic coefficient for the fluid"
    γ::T
    "Dynamic viscosity of the fluid"
    μ::T
    "Density of the fluid"
    ρ::T
    "Velocity of sound of the fluid"
    a₀::T
    "Prandtl number of the fluid"
    Pr::T
    "Fractional volumetric increase due to diaphragm deflection"
    σ::T
    "Polytropic constant for the volume"
    k::T
end
Base.broadcastable(tb::TubeSection) = Ref(tb)

Base.show(io::IO, tb::TubeSection) = print(io, "Tube section: L=$(tb.L*1000) mm, D=$(tb.R*2000) mm, V=$(tb.V*1e9) mm³")


"""
`TubeSection(;D=1.0e-3, L=0.2e-3, V=1e-9, sigma=0.0,  gamma=1.4, k=gamma,
                     mu=1.5e-5, rho=1.2, a0=343.0, Pr=0.707)`

Create a new `TubeSection` object with given parameters. A tube section is a part of
a tubing system or pressure line. It consists of a tube with internal diameter `D`,
length `L` connected to a volume `V`. 
"""
function TubeSection(;D=1.0e-3, L=0.2e-3, V=1e-9, sigma=0.0,  gamma=1.4, k=gamma,
                     mu=1.5e-5, rho=1.2, a0=343.0, Pr=0.707)
    
    
    R = D / 2

    At = π*R^2
    Vt = At * L

    Vr = V/Vt * (sigma + 1/k)

    return TubeSection(L,R,V,Vr,gamma,mu,rho,a0,Pr,sigma,k)
end

"""
`TubeSection(tb; D=1.5, ...)`

Create a new tube section from another one replacing one or more parameters.
The ones not specified as a keyword argument remain unchanged.
"""
function TubeSection(tb::TubeSection{T}; kw...) where {T}

    ks = keys(kw)
    
    D = :D ∈ ks ? T(kw[:D]) : tb.R*2
    L = :L ∈ ks ? T(kw[:L]) : tb.L
    V = :V ∈ ks ? T(kw[:V]) : tb.V
    γ = :gamma ∈ ks ? T(kw[:gamma]) : tb.γ
    μ = :mu ∈ ks ? T(kw[:mu]) : tb.μ
    ρ = :rho ∈ ks ? T(kw[:rho]) : tb.ρ
    a₀ = :a0 ∈ ks ? T(kw[:a0]) : tb.a₀
    Pr = :Pr ∈ ks ? T(kw[:Pr]) : tb.Pr
    σ = :sigma ∈ ks ? T(kw[:sigma]) : tb.σ
    k = :k ∈ ks ? T(kw[:k]) : tb.k
    
    TubeSection(D=D, L=L, V=V, sigma=σ, gamma=γ, k=k, mu=μ, rho=ρ, a0=a₀, Pr=Pr)
    
end

imsqrtim(T) = one(T)*im*sqrt(one(T)*im)

"""
`bessrat1(a)`

Retorna ``J₀(a⋅i√i) / J₂(a⋅i√i)`` 
"""
function bessrat1(a::T) where {T}
    ai = a*imsqrtim(T)
    besselj(2,ai) / besselj0(ai)
end

"""
`bessrat(a)`

Retorna ``J₂(a⋅i√i) / J₀(a⋅i√i)`` mas usa uma aproximação quando `a` é grande.

Quando `a` é grande, o numerador e denominador crescem bastante e existe forte possibilidade de perda
de precisão e, pior, overflow. Então esta funçãoi usa uma aproximação para `a` grande.

"""
function bessrat(a::T) where {T<:AbstractFloat}

    if abs(a) < 100
        return bessrat1(a)
    else
        return 1/(imsqrtim(T) * 2/a - one(T))
    end
        
end


"""
`tubecoefs(tb, ω)`

Return the basic coefficients for computing pressure correction for a given
tube section `tb` and angular frequency `ω`. The angular frequency is in
radians per second **not Hz**.
"""
function tubecoefs(tb::TubeSection{T}, ω) where {T}
    αᵣ =  tb.R * sqrt(tb.ρ*ω/tb.μ)
    s = αᵣ * sqrt(tb.Pr)
    γ = tb.γ

    
    n = 1 / ( one(T) + ((γ-1)/γ) * bessrat(s) )
    
    if ω < 1e-8
        ϕ = zero(Complex{T})
    else
        ϕ = ω/tb.a₀ * sqrt( 1/bessrat(αᵣ) * γ/n )
    end
    
    return αᵣ, ϕ, n
end


"""
`pressratio(tb, f)`

Compute the input output pressure ratio for a tube section at a given frequency (Hz).

The pressure ratio is the ratio between input pressure and output pressure in the
tubing system.

"""
function pressratio(tb::TubeSection{T}, f) where {T}
    ω = T(2π*f)
    α, ϕ, n = tubecoefs(tb, ω)
    ϕL = ϕ*tb.L
    return 1/( cosh(ϕL) + tb.Vr * n * ϕL * sinh(ϕL) )
end

"""
`phaseangle(r::AbstractVector)`

The pressure correction for a given frequency is a complex number so that
both amplitude and phase can be corrected. When plotting the phase
of the correction, it is often difficult to view the results since the
angle of the complex number (argument) is limited to ``-π ≤ angle(z) ≤ π``
so that a continuos change in frequency will result in a discontinuous
change in `angle(z)`. This function calculates the phase angle in sequence and
if a change in phase is angle is large, compensates this by adding or subtracting
2π.

"""
function phaseangle(r::AbstractVector{Complex{T}}) where {T}
    ϕ = angle.(r)
     
    s = zero(T)
    for i = firstindex(r)+1:lastindex(r)
        dϕ = (ϕ[i]+s) - ϕ[i-1]
        if dϕ > 6
            s -= 2π
        elseif dϕ < -6
            s += 2π
        end
        ϕ[i] += s
    end
    return ϕ
end


struct PressureLine{T}
    "`TubeSection`s that make up the pressure line"
    tubes::Vector{TubeSection{T}}
end
Base.broadcastable(tubes::PressureLine) = Ref(tubes)
"""
`PressureLine()`
`PressureLine( (t1, t2, ...) )`
`PressureLine(t1, t2, t3, ...)`
`PressureLine([t1, t2, t2, ...])`

Creates a `PressureLine` object from the different tube sections that make up
the pressure line. If no correction is necessary, a `PressureLine` with no tubes can
be used.

A tube section is a piece of tube with a given diameter and length. At the end of it
there is a volume. In the simplest case, the volume is the pressure chamber of the
pressure sensor.

If tubes with different diameters are connected together, this corresponds to
2 `TubeSection`s connected together where the middle volume is 0.

"""
PressureLine() = PressureLine(TubeSection{Float64}[])

PressureLine(tubes::NTuple{N,TubeSection{T}}) where {N,T} =
    PressureLine(collect(tubes))

PressureLine(tubes::Vararg{TubeSection{T}}) where {T} = PressureLine(collect(tubes))

Base.length(tubes::PressureLine) = length(tubes.tubes)

function Base.show(io::IO, tb::PressureLine)
    print(io, "Pressure line with $(length(tb.tubes)):")
    for t in tb.tubes
        println(io)
        show(io, t)
    end
end


function pressratio(tubes::PressureLine{T}, f) where {T}
    
    if length(tubes) == 0
        return one(Complex{T})
    end
    
    tb = tubes.tubes
    ω = T(2π*f)
    Nt = length(tb)

    # Temos Nt seções. Podemos imaginar que temos Nt+1 seções
    # mas a seção Nt+1 tem comprimento 0, volume 0 e diâmetro 0.
    # Assim r = P_Nt / P_Nt+1 = 1
    # Agora queremos calcular P_Nt / P_Nt-1
    α1, ϕ1, n1 = tubecoefs(tb[end], ω)
    L1 = tb[end].L
    ϕL1 = ϕ1*L1
    R1 = tb[end].R
    jrat_1 = bessrat(α1)
    
    ch1 = cosh(ϕL1); sh1 = sinh(ϕL1)
    r1 = 1 / ( ch1 + tb[end].Vr * n1 * ϕL1 * sh1 )
    rout = r1
    for i in lastindex(tb)-1:-1:firstindex(tb)
        R = tb[i].R
        L = tb[i].L
        α, ϕ, n = tubecoefs(tb[i], ω)
        ϕL = ϕ*L
        jrat = bessrat(α)
        ch = cosh(ϕL)
        sh = sinh(ϕL)

        ra = ch + tb[i].Vr * n * ϕL * sh

        # Outra contribuição:
        rb = ( R1*R1 * ϕ1 * jrat_1  * sh ) / ( R*R * ϕ * jrat *  sh1 ) *
            ( ch1 - r1)

        r = 1 / (ra + rb)

        r1, ϕ1, jrat_1, R1, sh1, ch1 = r, ϕ, jrat, R, sh, ch

        rout = rout * r
        
    end

    return rout
    
end

"""
`amplphase(tube,f)`

Bulds a table with amplitude ratio and phase angle for a given
pressure line `tube` at frequencies specified by the vector `f`.
"""
function amplphase(tube, f)
    r = pressratio.(tube, f)
    ampl = abs.(r)
    phase = phaseangle(r)*180/π

    return [f ampl phase]
end

            


struct PressCorrect{T}
    "`pressratio` for each frequency of a tubing system"
    r::Matrix{Complex{T}}
    "Sample rate of the measured pressure"
    fs::T
    "Number of pressure samples"
    n::Int
    ""
    indices::Vector{Int}
end
Base.broadcastable(corr::PressCorrect) = Ref(tb)


"""

"""
function presscorrect(tube::PressureLine, P, fs)
    corr = PressCorrect(size(P,1), fs, tube)
    return presscorrect(corr, P)
end

presscorrect(tube::TubeSection, P, fs) = presscorrect(PressureLine(tube), P, fs)

(tube::PressureLine)(P,fs) = presscorrect(tube, P, fs)
(tube::TubeSection)(P,fs) = presscorrect(tube, P, fs)


function PressCorrect(n, fs, tube::PressureLine{T}) where{T}

    f = rfftfreq(n, fs) 
    r = pressratio.(tube, f)

    PressCorrect(hcat(r), fs, n, Int[])
    
end


function PressCorrect(n, fs, indices::AbstractVector{Int},
                            tubes::AbstractVector{PressureLine{T}}) where {T}
    
    f = rfftfreq(n,fs)
    nf = length(f)
    nt = length(tubes)
    
    rall = zeros(Complex{T}, nf, nt)

    for (r,t) in zip(eachcol(rall), tubes)
        r .= pressratio.(t, f)
    end

    uind = extrema(indices)  # Verificar se os índices estão compatíveis com `tubes`

    if uind[1] < 0 || uind[2] > nt
        error("Indices make reference to the the tubes")
    end

    PressCorrect(rall, fs, n, indices)
    
end

(corr::PressCorrect)(P) = presscorrect(corr, P)


function presscorrect(corr::PressCorrect, P)
    if length(corr.indices) > 0
        # This case considers multiple pressure lines
        # Let's go to the more complicated function for this
        return presscorrect_mult(corr, P)
    end

    if size(P,1) !=  corr.n
        error("`PressCorrect` has n = $(corr.n) but P has $(size(P,1)) samples!")
    end

    PF = rfft(P, 1)
    for pf in eachcol(PF)
        pf .*= 1 ./ corr.r[:,1]
    end

    return irfft(PF, corr.n, 1)
end


function presscorrect_mult(corr::PressCorrect, P)

    if length(corr.indices) != size(P,2)
        error("The number of indices should the same as the number of columns of the pressure array!")
    end

    if size(P,1) !=  corr.n
        error("`PressCorrect` has n = $(corr.n) but P has $(size(P,1)) samples!")
    end
    

    PF = rfft(P, 1)
    for (k, pf) in enumerate(eachcol(PF))
        pf .*= 1 ./ corr.r[:,corr.indices[k]]
    end

    return irfft(PF, corr.n, 1)
    
end

#=
# Funções para plotar as linhas de pressão
function plot_pressline!(fig1, f; pr...)
    ks = keys(pr)
	
    fig = fig1[1,1] = GridLayout()
    ax = Axis(fig[1,1], ylabel="Amplitude")
    hidexdecorations!(ax, grid=false)

    r = [pressratio.(p, f) for (n,p) in pr]
    labels = [string(n) for (n,p) in pr]
    
    for (lab,r) in zip(labels,r)
        ampl = abs.(r)
	lines!(ax, f, ampl, label=lab)
    end
    
    ax1 = Axis(fig[2,1], xlabel="Frequência (Hz)", ylabel="Ângulo de fase (°)")
    linkxaxes!(ax1, ax)
    rowgap!(fig, 0)
    
    for (lab,r) in zip(labels,r)
        phase = phaseangle(r)*180/π
	lines!(ax1, f, phase, label=lab)
    end
    axislegend(ax1, position=:rt)	
    fig1
end

plot_pressline(f ; pr...) = plot_pressline!(Figure(), f; pr...)

=#
