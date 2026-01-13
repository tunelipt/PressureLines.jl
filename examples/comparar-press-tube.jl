### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ f5c1100c-319d-4472-accd-c8ecb776baf2
begin
	import Pkg
	Pkg.activate("..")
end

# ╔═╡ 20af59ce-f007-4050-ad95-8e5ed8aa8f57
begin
	using Revise
	using CairoMakie
	using FFTW
	using SpecialFunctions
	using Statistics
	using HDF5
end

# ╔═╡ b2d2bfe1-f1a9-4124-aa3c-a9c5cf33e7d0
using EnsaioCp

# ╔═╡ f1b98ec2-7670-11ef-31fa-e551e7cf4272
md"""
# Comparação da correção de pressão com dados publicados

A correção da pressão por influência dos tubos usados na medição é implementada usando a abordagem descrita no artigo

''Opttimization of dynamic-pressure-measurement systems. I. Single point measurements''. J. D. Holmes e R. E. Lewis, Journal of Wind Engineering and Industrial Aerodynamics, 25, p249-273, 1987.


Muitas funções e códigos que carregam os pacotes e preparam o ambiente estão no final deste notebook.
"""

# ╔═╡ 6716d060-d63f-436c-8f79-c1001e3821cc
md"""
## Carregar os dados

Eu digitalizei algumas curvas do artigo do Holmes de 87. Isso serve para comparar.

Um detalhe importante: as correções são bem sensíveis às propriedades termodinâmicas usadas no cálculo. 
"""

# ╔═╡ 1b969989-4051-408b-b3bf-b9ece8bffed6
md"""
## Comparação dos casos simples

### D = 1.5 mm, L = 500.0 e V = 257 mm³
"""

# ╔═╡ 98e57ba7-6ae8-4a67-a554-c4c595eb42eb
md"""
### D = 1.0 mm, L = 500.0 e V = 130 mm³
"""

# ╔═╡ 65d49156-d45b-46fb-93f7-b98afc3b1d47
md"""
## Casos compostos
## L₁=305 mm, D₁ = 1.5 mm, L₂ = 195 mm, D₂ = 0.5 mm, V₂ = 100 mm³

"""

# ╔═╡ 5c070315-bf55-4d82-beed-f326d96ed85c
md"""
## L₁=270 mm, D₁ = 1.5 mm, L₂ = 230 mm, D₂ = 0.55 mm, V₂ = 100 mm³

"""

# ╔═╡ e18c62e3-af90-4f9d-bc03-afd1df65107b
md"""
## L₁=344 mm, D₁ = 1.5 mm, L₂ = 156 mm, D₂ = 0.45 mm, V₂ = 100 mm³

"""

# ╔═╡ 70af2d77-15c0-4b37-b383-5bed3fd21a33
md"""
## Discussão e conclusões

 * Os casos simples caem como uma luva
 * Existem pequenas diferenças nos casos compostos
 * Não sabemos as condições de cálculo da correção então sabendo a provável cidade onde foram realizados os testes (Melbourne) e ajustando as condições ambientais obteve-se as curvas acima. 
 * Os casos compostos apresentam um termo complicado que pode ter algum errinho. Mais testes são necessários
"""

# ╔═╡ de665304-1b6a-4970-ab0b-5d58519268ad
begin
	const Ma = 28.9647
	const Ru = 8314.32 
	const Ra = Ru/Ma
end

# ╔═╡ 53ebc828-2c7e-4290-bb2b-2ccb1fc793b5
specmass(T,P=101.325) = P*1000 / (Ra*(T+273.15))

# ╔═╡ cb5216ba-e7ad-4d8d-86a0-e040c85ef96d
viscosity(T) = 1.716e-5 * ((T+273.15)/273.15)^1.5 * (273.15+110.4)/(T+273.15+110.4)

# ╔═╡ 864dfb41-d5c0-4ed1-9c85-9940d04f9290
begin
	const Ta = 15
	const Pa = 101.0
	const γ = 1.4
	const a₀ = sqrt(γ*Ra*(Ta+273.15))

	ρ = specmass(Ta, Pa)
	μ = viscosity(Ta)
end

# ╔═╡ 5d86bdee-7b76-4365-8b78-0299559dddc0
function compare_methods(fig1; pr...)
	ks = keys(pr)
	
	fig = fig1[1,1] = GridLayout()
	ax = Axis(fig[1,1], ylabel="Amplitude")
	hidexdecorations!(ax, grid=false)

	for (n,p) in pr
		lines!(ax, p[:,1], p[:,2], label=string(n))
	end

	ax1 = Axis(fig[2,1], xlabel="Frequency (Hz)")
	linkxaxes!(ax1, ax)
	rowgap!(fig, 0)
		
	for (n,p) in pr
		lines!(ax1, p[:,1], p[:,3], label=string(n))
	end
	axislegend(ax1, position=:rt)	
	fig1

end

# ╔═╡ b255bd19-1e1a-410f-9aa4-342728bb423a
function pratiotable(f, t)
	r = pressratio.(t, f)
	ampl = abs.(r)
	phase = phaseangle(r) * 180/π
	return [f ampl phase]
end

# ╔═╡ 822c6413-7e6d-4bdb-ba9e-ccd7e2f00681
function load_paper_data(fname, path; kw...)
	prat, D, L, V = h5open(fname, "r") do h
		prat = read(h["$path/press_ratio"])
		D = read(h["$path/D"])
		L = read(h["$path/L"])
		V = read(h["$path/V"])
		prat, D, L, V
	end
	f = prat[:,1]	
	ampl = prat[:,2]
	phase = prat[:,3] * π/180
	r = ampl .* cis.(phase)

	tubes = [TubeSection(D=D[i], L=L[i], V=V[i]; kw...) for i in eachindex(L)]
	return (f=f, r=r, tubes=tubes, prat=prat)
	
end

# ╔═╡ 0037891f-660c-4b6d-a4f0-daea8904b8cd
# Ler os dados dos gráficos do artigo do Holmes
begin
	# Figura 6 - tubos rígidos
	f06r = load_paper_data("press-tube-data.h5", "FIG06/rigid",rho=ρ,mu=μ,a0=a₀)
	# Figura 6 - tubos flexíveis
	f06f =  load_paper_data("press-tube-data.h5", "FIG06/flexible",rho=ρ,mu=μ,a0=a₀)
	# Figura 7 - tubos rígidos
	f07r = load_paper_data("press-tube-data.h5", "FIG07/rigid", rho=ρ, mu=μ, a0=a₀)
	# Figura 7 - tubos flexíveis
	f07f = load_paper_data("press-tube-data.h5", "FIG07/flexible", rho=ρ, mu=μ, a0=a₀)
	
	# Figura 8 - L = 195mm
	f08L195 = load_paper_data("press-tube-data.h5", "FIG08/L195",rho=ρ,mu=μ,a0=a₀)
	# Figura 8 - L = 230mm
	f08L230 = load_paper_data("press-tube-data.h5", "FIG08/L230",rho=ρ,mu=μ,a0=a₀)
	# Figura 8 - L = 156mm
	f08L156 = load_paper_data("press-tube-data.h5", "FIG08/L156",rho=ρ,mu=μ,a0=a₀)
end;

# ╔═╡ 2b04950e-d48c-48b4-8a53-565c5183f2ab
begin
	# Vamos criar as linhas de pressão correspondentes a cada caso.
	# Isso é a minha implementação do cálculo
	t06r = PressureLine(f06r.tubes)
	t06f = PressureLine(f06f.tubes)
	t07r = PressureLine(f07r.tubes)
	t07f = PressureLine(f07f.tubes)
	f1 = f06r.prat[:,1]

	t195 = PressureLine(f08L195.tubes)
	t230 = PressureLine(f08L230.tubes)
	t156 = PressureLine(f08L156.tubes)
	f2 = f08L195.prat[:,1]

end;
	

# ╔═╡ a4d11255-7b4f-4ca3-93fd-fd44569314b8
compare_methods(Figure(), Rigido=f06r.prat, EnsaioCp=pratiotable(f1, t06r))

# ╔═╡ 71a104db-a3f7-4a26-aae3-4b38e6aa1c53
compare_methods(Figure(), Rigido=f06r.prat, EnsaioCp=pratiotable(f1, t06f))

# ╔═╡ b341d7ce-12e8-4d18-9ff5-96c8486bcf04
compare_methods(Figure(), Rigido=f08L195.prat, EnsaioCp=pratiotable(f2, t195))

# ╔═╡ 4a1f9deb-ff44-4ecb-bdcf-79b777b8f344
compare_methods(Figure(), Rigido=f08L230.prat, EnsaioCp=pratiotable(f2, t230))

# ╔═╡ e7b88480-4c2f-49be-be54-63c5672b954a
compare_methods(Figure(), Rigido=f08L156.prat, EnsaioCp=pratiotable(f2, t156))

# ╔═╡ Cell order:
# ╟─f1b98ec2-7670-11ef-31fa-e551e7cf4272
# ╟─6716d060-d63f-436c-8f79-c1001e3821cc
# ╟─0037891f-660c-4b6d-a4f0-daea8904b8cd
# ╟─2b04950e-d48c-48b4-8a53-565c5183f2ab
# ╟─1b969989-4051-408b-b3bf-b9ece8bffed6
# ╟─a4d11255-7b4f-4ca3-93fd-fd44569314b8
# ╟─98e57ba7-6ae8-4a67-a554-c4c595eb42eb
# ╟─71a104db-a3f7-4a26-aae3-4b38e6aa1c53
# ╟─65d49156-d45b-46fb-93f7-b98afc3b1d47
# ╟─b341d7ce-12e8-4d18-9ff5-96c8486bcf04
# ╟─5c070315-bf55-4d82-beed-f326d96ed85c
# ╟─4a1f9deb-ff44-4ecb-bdcf-79b777b8f344
# ╟─e18c62e3-af90-4f9d-bc03-afd1df65107b
# ╟─e7b88480-4c2f-49be-be54-63c5672b954a
# ╟─70af2d77-15c0-4b37-b383-5bed3fd21a33
# ╠═f5c1100c-319d-4472-accd-c8ecb776baf2
# ╠═20af59ce-f007-4050-ad95-8e5ed8aa8f57
# ╠═b2d2bfe1-f1a9-4124-aa3c-a9c5cf33e7d0
# ╠═de665304-1b6a-4970-ab0b-5d58519268ad
# ╠═864dfb41-d5c0-4ed1-9c85-9940d04f9290
# ╠═53ebc828-2c7e-4290-bb2b-2ccb1fc793b5
# ╠═cb5216ba-e7ad-4d8d-86a0-e040c85ef96d
# ╠═5d86bdee-7b76-4365-8b78-0299559dddc0
# ╠═b255bd19-1e1a-410f-9aa4-342728bb423a
# ╠═822c6413-7e6d-4bdb-ba9e-ccd7e2f00681
