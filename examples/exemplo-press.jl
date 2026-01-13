### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 3e321160-12f9-4b6c-a4ec-fbf0e9cbc485
begin
	import Pkg
	Pkg.activate("..")
end

# ╔═╡ d2e1f483-68b0-4868-9693-7b728109c8c9
begin
	using Revise
	using CairoMakie
	using FFTW
	using Statistics
	using HDF5
	using EnsaioCp
end

# ╔═╡ 4509b682-76c0-11ef-2b72-45853055e84f
md"""
# Exemplo de correção de pressão
"""

# ╔═╡ 17f50762-5cb1-44a9-a351-2359057c2114
md"""
## Caso sintético

Vamos gerar funções conhecidas com frequências adequadas.
"""

# ╔═╡ 081e663e-ae0a-4879-9b3d-b2b7b04be4e0


# ╔═╡ 39b49862-2d55-4d9e-89f3-c88bb2ded8e2
begin
	fs = 500.0 # Taxa de amostragem
	dt = 1/fs
	N = 500 # 1s de amostragem
	f = rfftfreq(N, fs)

	f1 = f[101] # Primeira freqencua
	f2 = f[121] # Segunda frequencia
	f3 = f[161] 

	t = range(0.0, step=dt, length=N+1)[1:end-1] # Olha o detalhe do tempo...
end;
	

# ╔═╡ d0fde5ac-20b2-4c92-b4f3-3391d30d223e
begin
	x0 = 0 * t .+ 1  # Função constante
	x1 = sin.(2π*f1.*t)
	x2 = sin.(2π*f2.*t)
	x3 = sin.(2π*f3.*t)
	x4 = x0 + x1 + 2x2 + 3x3
end

# ╔═╡ ab2944e3-f701-42fd-96cf-ccf4a9660356
md"""
## Vamos gerar um tubo
"""

# ╔═╡ fbc31556-cd52-465a-abcf-3000c1eb0beb
tubo = PressureLine(TubeSection(D=1.0e-3, L=0.5, V=100e-9, rho=1.1, mu=1.8e-5))

# ╔═╡ 8491dc3a-f0e9-4499-a9c7-b3f093c05be3
pressratio(tubo, 200)

# ╔═╡ ecdca582-77c8-44ed-bf7d-c9eb34dcf7db
md"""
## Aplicar as equações do tubo

Começando pelo caso trivial: função com valor constante!
"""

# ╔═╡ 50fd7d08-b455-4d3d-b26d-0cb5b4508f91
P0 = tubo(x0, fs)

# ╔═╡ 8d22cddc-5052-436e-bc1d-8ca3e637c7f0
extrema(P0 - x0)

# ╔═╡ 47f52d19-cca0-4992-8298-271ca1fee838
# Caso mais complicado

# ╔═╡ b2a38828-95b3-4689-9ee7-f4de3b9a5593
P1 = tubo(x1, fs)

# ╔═╡ 0deec41d-68d3-4089-b532-212c616547f5
# Outro jeito de fazer isso é reaproveitando os cálculos
corr = PressCorrect(N, fs, tubo)

# ╔═╡ adfd7050-26ae-4c89-92a7-d0b65275877d
# Agora `corr` pode ser aplicado sucessivamente
P1b = corr(x1)

# ╔═╡ 6fbbd11f-0b56-4bfe-a011-d8714788abf3
# Tá certinho, como esperado
extrema(P1b - P1)

# ╔═╡ ce725342-96b1-4bb4-972e-23952ec93e4e
let fig=Figure()
	ax = Axis(fig[1,1], xlabel="Tempo (s)", ylabel="Pressão")
	xlims!(ax, (0.0, 0.1))
	lines!(ax, t, x1, label="Sinal original")
	lines!(ax, t, P1, label="Sinal corrigido")
	axislegend(ax, position=:rt)
	fig
end

# ╔═╡ 27ed5c02-7d28-4346-add9-d4b674680c67
# Comparando a amplitude do sinal com a correção
std(P1) / std(x1), 1/abs(corr.r[101,1])

# ╔═╡ 38e946e4-0bf9-468a-9cb3-704ce6b94900
begin
	P2 = corr(x2)
	P3 = corr(x3)
	P4 = corr(x4)
end

# ╔═╡ 1469aced-9865-4689-bfba-0aa434cda374
std(P2) / std(x2), 1/abs(corr.r[121,1])

# ╔═╡ ac68c125-1a56-40df-82b1-af4ce5a7031a
std(P3) / std(x3), 1/abs(corr.r[161,1])

# ╔═╡ 0789d3dc-4da2-4ce6-a058-b75f7596b52e
P4b = P0 + P1 + 2P2 + 3P3 

# ╔═╡ ba6ff409-f85b-4cd3-ba21-615c8d830781
let fig=Figure()
	ax = Axis(fig[1,1], xlabel="Tempo (s)", ylabel="Pressão")
	xlims!(ax, (0.5, 0.51))
	lines!(ax, t, x4, label="Sinal original")
	lines!(ax, t, P4, label="Sinal corrigido")
	lines!(ax, t, P4b, label="Sinal corrigido montado", linestyle=:dash, linewidth=2)
	axislegend(ax, position=:rt)
	fig
end

# ╔═╡ f0df8782-b50e-4ea8-be55-9709b227e5d6
md"""
## Usando os dados com um bloco de pressões

A estrutura de dados `PressCorrect` na forma de `corr` pode ser aplicado a uma matriz de dados onde cada coluna corresponde a uma pressão.
"""

# ╔═╡ 3f1feeb9-88b1-42be-9b27-2ff836a03635
# Vamos criar a matriz de dados
X = [x0 x1 x2 x3 x4]

# ╔═╡ a0684153-6664-41ed-bd07-f86b6c546954
P = corr(X)

# ╔═╡ dc152d7b-b337-4351-97cb-e34393b1a0eb
extrema(P - [P0 P1 P2 P3 P4]) # Batendo tudo

# ╔═╡ 9fe93266-21d1-4569-8438-58a99f479735
# Também podemos usar o tubo diretamente
Pb = tubo(X, fs)

# ╔═╡ c866fa56-1d3f-4a65-8261-9e01f1a569d3
extrema(Pb - P)

# ╔═╡ 4eb24f02-bb04-4335-a22e-81cc866b334b
md"""
## Casos com mais de uma mangueira

Também contemplamos esta situação
"""


# ╔═╡ 4311b308-ed1f-4420-b258-0ed33aadaf41
# Vamos fazer um segundo conjunto de tubulações que são mais compridas
tubo2 = PressureLine(TubeSection(D=1.0e-3, L=0.9, V=100e-9, rho=1.1, mu=1.8e-5))

# ╔═╡ d519e4de-c876-47f0-82eb-48b7368c2ddc
# Cada coluna das pressões corresponde a uma linha de pressão.
# Agora especificamos qual
indices = [1, 2, 1, 2, 2]  # Canais 1, 3 - tubo. Canais 2,4,5 - tubo2

# ╔═╡ 35e27df6-783a-4a5c-abce-1c407265ebe0
corr_bloco = PressCorrect(N, fs, indices, [tubo, tubo2])

# ╔═╡ ebde5e48-2d5c-4c07-9c85-8f39dfa99089
# Vamos usar `corr_bloco`
Y = corr_bloco(X)

# ╔═╡ f53a1654-4f38-4e93-a95d-0e1ba61712a4
# Vamos aplicar os tubos individualmente
begin
	Y0 = tubo(x0, fs)
	Y1 = tubo2(x1, fs)
	Y2 = tubo(x2, fs)
	Y3 = tubo2(x3, fs)
	Y4 = tubo2(x4, fs)
end

# ╔═╡ 09a61078-9426-4701-b7b2-0e3b5178cfb1
extrema(Y - [Y0 Y1 Y2 Y3 Y4])

# ╔═╡ 10df6479-fdc7-4825-8138-5a47c38c798a
# Vamos ver nos casos do tubo 2, se estamos calculando a amplitude corretamente:
std(Y1) / std(x1), 1/abs(corr_bloco.r[101,2])

# ╔═╡ a5ca1552-fcef-4146-b28c-7ba43d24a494
std(Y3) / std(x3), 1/abs(corr_bloco.r[161,2])

# ╔═╡ 08c0533f-8c1a-42d2-8812-d4e41f2d9d54


# ╔═╡ 43732ab1-d3f5-40cb-b334-f43d467ddd51
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

# ╔═╡ 62145dec-2b13-463b-a5bf-f2cad2188cad
# Vamos ver as características do tubo:
compare_methods(Figure(), L500=amplphase(tubo, f))	

# ╔═╡ f8b9e1a9-9a9a-4acf-a9c2-a43e92041346
# Vamos ver o que rola
compare_methods(Figure(), L900=amplphase(tubo2, f))	

# ╔═╡ e785928e-d788-4a9c-a1ed-1beb2751e25b


# ╔═╡ Cell order:
# ╠═4509b682-76c0-11ef-2b72-45853055e84f
# ╠═3e321160-12f9-4b6c-a4ec-fbf0e9cbc485
# ╠═d2e1f483-68b0-4868-9693-7b728109c8c9
# ╠═17f50762-5cb1-44a9-a351-2359057c2114
# ╠═081e663e-ae0a-4879-9b3d-b2b7b04be4e0
# ╠═39b49862-2d55-4d9e-89f3-c88bb2ded8e2
# ╠═d0fde5ac-20b2-4c92-b4f3-3391d30d223e
# ╠═ab2944e3-f701-42fd-96cf-ccf4a9660356
# ╠═fbc31556-cd52-465a-abcf-3000c1eb0beb
# ╠═62145dec-2b13-463b-a5bf-f2cad2188cad
# ╠═8491dc3a-f0e9-4499-a9c7-b3f093c05be3
# ╠═ecdca582-77c8-44ed-bf7d-c9eb34dcf7db
# ╠═50fd7d08-b455-4d3d-b26d-0cb5b4508f91
# ╠═8d22cddc-5052-436e-bc1d-8ca3e637c7f0
# ╠═47f52d19-cca0-4992-8298-271ca1fee838
# ╠═b2a38828-95b3-4689-9ee7-f4de3b9a5593
# ╠═0deec41d-68d3-4089-b532-212c616547f5
# ╠═adfd7050-26ae-4c89-92a7-d0b65275877d
# ╠═6fbbd11f-0b56-4bfe-a011-d8714788abf3
# ╠═ce725342-96b1-4bb4-972e-23952ec93e4e
# ╠═27ed5c02-7d28-4346-add9-d4b674680c67
# ╠═38e946e4-0bf9-468a-9cb3-704ce6b94900
# ╠═1469aced-9865-4689-bfba-0aa434cda374
# ╠═ac68c125-1a56-40df-82b1-af4ce5a7031a
# ╠═0789d3dc-4da2-4ce6-a058-b75f7596b52e
# ╟─ba6ff409-f85b-4cd3-ba21-615c8d830781
# ╠═f0df8782-b50e-4ea8-be55-9709b227e5d6
# ╠═3f1feeb9-88b1-42be-9b27-2ff836a03635
# ╠═a0684153-6664-41ed-bd07-f86b6c546954
# ╠═dc152d7b-b337-4351-97cb-e34393b1a0eb
# ╠═9fe93266-21d1-4569-8438-58a99f479735
# ╠═c866fa56-1d3f-4a65-8261-9e01f1a569d3
# ╟─4eb24f02-bb04-4335-a22e-81cc866b334b
# ╠═4311b308-ed1f-4420-b258-0ed33aadaf41
# ╠═f8b9e1a9-9a9a-4acf-a9c2-a43e92041346
# ╠═d519e4de-c876-47f0-82eb-48b7368c2ddc
# ╠═35e27df6-783a-4a5c-abce-1c407265ebe0
# ╠═ebde5e48-2d5c-4c07-9c85-8f39dfa99089
# ╠═f53a1654-4f38-4e93-a95d-0e1ba61712a4
# ╠═09a61078-9426-4701-b7b2-0e3b5178cfb1
# ╠═10df6479-fdc7-4825-8138-5a47c38c798a
# ╠═a5ca1552-fcef-4146-b28c-7ba43d24a494
# ╠═08c0533f-8c1a-42d2-8812-d4e41f2d9d54
# ╠═43732ab1-d3f5-40cb-b334-f43d467ddd51
# ╠═e785928e-d788-4a9c-a1ed-1beb2751e25b
