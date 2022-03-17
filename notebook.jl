### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 1fbab730-f787-4b52-9710-52ef5bb24d21
using DrWatson

# ╔═╡ 69c07bd8-6b8c-4e5c-b9ce-bfa98961126e
DrWatson.@quickactivate "TCellMemoryLineage"

# ╔═╡ 148d1c62-077f-4b14-b539-b0e95e1a4a93
begin
	using DifferentialEquations
	using Plots
end

# ╔═╡ 8c2700e1-9b0a-484b-82cb-2d8fa6e44056
begin
	pathname = @__DIR__
	foldername = basename(pathname)
end

# ╔═╡ a5541b99-4ace-47a8-9347-9e189e9fcb67
md"# $(foldername)"

# ╔═╡ aa68c3c8-f333-4c8d-8279-1a693b7f6477
md"## Libraries"

# ╔═╡ ff843cb4-2247-4903-b3bf-73168f726b06
md"## Parameters"

# ╔═╡ 90f9ad3d-2803-4766-8b25-fcfe69af97da
md"Empirical labelling parameters"

# ╔═╡ bc150210-2735-4f53-95a8-673b4cba359a
fr=0.026

# ╔═╡ 34645434-302b-4516-9137-898518bbdacc
delta=	0.1

# ╔═╡ 559980c9-2bd0-4a32-82f4-6fa9535a5ac2
tau= 5.0

# ╔═╡ d40360c1-2652-48b9-bafe-e835c7c94237
beta=0.002

# ╔═╡ b2f3506f-7410-4444-b5be-bd6c1f2d2e58
bw = 5.0

# ╔═╡ 3c377efc-375d-4faa-a196-368b74f33748
label_p_v = [tau, fr, delta, beta, bw]

# ╔═╡ 18063995-38f5-44f6-b456-fa005b7fd39d
md"Kinetic parameters"

# ╔═╡ 7dd830a9-65f3-4c64-8602-ac142d1ebae3
p_a = 0.015

# ╔═╡ 9471f592-6959-47d7-9298-77e35f9a72bb
d_a=0.05

# ╔═╡ 164b0544-b33a-4d62-9eb1-9229f55d7eed
p_vec = [p_a, d_a]

# ╔═╡ b36d44d2-972c-48bb-b2e2-e134df25b0a5
md"## ODE model"

# ╔═╡ ddbc3842-8a42-47ee-8070-98fd63d9e2a8
tspan = (0.0, 50.0)

# ╔═╡ d76aaf92-6401-427e-bff3-0bbba8352e33
u0 = [0.0]

# ╔═╡ b332cb33-c8d4-45bd-b92f-f7e8027839cb
function _ODEmodel(du,u,p,t, label_p)
    tau, fr, delta, beta, bw = label_p
   	p_a, d_a = p

    U = ifelse(t >= tau, (fr*(1-exp(-delta*tau))+beta*exp(-delta*tau))*exp(-delta*(t-tau)), fr*(1-exp(-delta*t))+beta*exp(-delta*t))

    du[1] = p_a * bw * U - d_a *u[1]
end

# ╔═╡ 0688ad66-e327-4189-a89a-6e22a6e0726b
ODE_model(du,u,p,t) = _ODEmodel(du,u,p,t, label_p_v)

# ╔═╡ d1b8d641-2fd4-4698-acdb-08930bf049be
ode_prob = ODEProblem(ODE_model, u0, tspan, p_vec)

# ╔═╡ 78c597dc-8caa-469f-8179-fdc7aa767bec
md"## Solve ODE model"

# ╔═╡ b9ba7ad1-6379-4393-98ce-72b190b07b7e
sol = solve(ode_prob)

# ╔═╡ 6ba2c6d3-abf3-48b1-bd5e-e5021dfba49b
md"## Plot simulation"

# ╔═╡ ec205b54-d3f6-4b2e-bc42-31c5317cf318
plot(sol, ylabel="label enrichment", label= "A")

# ╔═╡ Cell order:
# ╠═a5541b99-4ace-47a8-9347-9e189e9fcb67
# ╠═8c2700e1-9b0a-484b-82cb-2d8fa6e44056
# ╠═aa68c3c8-f333-4c8d-8279-1a693b7f6477
# ╠═1fbab730-f787-4b52-9710-52ef5bb24d21
# ╠═69c07bd8-6b8c-4e5c-b9ce-bfa98961126e
# ╠═148d1c62-077f-4b14-b539-b0e95e1a4a93
# ╠═ff843cb4-2247-4903-b3bf-73168f726b06
# ╠═90f9ad3d-2803-4766-8b25-fcfe69af97da
# ╠═bc150210-2735-4f53-95a8-673b4cba359a
# ╠═34645434-302b-4516-9137-898518bbdacc
# ╠═d40360c1-2652-48b9-bafe-e835c7c94237
# ╠═b2f3506f-7410-4444-b5be-bd6c1f2d2e58
# ╠═3c377efc-375d-4faa-a196-368b74f33748
# ╠═18063995-38f5-44f6-b456-fa005b7fd39d
# ╠═164b0544-b33a-4d62-9eb1-9229f55d7eed
# ╠═b36d44d2-972c-48bb-b2e2-e134df25b0a5
# ╠═ddbc3842-8a42-47ee-8070-98fd63d9e2a8
# ╠═d76aaf92-6401-427e-bff3-0bbba8352e33
# ╠═b332cb33-c8d4-45bd-b92f-f7e8027839cb
# ╠═0688ad66-e327-4189-a89a-6e22a6e0726b
# ╠═d1b8d641-2fd4-4698-acdb-08930bf049be
# ╠═78c597dc-8caa-469f-8179-fdc7aa767bec
# ╠═b9ba7ad1-6379-4393-98ce-72b190b07b7e
# ╠═6ba2c6d3-abf3-48b1-bd5e-e5021dfba49b
# ╠═7dd830a9-65f3-4c64-8602-ac142d1ebae3
# ╠═9471f592-6959-47d7-9298-77e35f9a72bb
# ╠═559980c9-2bd0-4a32-82f4-6fa9535a5ac2
# ╠═ec205b54-d3f6-4b2e-bc42-31c5317cf318
