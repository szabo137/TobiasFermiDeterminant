using LinearAlgebra
using Richardson
using SpecialFunctions
using BenchmarkTools


function _beta(rs::Real,theta::Real)
    return ( 9.0*pi/4 )^(-2.0/3.0) * 2.0*rs^2/theta
end

function _n(rs::Real)
    return 1.0 / ( 4.0/3.0 * pi * rs^3 )
end


function _L(rs::Real,N::Real)
    return ( 4.0/3.0*pi*rs^3*N*2)^(1.0/3.0 )
end

function Z_sum_index2(beta::Real, L::Real)

 	E = pi^2*2.0/(L^2)

	ans, err = extrapolate(1, x0=Inf) do N
            one(beta) + sum(n -> 2.0 * exp( -beta * (E * n^2 ) ), 1:Int(N))
        end

	return ans^3
end

function _fermi_mat(n,beta,L)
    Fermi = zeros(typeof(beta),n,n)
    @inbounds for idx in CartesianIndices(Fermi)

        if idx[1] == idx[2]+1
            Fermi[idx] = float(idx[1])
        end

        if idx[1] <= idx[2]
            multi = idx[2] + 1 - idx[1]

            Fermi[idx] = Z_sum_index2(multi*beta, L)
        end
    end

    return Fermi
end

function _fermi_logdet(n,beta,L)
    val,sig = logabsdet(UpperHessenberg(_fermi_mat(n,beta,L)))
    return val
end

function UEG_N(n,beta,L)
    fermi_logdet= _fermi_logdet(n,beta,L)
    gamma_logdet = loggamma(n+1)
    return (-fermi_logdet+gamma_logdet)/beta
end

function main()
    println("Calculation: UEG_N (Fermi)")
    N = 1000
    @show N
	rs = 2
    @show rs
	theta = 1.0
    @show theta
	beta = _beta(rs,theta)
    @show beta
	L = _L(rs,N)
    @show L
        
    res = UEG_N(N,beta,L)
    @show res


    return nothing
end

function bench()
    println("Benchmark: UEG_N (Fermi)")
    N = 1000
    @show N
	rs = 2
    @show rs
	theta = 1.0
    @show theta
	beta = _beta(rs,theta)
    @show beta
	L = _L(rs,N)
    @show L

    bench_result = @benchmark UEG_N($N,$beta,$L)

    display(bench_result)
    
    return nothing
end
