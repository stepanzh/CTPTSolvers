struct StabilityResult{T}
    converged::Bool
    molfrac::Vector{T}
    tpd::T
end

stability_success(molfrac, tpd) = StabilityResult(true, molfrac, tpd)
stability_fail(n::Int) = StabilityResult(false, fill(NaN, n), NaN)


function stability(
    mixture::Mixture,
    nmol,
    pressure,
    RT,
    Xinit,
    basephase::Symbol,
    testphase::Symbol,
)
    @assert basephase in (:gas, :liquid) "Неверная фаза (:gas, :liquid)"
    @assert testphase in (:gas, :liquid) "Неверная фаза (:gas, :liquid)"

    molfrac = nmol / sum(nmol)

    try
        # Метод простой итерации
        # Поиск грубого приближения корня для метода Ньютона
        base_term = (
            log.(molfrac)
            .+ log_fugacity_coef(mixture, molfrac, pressure, RT, basephase)
        )
        logX, logX_prev = log.(Xinit), fill(NaN, ncomponents(mixture))

        for iter in 1:100
            logX, logX_prev = logX_prev, logX

            X = exp.(logX_prev)
            molfrac_test = X / sum(X)
            log_φ = log_fugacity_coef(mixture, molfrac_test, pressure, RT, testphase)

            logX = base_term - log_φ

            norm(logX - logX_prev, 2) ≤ 1e-6 && break
        end

        Xrude = exp.(logX)

        # Метод Ньютона
        # Подготовить систему нелинейных уравнений
        target_function = __create_target(
            mixture,
            molfrac,
            pressure,
            RT,
            basephase,
            testphase,
        )
        # Начальное приближение
        # Якобиан
        jacobian(X) = ForwardDiff.jacobian(target_function, X)

        # Решить систему
        Xstationary = newtonsys(target_function, Xrude, jacobian;
            maxiter=50, xtol=1e-6, ftol=1e-6
        )

        # Вернуть ответ
        molfrac_test = Xstationary / sum(Xstationary)

        tpd = let x = molfrac_test, z = molfrac
            ln_φ_test = log_fugacity_coef(mixture, x, pressure, RT, testphase)
            ln_φ_base = log_fugacity_coef(mixture, z, pressure, RT, basephase)
            δμ = @. log(x) + ln_φ_test - (log(z) + ln_φ_base)
            RT * dot(x, δμ)
        end

        return stability_success(molfrac_test, tpd)
    catch e
        @warn e
        return stability_fail(ncomponents(mixture))
    end
end

"Создаёт систему уравнений."
function __create_target(mixture, molfrac_base, pressure, RT, basephase, testphase)
    base_term = (
        log.(molfrac_base)
        .+ log_fugacity_coef(mixture, molfrac_base, pressure, RT, basephase)
    )

    function closure(X)
        # ln Xᵢ + ln φᵢ - hᵢ
        molfrac = X / sum(X)
        ln_φ = log_fugacity_coef(mixture, molfrac, pressure, RT, testphase)
        return log.(X) .+ ln_φ .- base_term
    end

    return closure
end
