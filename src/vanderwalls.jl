struct VanDerWaalsComponent{T} <: Component
    a::T
    b::T
    critical_pressure::T
    critical_thermal_energy::T  # R * temperature
end

struct VanDerWaalsMixture{C,M} <: Mixture
    components::C
    aij::M
end

function VanDerWaalsMixture(components)
    a = getfield.(components, :a)
    aij = sqrt.(a * a')
    sc = StructVector(components)

    return VanDerWaalsMixture(sc, aij)
end

ncomponents(m::VanDerWaalsMixture) = length(m.components)

# ΣNᵢ bᵢ
covolume(m::VanDerWaalsMixture, molfrac) = dot(m.components.b, molfrac)
binary_interaction(m::VanDerWaalsMixture, molfrac) = molfrac' * m.aij * molfrac

function log_fugacity_coef(m::VanDerWaalsMixture, molfrac, pressure, RT, phase::Symbol)
    factors = zfactors(m, molfrac, pressure, RT)

    zgas = minimum(factors)
    zliq = maximum(factors)

    z = phase == :gas ? zgas : zliq

    return log_fugacity_coef(m, molfrac, pressure, RT, z)
end

function log_fugacity_coef(m::VanDerWaalsMixture, molfrac, pressure, RT, z::Real)
    b = covolume(m, molfrac)
    B = b * pressure / RT
    Aderiv = 2 * m.aij * molfrac

    return (
        .- log(z - B)
        .+ m.components.b .* pressure ./ (RT * (z - B))
        .- pressure .* Aderiv / (z * RT^2)
    )
end

function zfactors(m::VanDerWaalsMixture, molfrac, pressure, RT)
    A = binary_interaction(m, molfrac) * pressure / RT^2
    B = covolume(m, molfrac) * pressure / RT

    roots = solve_cubic(
        1,
        (- B - 1),
        A,
        - A * B
    )
    return filter(!isnan, roots)
end
