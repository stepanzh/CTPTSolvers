# Kᵢ
michelsen(c::Component, pressure, RT) = c.critical_pressure / pressure * exp(5.42 * (1 - c.critical_thermal_energy / RT))
