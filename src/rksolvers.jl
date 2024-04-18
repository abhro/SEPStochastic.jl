"""
Perform 1 step of the RK4 differential equation solver.

Rewritten to be more like MATLAB interface
"""
function rk4(odefun::Function, x0, y0, h; odefun_params)
    # Pretty much just taken from wikipedia
    @debug("Got called with params", odefun, x0, y0, h, odefun_params)
    k1 = odefun(x0,       y0,            odefun_params...)
    @debug("k1 = $k1")
    k2 = odefun(x0 + h/2, y0 + h/2 * k1, odefun_params...)
    @debug("k2 = $k2")
    k3 = odefun(x0 + h/2, y0 + h/2 * k2, odefun_params...)
    @debug("k3 = $k3")
    k4 = odefun(x0 + h,   y0 + h   * k3, odefun_params...)
    @debug("k4 = $k4")

    return y0 + h/6 * (k1 + 2k2 + 2k3 + k4)
end