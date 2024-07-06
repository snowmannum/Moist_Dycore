using Statistics
using Interpolations
export Compute_Corrections_Init, Compute_Corrections!, Four_In_One!, Spectral_Dynamics!, Get_Topography!, Spectral_Initialize_Fields!, Spectral_Dynamics_Physics!, Atmosphere_Update!



function Compute_Corrections_Init(vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    grid_u_p::Array{Float64, 3}, grid_v_p::Array{Float64, 3}, grid_ps_p::Array{Float64, 3}, grid_t_p::Array{Float64, 3}, 
    grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3}, grid_δt::Array{Float64, 3},  
    Δt::Int64, grid_energy_temp::Array{Float64, 3}, grid_tracers_p::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_δtracers::Array{Float64,3})
    
    do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
    
    sum_tracers_p = 0.

    if (do_mass_correction) 
        mean_ps_p = Area_Weighted_Global_Mean(mesh, grid_ps_p)
    end
    
    if (do_energy_correction) 
        # due to dissipation introduced by the forcing
        cp_air, grav       = atmo_data.cp_air, atmo_data.grav 
        grid_energy_temp  .= 0.5*((grid_u_p + Δt*grid_δu).^2 + (grid_v_p + Δt*grid_δv).^2) + cp_air*(grid_t_p + Δt*grid_δt)
        mean_energy_p      = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_p)
        ###
    end

    if (do_water_correction)
        # error("water correction has not implemented")
        mean_moisture_p    =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_p .+ grid_δtracers * Δt, grid_ps_p)

    end
    
    return mean_ps_p, mean_energy_p, mean_moisture_p 
end 

function Compute_Corrections!(semi_implicit::Semi_Implicit_Solver, vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    mean_ps_p::Float64, mean_energy_p::Float64, mean_moisture_p::Float64,
    grid_u_n::Array{Float64, 3}, grid_v_n::Array{Float64, 3},
    grid_energy_temp::Array{Float64, 3}, grid_ps_p::Array{Float64, 3},grid_ps_c::Array{Float64, 3},
    grid_ps_n::Array{Float64, 3}, spe_lnps_n::Array{ComplexF64, 3}, 
    grid_t_n::Array{Float64, 3}, spe_t_n::Array{ComplexF64, 3},
    grid_tracers_p::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_tracers_n::Array{Float64, 3}, 
    grid_t::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_p_half::Array{Float64, 3}, grid_z_full::Array{Float64, 3}, grid_u_p::Array{Float64, 3}, grid_v_p::Array{Float64, 3},
    grid_geopots::Array{Float64, 3}, grid_w_full::Array{Float64,3}, grid_t_p::Array{Float64, 3}, dyn_data::Dyn_Data, grid_δt::Array{Float64,3}, factor1::Array{Float64,3}, factor2::Array{Float64,3})#,  grid_tracers_diff::Array{Float64, 3})


    do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
    
    
    if (do_mass_correction) 
        mean_ps_n              = Area_Weighted_Global_Mean(mesh, grid_ps_n)
        mass_correction_factor = mean_ps_p/mean_ps_n
        grid_ps_n            .*= mass_correction_factor
        #P00 = 1 
        spe_lnps_n[1,1,1]     += log(mass_correction_factor)
    end
    
    if (do_energy_correction) 
        cp_air, grav           = atmo_data.cp_air, atmo_data.grav
        grid_energy_temp      .= 0.5*(grid_u_n.^2 + grid_v_n.^2) + cp_air*grid_t_n
        mean_energy_n          = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_n)
        
        temperature_correction = grav*(mean_energy_p - mean_energy_n)/(cp_air*mean_ps_p)
        #@info grav, mean_energy_p , mean_energy_n, cp_air, mean_ps_p
        grid_t_n             .+= temperature_correction
        spe_t_n[1,1,:]       .+= temperature_correction
    end

    # @info mean_ps_p, mean_energy_p, mass_correction_factor, temperature_correction
    ### By CJY 0517
    nλ         = mesh.nλ
    nθ         = mesh.nθ
    nd         = mesh.nd
    grav       = atmo_data.grav
    integrator = semi_implicit.integrator
    Δt         = Get_Δt(integrator)

    if (do_water_correction) 
        
        grid_tracers_n[grid_tracers_n .< 0.] .=  0.
        mean_moisture_n                       =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_n, grid_ps_n)
        grid_tracers_n                      .*=  mean_moisture_p ./ mean_moisture_n 
        mean_moisture_n                       =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_tracers_n, grid_ps_n)        
        ### 10/30 
        # @info "#### moisture correction:", (mean_moisture_n - mean_moisture_p)
        return mean_moisture_n
    end
    
end 


"""
compute vertical mass flux and velocity 
grid_M_half[:,:,k+1] = downward mass flux/per unit area across the K+1/2
grid_w_full[:,:,k]   = dp/dt vertical velocity 

update residuals
grid_δps[:,:,k]  += -∑_{r=1}^nd Dr = -∑_{r=1}^nd ∇(vrΔp_r)
grid_δt[:,:,k]  += κTw/p 
(grid_δu[:,:,k], grid_δv[:,:,k]) -= RT ∇p/p 

!  cell boundary. This is the "vertical velocity" in the hybrid coordinate system.
!  When vertical coordinate is pure sigma: grid_M_half = grid_ps*d(sigma)/dt
"""

function Four_In_One!(vert_coord::Vert_Coordinate, atmo_data::Atmo_Data, 
    grid_div::Array{Float64,3}, grid_u::Array{Float64,3}, grid_v::Array{Float64,3}, 
    grid_ps::Array{Float64,3},  grid_Δp::Array{Float64,3}, grid_lnp_half::Array{Float64,3}, grid_lnp_full::Array{Float64,3}, grid_p_full::Array{Float64,3},
    grid_dλ_ps::Array{Float64,3}, grid_dθ_ps::Array{Float64,3}, 
    grid_t::Array{Float64,3}, 
    grid_M_half::Array{Float64,3}, grid_w_full::Array{Float64,3}, 
    grid_δu::Array{Float64,3}, grid_δv::Array{Float64,3}, grid_δps::Array{Float64,3}, grid_δt::Array{Float64,3}, grid_δtracers::Array{Float64,3})
    
    rdgas, cp_air          = atmo_data.rdgas, atmo_data.cp_air
    nd, bk                 = vert_coord.nd, vert_coord.bk
    Δak, Δbk               = vert_coord.Δak, vert_coord.Δbk
    vert_difference_option = vert_coord.vert_difference_option
    
    kappa                  = rdgas / cp_air
    
    # dmean_tot = ∇ ∑_{k=1}^{nd} vk Δp_k = ∑_{k=1}^{nd} Dk
    nλ, nθ, _              = size(grid_ps)
    dmean_tot              = zeros(Float64, nλ, nθ)
    Δlnp_p                 = zeros(Float64, nλ, nθ)
    Δlnp_m                 = zeros(Float64, nλ, nθ)
    Δlnp                   = zeros(Float64, nλ, nθ)
    x1                     = zeros(Float64, nλ, nθ)
    dlnp_dλ                = zeros(Float64, nλ, nθ)
    dlnp_dθ                = zeros(Float64, nλ, nθ)
    dmean                  = zeros(Float64, nλ, nθ)
    x5                     = zeros(Float64, nλ, nθ)
        
    if (vert_difference_option == "simmons_and_burridge") 
        for k = 1:nd
        Δp       = grid_Δp[:,:,k]
        
        Δlnp_p  .= grid_lnp_half[:,:,k + 1] - grid_lnp_full[:,:,k]
        Δlnp_m  .= grid_lnp_full[:,:,k]     - grid_lnp_half[:,:,k]
        Δlnp    .= grid_lnp_half[:,:,k + 1] - grid_lnp_half[:,:,k]
        
        # angular momentum conservation 
        #    ∇p_k/p =  [(lnp_k - lnp_{k-1/2})∇p_{k-1/2} + (lnp_{k+1/2} - lnp_k)∇p_{k+1/2}]/Δpk
        #         =  [(lnp_k - lnp_{k-1/2})B_{k-1/2} + (lnp_{k+1/2} - lnp_k)B_{k+1/2}]/Δpk * ∇ps
        #         =  x1 * ∇ps
        x1      .= (bk[k] * Δlnp_m + bk[k + 1] * Δlnp_p ) ./ Δp
        
        dlnp_dλ .= x1 .* grid_dλ_ps[:,:,1]
        dlnp_dθ .= x1 .* grid_dθ_ps[:,:,1]
        
        
        # (grid_δu, grid_δv) -= RT ∇p/p 
        grid_δu[:,:,k] .-=  rdgas * grid_t[:,:,k] .* dlnp_dλ
        grid_δv[:,:,k] .-=  rdgas * grid_t[:,:,k] .* dlnp_dθ
        
        # dmean = ∇ (vk Δp_k) =  divk Δp_k + vk  Δbk[k] ∇ p_s
        dmean           .= grid_div[:,:,k] .* Δp + Δbk[k] * (grid_u[:,:,k] .* grid_dλ_ps[:,:,1] + grid_v[:,:,k] .* grid_dθ_ps[:,:,1])
        
    
        # energy conservation for temperature
        # w/p = dlnp/dt = ∂lnp/∂t + dσ ∂lnp/∂σ + v∇lnp
        # dσ ∂ξ_k/∂σ = [M_{k+1/2}(ξ_k+1/2 - ξ_k) + M_{k-1/2}(ξ_k - ξ_k-1/2)]/Δp_k
        # weight the same way (TODO)
        # vertical advection operator (M is the downward speed)
        # dσ ∂lnp_k/∂σ = [M_{k+1/2}(lnp_k+1/2 - lnp_k) + M_{k-1/2}(lnp_k - lnp_k-1/2)]/Δp_k
        # ∂lnp/∂t = 1/p ∂p/∂t = [∂p/∂t_{k+1/2}(lnp_k+1/2 - lnp_k) + ∂p/∂t_{k-1/2}(lnp_k - lnp_k-1/2)]/Δp_k
        # As we know
        # ∂p/∂t_{k+1/2} = -∑_{r=1}^k Dr - M_{k+1/2}
        
        # ∂lnp/∂t + dσ ∂lnp/∂σ =  [(-∑_{r=1}^k Dr)(lnp_k+1/2 - lnp_k) + (-∑_{r=1}^{k-1} Dr)(lnp_k - lnp_k-1/2)]/Δp_k
        #                      = -[(∑_{r=1}^{k-1} Dr)(lnp_k+1/2 - lnp_k-1/2) + D_k(lnp_k+1/2 - lnp_k)]/Δp_k
        
        x5                     .= -(dmean_tot .* Δlnp + dmean .* Δlnp_p) ./ Δp .+ grid_u[:,:,k] .* dlnp_dλ + grid_v[:,:,k] .* dlnp_dθ
        # grid_δt += κT w/p
        grid_δt[:,:,k]        .+=  kappa * grid_t[:,:,k] .* x5
        # grid_w_full = w
        grid_w_full[:,:,k]     .= x5 .* grid_p_full[:,:,k]
        # update dmean_tot to ∑_{r=1}^k ∇(vrΔp_r)
        dmean_tot             .+= dmean
        # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
        grid_M_half[:,:,k + 1] .= -dmean_tot
        end
        
    else
        error("vert_difference_option ", vert_difference_option, " is not a valid value for option")
        
    end
    # ∂ps/∂t = -∑_{r=1}^nd ∇(vrΔp_r) = -dmean_tot
    grid_δps[:,:,1]        .-= dmean_tot
    
    for k = 1:nd-1
        # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
        grid_M_half[:,:,k+1] .+= dmean_tot * bk[k+1]
    end
    
    grid_M_half[:,:,1]      .= 0.0
    grid_M_half[:,:,nd + 1] .= 0.0
end 



"""
The governing equations are
∂div/∂t = ∇ × (A, B) - ∇^2E := f^d                    
∂lnps/∂t= (-∑_k div_k Δp_k + v_k ∇ Δp_k)/ps := f^p    
∂T/∂t = -(u,v)∇T - dσ∂T∂σ + κTw/p + J:= f^t           
Φ = f^Φ                                               

implicit part: -∇^2Φ - ∇(RT∇lnp) ≈ I^d = -∇^2(γT + H2 ps_ref lnps) - ∇^2 H1 ps_ref lnps, here RT∇lnp ≈  H1 ps_ref ∇lnps
implicit part:  f^p              ≈ I^p = -ν div / ps_ref
implicit part:  - dσ∂T∂σ + κTw/p ≈ I^t = -τ div  
implicit part:  f^Φ              ≈ I^Φ = γT + H2 ps_ref lnps 

We have 
δdiv = f^d - I^d + I^d
δlnps = f^p - I^p + I^p
δT = f^t - I^t + I^t

"""

function Spectral_Dynamics!(mesh::Spectral_Spherical_Mesh,  vert_coord::Vert_Coordinate, 
    atmo_data::Atmo_Data, dyn_data::Dyn_Data, 
    semi_implicit::Semi_Implicit_Solver, L::Float64 = 0.1)
    
    # spectral equation quantities
    spe_lnps_p, spe_lnps_c, spe_lnps_n, spe_δlnps = dyn_data.spe_lnps_p, dyn_data.spe_lnps_c, dyn_data.spe_lnps_n, dyn_data.spe_δlnps
    spe_vor_p, spe_vor_c, spe_vor_n, spe_δvor     = dyn_data.spe_vor_p, dyn_data.spe_vor_c, dyn_data.spe_vor_n, dyn_data.spe_δvor
    spe_div_p, spe_div_c, spe_div_n, spe_δdiv     = dyn_data.spe_div_p, dyn_data.spe_div_c, dyn_data.spe_div_n, dyn_data.spe_δdiv
    spe_t_p, spe_t_c, spe_t_n, spe_δt             = dyn_data.spe_t_p, dyn_data.spe_t_c, dyn_data.spe_t_n, dyn_data.spe_δt
    
    # grid quantities
    grid_u_p, grid_u, grid_u_n    = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
    grid_v_p, grid_v, grid_v_n    = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
    grid_ps_p, grid_ps, grid_ps_n = dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
    grid_t_p, grid_t, grid_t_n    = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n


    # related quanties
    grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
    grid_dλ_ps, grid_dθ_ps                                 = dyn_data.grid_dλ_ps, dyn_data.grid_dθ_ps
    grid_lnps                                              = dyn_data.grid_lnps
    
    grid_div, grid_absvor, grid_vor                        = dyn_data.grid_div, dyn_data.grid_absvor, dyn_data.grid_vor
    grid_w_full, grid_M_half                               = dyn_data.grid_w_full, dyn_data.grid_M_half
    grid_geopots, grid_geopot_full, grid_geopot_half       = dyn_data.grid_geopots, dyn_data.grid_geopot_full, dyn_data.grid_geopot_half
    
    grid_energy_full, spe_energy                           = dyn_data.grid_energy_full, dyn_data.spe_energy
    
    # By CJY2
    spe_tracers_n     = dyn_data.spe_tracers_n
    spe_tracers_c     = dyn_data.spe_tracers_c
    spe_tracers_p     = dyn_data.spe_tracers_p 
    
    grid_tracers_n    = dyn_data.grid_tracers_n
    grid_tracers_c    = dyn_data.grid_tracers_c
    grid_tracers_p    = dyn_data.grid_tracers_p 
    
    grid_tracers_diff = dyn_data.grid_tracers_diff
    
    spe_δtracers      = dyn_data.spe_δtracers  
    grid_δtracers     = dyn_data.grid_δtracers 
    
    ### 11/07
    grid_z_full       = dyn_data.grid_z_full
    grid_z_half       = dyn_data.grid_z_half
    ###
    grid_w_full       = dyn_data.grid_w_full
    # todo !!!!!!!!
    #  grid_q = grid_t
    ###
    grav              = atmo_data.grav
    integrator        = semi_implicit.integrator
    Δt                = Get_Δt(integrator)
    factor1           = dyn_data.factor1 
    factor2           = dyn_data.factor2 
    factor3           = dyn_data.factor3  
    # factor4 = dyn_data.factor4  

    grid_z_full       = dyn_data.grid_z_full
    grid_z_half       = dyn_data.grid_z_half
    grid_δtracers     = dyn_data.grid_δtracers 

    K_E               = dyn_data.K_E
    ###############################################################################
    ###
    # original 
    # pressure difference
    grid_Δp             = dyn_data.grid_Δp
    # temporary variables
    grid_δQ             = dyn_data.grid_d_full1
        
    # incremental quantities
    grid_δu, grid_δv, grid_δps, grid_δlnps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δlnps, dyn_data.grid_δt
    integrator          = semi_implicit.integrator
    Δt                  = Get_Δt(integrator)


    # ###
    # V_c      = zeros(((128,64,20)))
    # za       = zeros(((128,64,20)))
    # rho      = zeros(((128,64,20)))
    # ##
    # C_E = 0.0044
    # Lv  = 2.5*10^6.
    # Rv  = atmo_data.rvgas  # 461.
    # Rd  = atmo_data.rdgas  # 287.
    # cp  = atmo_data.cp_air # 1004.
    # # # ### factor3
    # ### use n

    # # grid_δtracers .-= factor3 ./(2. .* Δt)
    # ### try
    # # @info maximum(grid_u)
    # """
    # # Cal V_c and za
    # """
    # V_c_loc, za_loc, rho_loc = Calculate_V_c_za_rho!(dyn_data, atmo_data, grid_p_half, grid_p_full, grid_ps, grid_t, grid_u, grid_v, grid_tracers_c)

    
    # V_c     .= V_c_loc
    # za      .= za_loc
    # rho     .= rho_loc
    
    # """
    # ## large-scale precipitation
    # """
    # grid_tracers_diff_new                  = HS_forcing_water_vapor!(semi_implicit, dyn_data, grid_tracers_n,  grid_t_n, grid_δt, grid_p_full, grid_u, grid_v, factor3, grid_δtracers, grid_tracers_c, grid_t, L)
    # grid_tracers_diff                     .= grid_tracers_diff_new
    # grid_tracers_c[grid_tracers_c .< 0]   .= 0     
    # grid_tracers_c .= grid_tracers_c .+ grid_δtracers .* Δt *2
    # grid_t         .= grid_t .+ grid_δt .* Δt*2
    

    # # Sensible_heat_fluxes!(mesh, atmo_data, grid_t, grid_tracers_c, grid_δt, V_c, Δt, za)

    # # Latent_heat_flux!-> == surface evaporation
    
    # Surface_evaporation!(mesh, atmo_data, grid_t, grid_tracers_c, grid_δtracers, grid_ps, V_c, za, Δt, factor1)
    
    # Implicit_PBL_Scheme!(atmo_data, grid_t, grid_t_n, grid_tracers_c, grid_tracers_n, grid_δtracers, grid_δt, grid_p_full, grid_p_half, V_c, za, Δt, factor2, K_E, rho)
    
    


    # Calculate latent heat and modify qv_current
    # HS_forcing_water_vapor!(grid_tracers_c,  grid_δtracers, grid_t, grid_δt, grid_p_full)

    mean_ps_p, mean_energy_p, mean_moisture_p = Compute_Corrections_Init(vert_coord, mesh, atmo_data,
    grid_u_p, grid_v_p, grid_ps_p, grid_t_p, 
    grid_δu, grid_δv, grid_δt,  
    Δt, grid_energy_full, grid_tracers_p, grid_tracers_c, grid_δtracers)
    
    # compute pressure based on grid_ps -> grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full 
    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)
    ###

    ###
    # compute ∇ps = ∇lnps * ps
    Compute_Gradients!(mesh, spe_lnps_c,  grid_dλ_ps, grid_dθ_ps)
    grid_dλ_ps .*= grid_ps
    grid_dθ_ps .*= grid_ps


    
    # compute grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, 
    # except the contributions from geopotential or vertical advection
    Four_In_One!(vert_coord, atmo_data, grid_div, grid_u, grid_v, grid_ps, 
    grid_Δp, grid_lnp_half, grid_lnp_full, grid_p_full,
    grid_dλ_ps, grid_dθ_ps, 
    grid_t, 
    grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, grid_δtracers)

    Compute_Geopotential!(vert_coord, atmo_data, 
    grid_lnp_half, grid_lnp_full,  
    grid_t, 
    grid_geopots, grid_geopot_full, grid_geopot_half, grid_tracers_c)

    grid_δlnps .= grid_δps ./ grid_ps
    Trans_Grid_To_Spherical!(mesh, grid_δlnps, spe_δlnps)

    # compute vertical advection, todo  finite volume method 
    Vert_Advection!(vert_coord, grid_u, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δu  .+= grid_δQ
    Vert_Advection!(vert_coord, grid_v, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δv  .+= grid_δQ
    Vert_Advection!(vert_coord, grid_t, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δt  .+= grid_δQ
    ### By CJY2 spectral tracers need to be done first 
    Vert_Advection!(vert_coord, grid_tracers_c, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme,  grid_δQ)
    grid_δtracers .+= grid_δQ 
    Add_Horizontal_Advection!(mesh, spe_tracers_c, grid_u, grid_v, grid_δtracers) 
    Trans_Grid_To_Spherical!(mesh, grid_δtracers, spe_δtracers)
    Compute_Spectral_Damping!(integrator, spe_tracers_c, spe_tracers_p, spe_δtracers)
    Filtered_Leapfrog!(integrator, spe_δtracers, spe_tracers_p, spe_tracers_c, spe_tracers_n)
    Trans_Spherical_To_Grid!(mesh, spe_tracers_n, grid_tracers_n)
    ###################################################
    Add_Horizontal_Advection!(mesh, spe_t_c, grid_u, grid_v, grid_δt)
    Trans_Grid_To_Spherical!(mesh, grid_δt, spe_δt)
   
    
    grid_absvor = dyn_data.grid_absvor
    Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)
    
    
    grid_δu .+=  grid_absvor .* grid_v
    grid_δv .-=  grid_absvor .* grid_u
    
    
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)

    grid_energy_full .= grid_geopot_full .+ 0.5 * (grid_u.^2 + grid_v.^2)
    Trans_Grid_To_Spherical!(mesh, grid_energy_full, spe_energy)
    Apply_Laplacian!(mesh, spe_energy)
    spe_δdiv .-= spe_energy
    
    
    
    Implicit_Correction!(semi_implicit, vert_coord, atmo_data,
    spe_div_c, spe_div_p, spe_lnps_c, spe_lnps_p, spe_t_c, spe_t_p, 
    spe_δdiv, spe_δlnps, spe_δt)


    
    Compute_Spectral_Damping!(integrator, spe_vor_c, spe_vor_p, spe_δvor)
    Compute_Spectral_Damping!(integrator, spe_div_c, spe_div_p, spe_δdiv)
    Compute_Spectral_Damping!(integrator, spe_t_c, spe_t_p, spe_δt)
    # ### By CJY2
    # Compute_Spectral_Damping!(integrator, spe_tracers_c, spe_tracers_p, spe_δtracers)
    # ###

        
    Filtered_Leapfrog!(integrator, spe_δvor, spe_vor_p, spe_vor_c, spe_vor_n)
    Filtered_Leapfrog!(integrator, spe_δdiv, spe_div_p, spe_div_c, spe_div_n)
    Filtered_Leapfrog!(integrator, spe_δlnps, spe_lnps_p, spe_lnps_c, spe_lnps_n)
    Filtered_Leapfrog!(integrator, spe_δt, spe_t_p, spe_t_c, spe_t_n)
    
    
    Trans_Spherical_To_Grid!(mesh, spe_vor_n, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_n, grid_div)
    UV_Grid_From_Vor_Div!(mesh, spe_vor_n, spe_div_n, grid_u_n, grid_v_n)
    Trans_Spherical_To_Grid!(mesh, spe_lnps_n, grid_lnps)
    grid_ps_n .= exp.(grid_lnps)
    Trans_Spherical_To_Grid!(mesh, spe_t_n, grid_t_n) 


    
    mean_moisture_n_loc = Compute_Corrections!(semi_implicit, vert_coord, mesh, atmo_data, mean_ps_p, mean_energy_p,mean_moisture_p, 
        grid_u_n, grid_v_n,
        grid_energy_full, grid_ps_p,grid_ps,
        grid_ps_n, spe_lnps_n, 
        grid_t_n, spe_t_n, 
        grid_tracers_p, grid_tracers_c, grid_tracers_n,
        grid_t, grid_p_full, grid_p_half, grid_z_full, grid_u_p, grid_v_p, grid_geopots, grid_w_full, grid_t_p, dyn_data, grid_δt, factor1, factor2)

   
    day_to_sec = 86400
    if (integrator.time%(day_to_sec/4) == 0)
        # dyn_data.grid_tracers_c[dyn_data.grid_tracers_c .< 0] .= 0
        @info "Day: ", (integrator.time/ day_to_sec), " Max |U|,|V|,|P|,|T|,|qv|: ", maximum(abs.(dyn_data.grid_u_c)), maximum(abs.(dyn_data.grid_v_c)), maximum(dyn_data.grid_p_full), maximum(dyn_data.grid_t_c), maximum(dyn_data.grid_tracers_c), maximum(dyn_data.grid_tracers_diff)
        @info "Day: ", (integrator.time/ day_to_sec), " Min |U|,|V|,|P|,|T|,|qv|: ", minimum(abs.(dyn_data.grid_u_c)), minimum(abs.(dyn_data.grid_v_c)), minimum(dyn_data.grid_p_full), minimum(dyn_data.grid_t_c), minimum(dyn_data.grid_tracers_c)
    end

    Time_Advance!(dyn_data)

    # @info "min dyn.grid_tracers_c" minimum(dyn_data.grid_tracers_n)
    
    #@info "sec: ", integrator.time+1200, sum(abs.(grid_u_n)), sum(abs.(grid_v_n)), sum(abs.(grid_t_n)) , sum(abs.(grid_ps_n))
    #@info "max: ", maximum(abs.(grid_u_n)), maximum(abs.(grid_v_n)), maximum(abs.(grid_t_n)) , maximum(abs.(grid_ps_n))
    #@info "loc", grid_u_n[100,30,10],  grid_t_n[100,30,10], grid_u_n[1,32,1],  grid_t_n[1,32,1]
    
    #@assert(maximum(grid_u) <= 100.0 && maximum(grid_v) <= 100.0)


    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)
    
    return 
end 

function Get_Topography!(grid_geopots::Array{Float64, 3}, warm_start_file_name::String = "None", initial_day::Int64 = 5)
    # original_start = true
    # load_old_file  = false
    if warm_start_file_name == "None" # load warm start file
        grid_geopots .= 0.0
    end
    ### 2023/10/25
    # read_file     = load("/work/kaichiht/Colab/2023_research/annular_mode/300day_dry_run_all.dat")
    # grid_geopots .= read_file["grid_geopots_xyzt"][:,:,1,300]
    if warm_start_file_name != "None" # load warm start file
        read_file     = load(warm_start_file_name)
        grid_geopots .= read_file["grid_geopots_final"][:,:,:,1]
    end
    
    return
end 

function Spectral_Initialize_Fields!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, vert_coord::Vert_Coordinate, sea_level_ps_ref::Float64, init_t::Float64, grid_geopots::Array{Float64,3}, dyn_data::Dyn_Data, Δt::Int64, warm_start_file_name::String = "None", initial_day::Int64 = 5)
    # load_old_file  = false
    # original_start = true
    if warm_start_file_name != "None" # load warm start file
        spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c = dyn_data.spe_vor_c, dyn_data.spe_div_c, dyn_data.spe_lnps_c, dyn_data.spe_t_c
        spe_vor_p, spe_div_p, spe_lnps_p, spe_t_p = dyn_data.spe_vor_p, dyn_data.spe_div_p, dyn_data.spe_lnps_p, dyn_data.spe_t_p
        grid_u, grid_v, grid_ps, grid_t           = dyn_data.grid_u_c, dyn_data.grid_v_c, dyn_data.grid_ps_c, dyn_data.grid_t_c
        grid_u_p, grid_v_p, grid_ps_p, grid_t_p   = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_ps_p, dyn_data.grid_t_p
        
        grid_lnps,  grid_vor, grid_div            = dyn_data.grid_lnps, dyn_data.grid_vor, dyn_data.grid_div
        
        grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_Δp, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
        nλ, nθ, nd                                = mesh.nλ, mesh.nθ, mesh.nd
        
        ### By CJY2
        grid_t_n          = dyn_data.grid_t_n
        spe_tracers_c     = dyn_data.spe_tracers_c
        spe_tracers_p     = dyn_data.spe_tracers_p 

        grid_tracers_n    = dyn_data.grid_tracers_n
        grid_tracers_c    = dyn_data.grid_tracers_c
        grid_tracers_p    = dyn_data.grid_tracers_p 

        grid_u_n      = dyn_data.grid_u_n
        grid_v_n      = dyn_data.grid_v_n
        ########################################################
        # Tendency 
        grid_δu = dyn_data.grid_δu
        grid_δv = dyn_data.grid_δv

        grid_δtracers = dyn_data.grid_δtracers
        ########################################################
        @info warm_start_file_name # to make sure get the correct warmstart_PR.dat
        read_file     = load(warm_start_file_name)        
        grid_u[:,:,:]    .= read_file["grid_u_c_final"][:,:,:,1]
        grid_v[:,:,:]    .= read_file["grid_v_c_final"][:,:,:,1]  
        grid_t       .= read_file["grid_t_c_final"][:,:,:,1] 
        
        grid_lnps    .= log.(read_file["grid_ps_c_final"][:,:,1,1])
        grid_ps      .= read_file["grid_ps_c_final"][:,:,1,1]

        # grid_ps -> grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full
        Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp,
        grid_lnp_half, grid_p_full, grid_lnp_full)
        ########################################################
        # By CJY
        num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical

        spe_t_c       .= read_file["spe_t_c_final"][:,:,:,1] 
        
        
        spe_vor_c[:,:,:] .= read_file["spe_vor_c_final"][:,:,:,1]
        spe_div_c[:,:,:] .= read_file["spe_div_c_final"][:,:,:,1]
        
        spe_lnps_c    .= (read_file["spe_lnps_c_final"][:,:,1,1])
        

        grid_vor .= read_file["grid_vor_final"][:,:,:,1] # Compute_Abs_Vor! need it 
        grid_div .= read_file["grid_div_final"][:,:,:,1] # Four_in_one! need it
        ########################################################
        spe_vor_p   .= read_file["spe_vor_p_final"][:,:,:,1]
        spe_div_p   .= read_file["spe_div_p_final"][:,:,:,1]
        spe_lnps_p  .= read_file["spe_lnps_p_final"][:,:,:,1]
        spe_t_p     .= read_file["spe_t_p_final"][:,:,:,1]

        grid_u_p    .= read_file["grid_u_p_final"][:,:,:,1]
        grid_v_p    .= read_file["grid_v_p_final"][:,:,:,1]
        grid_ps_p   .= read_file["grid_ps_p_final"][:,:,:,1]
        grid_t_p    .= read_file["grid_t_p_final"][:,:,:,1]
        ########################################################
        # Tracer initialization
        grid_tracers_n .= read_file["grid_tracers_n_final"][:,:,:,1] # large precipitation need next DO NOT REMOVE IT !!!
        grid_tracers_c .= read_file["grid_tracers_c_final"][:,:,:,1]
        grid_tracers_p .= read_file["grid_tracers_p_final"][:,:,:,1]
        
        # Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        # Trans_Grid_To_Spherical!(mesh, grid_tracers_p, spe_tracers_p)
        spe_tracers_c  .= read_file["spe_tracers_c_final"][:,:,:,1]
        spe_tracers_p  .= read_file["spe_tracers_p_final"][:,:,:,1]

        # grid_δu .= read_file["grid_δu_xyzt"][:,:,:,initial_day] # Rayleigh_Damping! has given it value 
        # grid_δv .= read_file["grid_δv_xyzt"][:,:,:,initial_day] # Rayleigh_Damping! has given it value

        ####################################################################
        # Correction_Init! would use these!!!
        # grid_t_n    .= read_file["grid_t_n_xyzt"][:,:,:,initial_day] 
        # grid_u_n   .= read_file["grid_u_n_xyzt"][:,:,:,initial_day]
        # grid_v_n   .= read_file["grid_v_n_xyzt"][:,:,:,initial_day]
        # grid_δtracers .= read_file["grid_δtracers_xyzt"][:,:,:,initial_day]

    end

    if warm_start_file_name == "None" # then use original start
        spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c = dyn_data.spe_vor_c, dyn_data.spe_div_c, dyn_data.spe_lnps_c, dyn_data.spe_t_c
        spe_vor_p, spe_div_p, spe_lnps_p, spe_t_p = dyn_data.spe_vor_p, dyn_data.spe_div_p, dyn_data.spe_lnps_p, dyn_data.spe_t_p
        grid_u, grid_v, grid_ps, grid_t = dyn_data.grid_u_c, dyn_data.grid_v_c, dyn_data.grid_ps_c, dyn_data.grid_t_c
        grid_u_p, grid_v_p, grid_ps_p, grid_t_p = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_ps_p, dyn_data.grid_t_p
        
        grid_lnps,  grid_vor, grid_div =  dyn_data.grid_lnps, dyn_data.grid_vor, dyn_data.grid_div
        
        grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_Δp, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
        nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd
                
        ### By CJY2
        spe_tracers_c     = dyn_data.spe_tracers_c
        spe_tracers_p     = dyn_data.spe_tracers_p 
            
        grid_tracers_c    = dyn_data.grid_tracers_c
        grid_tracers_p    = dyn_data.grid_tracers_p 

        rdgas = atmo_data.rdgas
        grid_t         .=  init_t 

        # 𝚪 = 0.005
        # a = 6.371E6
        # b = 2
        # k = 3
        # p0 = 100000
        # Rd = 287
        # g  = 9.81
        # T0P = 240
        # T0E = 310
        # T0 = (T0E + T0P) * 0.5
        # H = Rd * T0/g
        
        # grid_z_full = zeros(((128,64,20)))
        # dry_run_file = load("test_final.dat")
        # grid_z_full .= dry_run_file["grid_z_full_xyzt"][:,:,:,5]
        
        # A = 1/𝚪 
        # B = (T0E - T0P) / (T0E + T0P) / T0P
        # C = (k+2)/2 * (T0E - T0P) / (T0E * T0P)
        
        # τ1 = zeros(((128,64,20)))
        # τ2 = zeros(((128,64,20)))
        
        # τ1 .= A * 𝚪 / T0 .* exp.(𝚪/T0 .* grid_z_full) .+ B .* (1 .- 2 .* (grid_z_full./(b*H)).^2) .* exp.(-1 .* (grid_z_full./(b*H)).^2) 
        
        # τ2 .= C .* (1 .- 2 .* (grid_z_full./(b*H)).^2) .* exp.(-1 .* (grid_z_full./(b*H)).^2) 
        
        # θc2  = LinRange(-90,90,64)
        # θc  = deg2rad.(θc2)
        # for j in 1:64
        #     grid_t[:,j,:] .= (τ1[:,j,:] .- τ2[:,j,:] .* ((cos(θc[j]))^k - (k/(k+2)) .* (cos(θc[j]))^(k+2))).^-1
        # end

        
        # dΦ/dlnp = -RT    Δp = -ΔΦ/RT
        grid_lnps[:,:,1] .= log(sea_level_ps_ref) .- grid_geopots[:,:,1] ./ (rdgas * init_t) 
        grid_ps   .= exp.(grid_lnps)
        
        
        spe_div_c .= 0.0
        spe_vor_c .= 0.0
      
        # # initial perturbation
        num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical
        
        initial_perturbation = 1.0e-7/sqrt(2.0)
        # initial vorticity perturbation used in benchmark code
        # In gfdl spe[i,j] =  myspe[i, i+j-1]*√2
        for k = nd-2:nd
          spe_vor_c[2,5,k] = initial_perturbation
          spe_vor_c[6,9,k] = initial_perturbation
          spe_vor_c[2,4,k] = initial_perturbation  
          spe_vor_c[6,8,k] = initial_perturbation
        end
      
        UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
        
      
        # initial spectral fields (and spectrally-filtered) grid fields
        
        Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)
      
        Trans_Grid_To_Spherical!(mesh, grid_lnps, spe_lnps_c)
        Trans_Spherical_To_Grid!(mesh, spe_lnps_c,  grid_lnps)
        grid_ps .= exp.(grid_lnps)
        
      
        Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c)
      
        UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
        
        Trans_Spherical_To_Grid!(mesh, spe_vor_c, grid_vor)
        Trans_Spherical_To_Grid!(mesh, spe_div_c, grid_div)
        
        #update pressure variables for hs forcing
        Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp,
        grid_lnp_half, grid_p_full, grid_lnp_full)
        
        
        spe_vor_p  .= spe_vor_c
        spe_div_p  .= spe_div_c
        spe_lnps_p .= spe_lnps_c
        spe_t_p    .= spe_t_c
      
      
        grid_u_p   .= grid_u
        grid_v_p   .= grid_v
        grid_ps_p  .= grid_ps
        grid_t_p   .= grid_t

        # Tracer initialization
        initial_RH      = 0.8
        Lv              = 2.5*10^6.
        Rv              = 461.
        qv0             = 0.018
        θc              = mesh.θc # lat
        phi_hw          = 2 * pi / 9 * deg2rad(40)
        p_hw            = 30000.
        phi             = LinRange(-90,90,64)
        p0              = 100000.
        for k in 1:20
            for j in 1:64
               for i in 1:128
                   # grid_tracers_c[i,j,k] = qv0 * exp(-((grid_p_full[i,j,k]/grid_ps[i,j,1] - 1.)*(p0/p_hw))^2) * exp(-((deg2rad(phi[j]))/phi_hw)^4) 
                   grid_tracers_c[i,j,k] = qv0 * exp(-((grid_p_full[i,j,k]/grid_ps[i,j,1] - 1.)*(p0/p_hw))^2) * exp(-((θc[j])/phi_hw)^4) 
                    
               end            
            end
        end
        grid_tracers_c[:,:,1] .= 0.
        
        Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)

        grid_tracers_p .= grid_tracers_c
        spe_tracers_p  .= spe_tracers_c
    end
      

     
end 


function Spectral_Dynamics_Physics!(semi_implicit::Semi_Implicit_Solver, atmo_data::Atmo_Data, mesh::Spectral_Spherical_Mesh, dyn_data::Dyn_Data, Δt::Int64, physics_params::Dict{String, Float64}, L::Float64)
    grid_δu, grid_δv, grid_δps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δt
    grid_u_p, grid_v_p,  grid_t_p       = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_t_p
    grid_p_half, grid_p_full            = dyn_data.grid_p_half, dyn_data.grid_p_full
    grid_t_eq                           = dyn_data.grid_t_eq

    grid_δtracers                       = dyn_data.grid_δtracers
    spe_δtracers                        = dyn_data.spe_δtracers

    
    grid_δps .= 0.0

    spe_δtracers   .= 0.
    grid_δtracers  .= 0.

    #####################################################################################################
     # spectral equation quantities
    spe_lnps_p, spe_lnps_c, spe_lnps_n, spe_δlnps = dyn_data.spe_lnps_p, dyn_data.spe_lnps_c, dyn_data.spe_lnps_n, dyn_data.spe_δlnps
    spe_vor_p, spe_vor_c, spe_vor_n, spe_δvor     = dyn_data.spe_vor_p, dyn_data.spe_vor_c, dyn_data.spe_vor_n, dyn_data.spe_δvor
    spe_div_p, spe_div_c, spe_div_n, spe_δdiv     = dyn_data.spe_div_p, dyn_data.spe_div_c, dyn_data.spe_div_n, dyn_data.spe_δdiv
    spe_t_p, spe_t_c, spe_t_n, spe_δt             = dyn_data.spe_t_p, dyn_data.spe_t_c, dyn_data.spe_t_n, dyn_data.spe_δt
    
    # grid quantities
    grid_u_p, grid_u, grid_u_n    = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
    grid_v_p, grid_v, grid_v_n    = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
    grid_ps_p, grid_ps, grid_ps_n = dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
    grid_t_p, grid_t, grid_t_n    = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n


    # related quanties
    grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
    grid_dλ_ps, grid_dθ_ps                                 = dyn_data.grid_dλ_ps, dyn_data.grid_dθ_ps
    grid_lnps                                              = dyn_data.grid_lnps
    
    grid_div, grid_absvor, grid_vor                        = dyn_data.grid_div, dyn_data.grid_absvor, dyn_data.grid_vor
    grid_w_full, grid_M_half                               = dyn_data.grid_w_full, dyn_data.grid_M_half
    grid_geopots, grid_geopot_full, grid_geopot_half       = dyn_data.grid_geopots, dyn_data.grid_geopot_full, dyn_data.grid_geopot_half
    
    grid_energy_full, spe_energy                           = dyn_data.grid_energy_full, dyn_data.spe_energy
    
    # By CJY2
    spe_tracers_n     = dyn_data.spe_tracers_n
    spe_tracers_c     = dyn_data.spe_tracers_c
    spe_tracers_p     = dyn_data.spe_tracers_p 
    
    grid_tracers_n    = dyn_data.grid_tracers_n
    grid_tracers_c    = dyn_data.grid_tracers_c
    grid_tracers_p    = dyn_data.grid_tracers_p 
    
    grid_tracers_diff = dyn_data.grid_tracers_diff
    
    spe_δtracers      = dyn_data.spe_δtracers  
    grid_δtracers     = dyn_data.grid_δtracers 
    
    ### 11/07
    grid_z_full       = dyn_data.grid_z_full
    grid_z_half       = dyn_data.grid_z_half
    ###
    grid_w_full       = dyn_data.grid_w_full
    # todo !!!!!!!!
    #  grid_q = grid_t
    ###
    grav              = atmo_data.grav
    integrator        = semi_implicit.integrator
    Δt                = Get_Δt(integrator)
    factor1           = dyn_data.factor1 
    factor2           = dyn_data.factor2 
    factor3           = dyn_data.factor3  
    # factor4 = dyn_data.factor4  

    grid_z_full       = dyn_data.grid_z_full
    grid_z_half       = dyn_data.grid_z_half
    grid_δtracers     = dyn_data.grid_δtracers 

    K_E               = dyn_data.K_E
    ###############################################################################
    ###
    # original 
    # pressure difference
    grid_Δp             = dyn_data.grid_Δp
    # temporary variables
    grid_δQ             = dyn_data.grid_d_full1
        
    # incremental quantities
    grid_δu, grid_δv, grid_δps, grid_δlnps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δlnps, dyn_data.grid_δt
    integrator          = semi_implicit.integrator
    Δt                  = Get_Δt(integrator)

    spe_tracers_c = dyn_data.spe_tracers_c
    spe_t_c = dyn_data.spe_t_c
    

    ###
    V_c      = zeros(((128,64,20)))
    za       = zeros(((128,64,20)))
    rho      = zeros(((128,64,20)))
    ##
    C_E = 0.0044
    Lv  = 2.5*10^6.
    Rv  = atmo_data.rvgas  # 461.
    Rd  = atmo_data.rdgas  # 287.
    cp  = atmo_data.cp_air # 1004.
    # # ### factor3
    ### use n

    # grid_δtracers .-= factor3 ./(2. .* Δt)
    ### try
    # @info maximum(grid_u)
    """
    # Cal V_c and za
    """
    V_c_loc, za_loc, rho_loc = Calculate_V_c_za_rho!(dyn_data, atmo_data, grid_p_half, grid_p_full, grid_ps, grid_t, grid_u, grid_v, grid_tracers_c)

    V_c     .= V_c_loc
    za      .= za_loc
    rho     .= rho_loc
    
    """
    ## large-scale precipitation
    """
    do_large_scale_precipitation = true
    do_Sensible_heat_fluxes      = true
    do_Surface_evaporation       = true
    do_Implicit_PBL_Scheme       = true
    
    if do_large_scale_precipitation == true
        HS_forcing_water_vapor!(semi_implicit, dyn_data, grid_tracers_n,  grid_t_n, grid_δt, grid_p_full, grid_u, grid_v, grid_δtracers, grid_tracers_c, grid_t, grid_tracers_diff, factor3, L)
        grid_tracers_c[grid_tracers_c .< 0]   .= 0     
    
        
        grid_tracers_c .= grid_tracers_c .- grid_δtracers .* (2*Δt)
        grid_t         .= grid_t         .+ grid_δt       .* (2*Δt)
    
        Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)
        
        Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)
        
        grid_δtracers .= 0.
        grid_δt       .= 0.
    end

    # Calculate grid_δt(.+=) and grid_t(.=)
    if do_Sensible_heat_fluxes == true
        Sensible_heat_fluxes!(mesh, atmo_data, grid_t, grid_t_n, grid_tracers_c, grid_δt, V_c, Δt, za)
        Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)
    end

    if do_Surface_evaporation == true
        # Calculate grid_δtracers(.+=) and grid_tracers_c(.=)  (Latent_heat_flux! == Surface_evaporation!)
        Surface_evaporation!(mesh, atmo_data, grid_t, grid_tracers_c, grid_tracers_n, grid_δtracers, grid_ps, V_c, za, Δt, factor1)
        Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)
    end

    # Calculate {grid_δtracers(.+=) and grid_tracers_c(.=)} and {grid_δt(.+=) and grid_t(.=)}
    if do_Implicit_PBL_Scheme == true
        Implicit_PBL_Scheme!(atmo_data, grid_t, grid_t_n, grid_tracers_c, grid_tracers_n, grid_δtracers, grid_δt, grid_p_full, grid_p_half, V_c, za, Δt, factor2, K_E, rho)
    
        Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)
        
        Trans_Grid_To_Spherical!(mesh, grid_tracers_c, spe_tracers_c)
        Trans_Spherical_To_Grid!(mesh, spe_tracers_c, grid_tracers_c)
    end
    
    
    
    ######################################################################################################
    HS_Forcing!(atmo_data, Δt, mesh.sinθ, grid_u_p, grid_v_p, grid_p_half, grid_p_full, grid_t_p, grid_δu, grid_δv,
    grid_t_eq, grid_δt, physics_params)

    
end


function Atmosphere_Update!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, vert_coord::Vert_Coordinate, semi_implicit::Semi_Implicit_Solver, 
                            dyn_data::Dyn_Data, physcis_params::Dict{String, Float64}, L::Float64)

    Δt = Get_Δt(semi_implicit.integrator)
    Spectral_Dynamics_Physics!(semi_implicit, atmo_data, mesh,  dyn_data, Δt, physcis_params, L) # HS forcing
    Spectral_Dynamics!(mesh,  vert_coord , atmo_data, dyn_data, semi_implicit, L) # dynamics 
    ### original
    grid_ps , grid_Δp, grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_ps_c,  dyn_data.grid_Δp, dyn_data.grid_p_half, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full 
    
    grid_t = dyn_data.grid_t_c
    grid_geopots, grid_z_full, grid_z_half = dyn_data.grid_geopots, dyn_data.grid_z_full, dyn_data.grid_z_half

    ### 1201
    grid_tracers_c = dyn_data.grid_tracers_c
        
    Compute_Pressures_And_Heights!(atmo_data, vert_coord,     
    grid_ps, grid_geopots, grid_t, 
    grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full, grid_z_full, grid_z_half, grid_tracers_c)

    return
end 


function HS_forcing_water_vapor!(semi_implicit::Semi_Implicit_Solver, dyn_data::Dyn_Data, grid_tracers_n::Array{Float64, 3},  grid_t_n::Array{Float64, 3}, grid_δt::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_u::Array{Float64, 3},  grid_v::Array{Float64, 3}, grid_δtracers::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_t::Array{Float64, 3}, grid_tracers_diff::Array{Float64, 3}, factor3::Array{Float64, 3}, L::Float64)

    integrator = semi_implicit.integrator
    Δt         = Get_Δt(integrator)
    cp         = 1004.
    Lv         = 2.5*10^6.
    Rd         = 287.04
    Rv         = 461.
    # C_E = 0.0044
    # grid_tracers_diff      = zeros(size(grid_tracers_c)...)
    grid_tracers_c_max     = zeros(size(grid_tracers_c)...)
    
    #grid_tracers_c_max     = deepcopy(grid_tracers_c)
    grid_tracers_c_max    .= (0.622 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t)) )) ./ (grid_p_full .- 0.378 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t)) )) 

    dq_sat_dT              = zeros(size(grid_tracers_c)...)
    dq_sat_dT             .= Lv.*grid_tracers_c_max./ (Rv .*grid_t.^2)
    #@info "max: ", maximum(dq_sat_dT)
    # grid_d_full2 = dyn_data.grid_d_full2
    # @info maximum(grid_d_full2), minimum(grid_d_full2)

    ### Condensation_rate == grid_tracers_diff
    grid_tracers_diff     .= (max.(grid_tracers_c, grid_tracers_c_max) .- grid_tracers_c_max) ./ (1 .+ (Lv / cp) .* dq_sat_dT) ./(2 .* Δt)
    # grid_tracers_c       .-= (max.(grid_tracers_c, grid_tracers_c_max) .- grid_tracers_c_max) ./ (1 .+ (Lv / cp) .* dq_sat_dT)
    ### error ###
    grid_δtracers       .= (max.(grid_tracers_c, grid_tracers_c_max) .- grid_tracers_c_max) ./ (1 .+ (Lv / cp) .* dq_sat_dT) /(2 .* Δt)
    #############

    # latent heat feedback to temperature tendency 
    # day_to_sec        = 86400.
    # L                 = 0.05
    # @info "L:", L
    factor3          .= grid_tracers_diff
    # diabatic_heating  = deepcopy(grid_tracers_diff)
    # diabatic_heating .= (grid_tracers_diff .* Lv ./ cp) ./day_to_sec .* L 
    # @info "max: ", maximum(diabatic_heating)
    
    grid_δt         .= (grid_tracers_diff .* Lv ./ cp) .* L 
    # @info "L=", L
    
    # grid_t         .+= (grid_tracers_diff .* Lv ./ cp) .* L 

    
    ###
end

function Calculate_V_c_za_rho!(dyn_data::Dyn_Data, atmo_data::Atmo_Data, grid_p_half::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_ps::Array{Float64, 3}, grid_t::Array{Float64, 3}, grid_u::Array{Float64, 3}, grid_v::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3})
    ##
    C_E = 0.0044
    Lv  = 2.5*10^6.
    Rv  = atmo_data.rvgas  # 461.
    Rd  = atmo_data.rdgas  # 287.
    cp  = atmo_data.cp_air # 1004.
    grav = atmo_data.grav
    # # ### factor3
    ### use n

    # grid_δtracers .-= factor3 ./(2. .* Δt)
    ### try
    # @info maximum(grid_u)
    """
    # Cal V_c and za
    """
    # V_c  = zeros(((128,64,20)))
    # V_c .= (grid_u[:,:,:].^2 .+ grid_v[:,:,:].^2).^0.5
    
    # Calculate V_c
    V_c = sqrt.(grid_u .^ 2 .+ grid_v .^ 2)
    
    ### add moisture at surface following paper
    ### ∂q_a/∂t = C_E * V_a * (q_sat,a - q_a) ./ z_a 

    # rho_s = zeros(((128,64,1)))
    # rho_s[:,:,1] .=  grid_ps_n[:,:,1] ./ Rd ./ (grid_t[:,:,20])

    # cal za
    tv         = zeros(((128,64,1)))
    za         = zeros(((128,64,1)))
    tv[:,:,1] .= grid_t[:,:,20] .* (1. .+ 0.608 .* grid_tracers_c[:,:,20])
    za[:,:,1] .= Rd .* tv[:,:,1] ./grav .* (log.(grid_ps[:,:,1] ./ ((grid_p_full[:,:,20] .+ grid_p_half[:,:,21]) ./ 2.) )) ./2
    # za[:,:,1] .= Rd .* tv[:,:,1] ./ grav .* (log.(grid_ps[:,:,1]) .- log.(grid_p_half[:,:,20])) .* 0.5

    # cal rho
    rho = zeros(((128,64,20)))
    for i in 1:20
        rho[:,:,i] .=  grid_p_full[:,:,i] ./ Rd ./ (grid_t[:,:,i].* (1. .+ 0.608 .* grid_tracers_c[:,:,20]))
    end
    
    # @info "#### za global minimum, maximum:" minimum(za), maximum(za)
    return V_c, za, rho
end

function Sensible_heat_fluxes!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, grid_t::Array{Float64, 3}, grid_t_n::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_δt::Array{Float64, 3}, V_c::Array{Float64, 3}, Δt::Int64, za::Array{Float64, 3})

    C_E = 0.0044
    θc  = mesh.θc
    Tsurf = zeros((128,64))
    Tsurf = deepcopy(grid_t[:,:,20]) .*0
    for i in 1:64
         Tsurf[:,i] .= 29. .* exp.(-(θc[i] .^2. ./ (2 * (26. * pi / 180.)^2.))) .+ 271.
    end

   # grid_δt[:,:,20] .+= (((C_E .* V_c[:,:,20] .* (Tsurf[:,:,1] .- min.(grid_t[:,:,20], Tsurf[:,:,1])) .* Δt ./ za[:,:,1])
                         # ./ (1. .+ C_E .* V_c[:,:,20] .* Δt ./ za[:,:,1])) ./ (2. * Δt))
   # grid_t[:,:,20]  .= ((grid_t[:,:,20] .+ C_E .* V_c[:,:,20] .* max.(grid_t[:,:,20],Tsurf[:,:,1]) .* Δt ./ za[:,:,1]) 
                         # ./ (1. .+ C_E .* V_c[:,:,20] .* Δt ./ za[:,:,1]))
#    @info "max: ", maximum(grid_δt[:,:,20])
#    @info "min: ", minimum(grid_δt[:,:,20])
#    max.(grid_t[:,:,20],Tsurf[:,:,1])
    
    # grid_δt[:,:,20] .+= (((grid_t[:,:,20] .+ C_E .* V_c[:,:,20] .* Tsurf[:,:,1] .* Δt ./ za[:,:,1])
                        # ./ (1. .+ C_E .* V_c[:,:,20] .* Δt ./ za[:,:,1]) .- grid_t[:,:,20])  ./ (2. * Δt))
    grid_t[:,:,20]  .= ((grid_t[:,:,20] .+ C_E .* V_c[:,:,20] .* Tsurf[:,:,1] .* Δt ./ za[:,:,1]) 
                        ./ (1. .+ C_E .* V_c[:,:,20] .* Δt ./ za[:,:,1]))

end

function Surface_evaporation!(mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, grid_t::Array{Float64, 3}, grid_tracers_c::Array{Float64,3}, grid_tracers_n::Array{Float64,3}, grid_δtracers::Array{Float64,3}, grid_ps::Array{Float64, 3}, V_c::Array{Float64, 3}, za::Array{Float64, 3}, Δt::Int64, factor1::Array{Float64, 3})

    C_E = 0.0044
    Lv  = 2.5*10^6.
    Rv  = atmo_data.rvgas  # 461.
    Rd  = atmo_data.rdgas  # 287.
    
    θc = mesh.θc
    Tsurf = zeros((128,64))
    Tsurf = deepcopy(grid_t[:,:,20]) .*0
    for i in 1:64
        Tsurf[:,i] .= 29. .* exp.(-(θc[i] .^2. ./ (2 * (26. * pi / 180.)^2.))) .+ 271.
    end
    surface_evaporation         = deepcopy(grid_δtracers).*0
    
    grid_tracers_c_ps_max           = zeros(((128,64,1))) 
    grid_tracers_c_ps_max          .= (0.622 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ Tsurf[:,:])) )) ./ (grid_ps[:,:,1] .- 0.378 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ Tsurf[:,:])))) 
    # 
    
    ###########################################################
    surface_evaporation[:,:,20] .= ((C_E .* V_c[:,:,20] .* Δt ./ za[:,:,1] .*  (grid_tracers_c_ps_max[:,:,1] .- min.(grid_tracers_c[:,:,20], grid_tracers_c_ps_max[:,:,1]))) ./ (1. .+ C_E .* V_c[:,:,20] .* Δt ./ za[:,:,1])) 
    
    # grid_δtracers[:,:,20]     .+= surface_evaporation[:,:,20] ./(2. .* Δt) 
    
    grid_tracers_c[:,:,20]      .= ((grid_tracers_c[:,:,20] .+ C_E .* V_c[:,:,20] .* max.(grid_tracers_c[:,:,20],grid_tracers_c_ps_max[:,:,1]) .* Δt ./ za[:,:,1]) ./ (1. .+ C_E .* V_c[:,:,20]  .* Δt ./ za[:,:,1]))
    #########################################################

    ### try original code ###
    # grid_δtracers[:,:,20] .= (((grid_tracers_c[:,:,20] .+ C_E .* V_c[:,:,20] .* grid_tracers_c_ps_max[:,:,1] .* Δt ./ za[:,:,1]) 
                                    # ./ (1. .+ C_E .* V_c[:,:,20] .* Δt ./ za[:,:,1]) .- grid_tracers_c[:,:,20])./ 2. .* Δt)

    # grid_tracers_n[:,:,20]      .= ((grid_tracers_c[:,:,20] .+ C_E .* V_c[:,:,20] .* grid_tracers_c_ps_max[:,:,1] .* Δt ./ za[:,:,1]) 
                                    # ./ (1. .+ C_E .* V_c[:,:,20]  .* Δt ./ za[:,:,1]))
    ##########################################################
    factor1[:,:,20]              .= grid_tracers_c[:,:,20] ./(2. .* Δt) 
    

    
end 

function Implicit_PBL_Scheme!(atmo_data::Atmo_Data,grid_t::Array{Float64, 3}, grid_t_n::Array{Float64, 3}, grid_tracers_c::Array{Float64, 3}, grid_tracers_n::Array{Float64, 3}, grid_δtracers::Array{Float64, 3},  grid_δt::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_p_half::Array{Float64, 3}, V_c::Array{Float64, 3}, za::Array{Float64, 3}, Δt::Int64, factor2::Array{Float64, 3}, K_E::Array{Float64, 3}, rho::Array{Float64, 3})

    """
    Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
    
    Step1. Calculate K_E = C_E * V_c * za
    Step2. Calculate A, B, C --> E, F
    Step3. Calculate new grid_tracers_c and grid_t
    """
    C_E = 0.0044
    Lv  = 2.5*10^6.
    Rv  = atmo_data.rvgas  # 461.
    Rd  = atmo_data.rdgas  # 287.    
    cp  = atmo_data.cp_air # 1004.
    
    V_a = V_c[:,:,20]
    
    grav = atmo_data.grav
    
    ### 12/28 upgrade output_manager
    # oringin is K_E = dyn_data.K_E
    # K_E = zeros(((128,64,20+1)))
    ###
    for i in 17:21
        K_E[:,:,i] .= C_E .* V_a .* za[:,:,1]
    end
    K_E[:,:, 1:16] .= C_E .* V_a .* za[:,:,1] .* exp.(-((85000. .- grid_p_half[:,:,1:16]) ./ 10000.).^2)
    ### cal PBL Scheme
    rpdel  = zeros(((128,64,20))) ### = 1 / (p^n_{+} - p^n_{-}) , which p^_{-} mean upper layer
    for i in 1:20
        rpdel[:,:,i] .= 1. ./ (grid_p_half[:,:,i+1] .- grid_p_half[:,:,i])
    end

    CA     = zeros(((128,64,20)))
    CC     = zeros(((128,64,20)))
    CE     = zeros(((128,64,20+1)))
    CF     = zeros(((128,64,20+1)))
    CFt    = zeros(((128,64,20+1)))
    

    for k in 1:19 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CA[:,:,k]   .= (rpdel[:,:,k]   .* Δt .* grav.^2 .* K_E[:,:,k+1]  .* rho[:,:,k+1].^2 
                       ./ (grid_p_full[:,:,k+1] .- grid_p_full[:,:,k]))
        CC[:,:,k+1] .= (rpdel[:,:,k+1] .* Δt .* grav.^2 .* K_E[:,:,k+1]  .* rho[:,:,k+1].^2
                       ./ (grid_p_full[:,:,k+1] .- grid_p_full[:,:,k]))
    end
    # @info maximum(rpdel) ## OK
    # @info maximum(rho) ## OK
    # @info maximum(grid_p_full[:,:,18+1] .- grid_p_full[:,:,18]) # # OK
    
    CA[:,:,20]   .= 0.
    CC[:,:, 1]   .= 0.
    CE[:,:, 1]   .= 0.
    CE[:,:,21]   .= 0.
    CF[:,:,21]   .= 0.
    CFt[:,:,21]  .= 0.
    
    # @info minimum(CA)
    # @info minimum(CC)

    p0 = 100000.
    for k in 20:-1:1
        CE[:,:,k]    .= CC[:,:,k] ./ (1. .+ CA[:,:,k] .+ CC[:,:,k] .- CA[:,:,k] .* CE[:,:,k+1])
        CF[:,:,k]    .= ((grid_tracers_c[:,:,k] .+ CA[:,:,k] .* CF[:,:,k+1])
                        ./ (1. .+ CA[:,:,k] .+ CC[:,:,k] .- CA[:,:,k] .* CE[:,:,k+1]))

        CFt[:,:,k]   .= (((p0./grid_p_full[:,:,k]).^(Rd/cp).*grid_t[:,:,k] .+ CA[:,:,k] .* CFt[:,:,k+1])
                        ./ (1. .+ CA[:,:,k] .+ CC[:,:,k] .- CA[:,:,k] .* CE[:,:,k+1]))
    end
    # @info maximum(CE)
    # @info maximum(CF)

    # first calculate the updates at the top model level
    # grid_δtracers[:,:,1] .+= (CF[:,:,1] .- grid_tracers_c[:,:,1]) ./ (2. .* Δt)
    ### WARNING factor1 just factor, so it did  ./ ./ (2. .* Δt). 
    ### So did factor2
    factor2[:,:,1]        .= (CF[:,:,1] .- grid_tracers_c[:,:,1]) ./ (2. .* Δt)  # because CE at top = 0
    grid_tracers_c[:,:,1] .= CF[:,:,1] 
    ##########################################################################################
    # grid_δt[:,:,1]   .+= (CFt[:,:,1] .* (grid_p_full[:,:,1]./p0).^(Rd/cp) .- grid_t[:,:,1]) ./ (2. .* Δt)
    grid_t[:,:,1]     .= (CFt[:,:,1] .* (grid_p_full[:,:,1]./p0).^(Rd/cp))

    
    # Loop over the remaining level
    for k in 2:20
        # grid_δtracers[:,:,k]  .+= (CE[:,:,k] .* grid_tracers_c[:,:,k-1] .+ CF[:,:,k] .- grid_tracers_c[:,:,k]) ./ (2. .* Δt)
        factor2[:,:,k]         .= (CE[:,:,k] .* grid_tracers_c[:,:,k-1] .+ CF[:,:,k] .- grid_tracers_c[:,:,k]) ./ (2. .* Δt)
        grid_tracers_c[:,:,k]  .=  CE[:,:,k] .* grid_tracers_c[:,:,k-1] .+ CF[:,:,k]
    #######################################################################################
        # grid_δt[:,:,k]    .+= ((CE[:,:,k] .* grid_t[:,:,k-1] .* (p0./grid_p_full[:,:,k-1]).^(Rd/cp) .+ CFt[:,:,k]) .* (grid_p_full[:,:,k]./p0).^(Rd/cp) .- grid_t[:,:,k]) ./ (2. .* Δt)
        grid_t[:,:,k]      .= ((CE[:,:,k] .* grid_t[:,:,k-1] .* (p0./grid_p_full[:,:,k-1]).^(Rd/cp) .+ CFt[:,:,k]) .* (grid_p_full[:,:,k]./p0).^(Rd/cp))
    end
    # @info maximum(CE), minimum(CE)
    # @info maximum(CF), minimum(CF)
end

# function Latent_heat_flux!
#    grid_tracers_c_ps_max           = zeros(((128,64,1)))
#    # grid_tracers_c_ps_max          .= (0.622 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t[:,:,20])) )) ./ (grid_ps[:,:,1] .- 0.378 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ grid_t[:,:,20])))) 
#    grid_tracers_c_ps_max          .= (0.622 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ Tsurf[:,:])) )) ./ (grid_ps[:,:,1] .- 0.378 .* (611.12 .* exp.(Lv ./ Rv .* (1. ./ 273.15 .- 1. ./ Tsurf[:,:])))) 

#    grid_δtracers[:,:,20] .+= (((C_E .* V_c[:,:,20] .* (grid_tracers_c_ps_max[:,:,1] .- min.(grid_tracers_c[:,:,20], grid_tracers_c_ps_max[:,:,1])) .* (2. * Δt) ./ za[:,:,1])
#                          ./ (1. .+ C_E .* V_c[:,:,20] .* (2. * Δt) ./ za[:,:,1])) ./ (2. * Δt))
#    # @info "max: ", maximum(grid_δtracers[:,:,20])
#    # @info "min: ", minimum(grid_δtracers[:,:,20])

#    grid_tracers_c[:,:,20]  .= ((grid_tracers_c[:,:,20] .+ C_E .* V_c[:,:,20] .* max.(grid_tracers_c[:,:,20],grid_tracers_c_ps_max[:,:,1]) .* (2. * Δt) ./ za[:,:,1]) 
#                          ./ (1. .+ C_E .* V_c[:,:,20] .* (2. * Δt) ./ za[:,:,1]))
# end

# function Spectral_Dynamics_Main()
#   # the decay of a sinusoidal disturbance to a zonally symmetric flow 
#   # that resembles that found in the upper troposphere in Northern winter.
#   name = "Spectral_Dynamics"
#   #num_fourier, nθ, nd = 63, 96, 20
#   num_fourier, nθ, nd = 42, 64, 20
#   #num_fourier, nθ, nd = 21, 32, 20
#   num_spherical = num_fourier + 1
#   nλ = 2nθ
  
#   radius = 6371000.0
#   omega = 7.292e-5
#   sea_level_ps_ref = 1.0e5
#   init_t = 264.0
  
#   # Initialize mesh
#   mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nθ, nλ, nd, radius)
#   θc, λc = mesh.θc,  mesh.λc
#   cosθ, sinθ = mesh.cosθ, mesh.sinθ
  
#   vert_coord = Vert_Coordinate(nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts", sea_level_ps_ref)
#   # Initialize atmo_data
#   do_mass_correction = true
#   do_energy_correction = true
#   do_water_correction = false
  
#   use_virtual_temperature = false
#   atmo_data = Atmo_Data(name, nλ, nθ, nd, do_mass_correction, do_energy_correction, do_water_correction, use_virtual_temperature, sinθ, radius,  omega)
  
#   # Initialize integrator
#   damping_order = 4
#   damping_coef = 1.15741e-4
#   robert_coef  = 0.04 
  
#   implicit_coef = 0.5
#   day_to_sec = 86400
#   start_time = 0
#   end_time = 2*day_to_sec  #
#   Δt = 1200
#   init_step = true
  
#   integrator = Filtered_Leapfrog(robert_coef, 
#   damping_order, damping_coef, mesh.laplacian_eig,
#   implicit_coef, Δt, init_step, start_time, end_time)
  
#   ps_ref = sea_level_ps_ref
#   t_ref = fill(300.0, nd)
#   wave_numbers = mesh.wave_numbers
#   semi_implicit = Semi_Implicit_Solver(vert_coord, atmo_data,
#   integrator, ps_ref, t_ref, wave_numbers)
  
#   # Initialize data
#   dyn_data = Dyn_Data(name, num_fourier, num_spherical, nλ, nθ, nd)
  
  
#   NT = Int64(end_time / Δt)
  
#   Get_Topography!(dyn_data.grid_geopots)
  
#   Spectral_Initialize_Fields!(mesh, atmo_data, vert_coord, sea_level_ps_ref, init_t,
#   dyn_data.grid_geopots, dyn_data)
  

#   Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data)

#   Update_Init_Step!(semi_implicit)
#   integrator.time += Δt 
#   for i = 2:NT

#     Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data)

#     integrator.time += Δt
#     @info integrator.time

#   end
  
# end


# #Spectral_Dynamics_Main()