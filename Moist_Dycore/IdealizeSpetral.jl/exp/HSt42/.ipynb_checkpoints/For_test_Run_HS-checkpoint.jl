using JGCM

include("HS.jl")


end_day = 300
spinup_day = 0

L = 0 / 100
warm_start_file_name = "None"
initial_day = 0 # warm start day

physics_params = Dict{String,Float64}("σ_b"=>0.7, "k_f" => 1.0, "k_a" => 1.0/40.0, "k_s" => 1.0/4.0, "ΔT_y" => 60.0, "Δθ_z" => 10.0) ### now 65.0!!!
op_man = Atmos_Spectral_Dynamics_Main(physics_params, end_day, spinup_day, L)

# Finalize_Output!(op_man, "200day_convection_all.dat", "200day_convection_final.dat")
# Finalize_Output!(op_man, "15day_test_all.dat", "15day_test_final.dat")
# Finalize_Output!(op_man, "10day_test_all.dat", "10day_test_final.dat")
# Finalize_Output!(op_man, "test_all.dat", "test_final.dat")
# Finalize_Output!(op_man, "50day_test_all.dat", "50day_test_final.dat")
Finalize_Output!(op_man, "300day_PR10_for_check_warmstart_all.dat", "300day_PR10_for_check_warmstart_final.dat")



# Finalize_Output!(op_man, "200day_edit_factor12_all.dat", "200day_edit_factor12_final.dat")
# Finalize_Output!(op_man, "200day_edit_factor12_delta_t_all.dat", "200day_edit_factor12_delta_t_final.dat")




