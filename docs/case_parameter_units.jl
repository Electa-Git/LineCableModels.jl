# Units of the declared nominal inputs, before conversion by a case builder.
# Keep this explicit: for example, core_r20 is in Ω/km while core_rho is in Ω·m.
const CASE_PARAMETER_UNITS = Dict(
    [id => unit for (unit, ids) in (
        "m" => (
            :aluminum_foil_thickness, :aluminum_sheath_thickness,
            :aluminum_tape_thickness, :armor_wire_diameter, :bedding_thickness,
            :binder_thickness, :cable_x, :cable_y, :copper_tape_thickness,
            :copper_tape_width, :core_diameter, :core_outer_radius, :core_radius,
            :core_strand_diameter, :core_strand_radius, :core_wire_radius,
            :first_x, :formation_clearance, :formation_x, :formation_y,
            :inner_pe_thickness, :inner_semicon_thickness, :inner_sheath_thickness,
            :insulation_thickness, :jacket_thickness, :lead_screen_thickness,
            :lead_sheath_thickness, :line_length, :outer_semicon_skin,
            :outer_semicon_thickness, :pe_face_thickness, :post_sheath_filler,
            :pre_sheath_filler, :screen_inner_bedding_thickness,
            :screen_outer_bedding_thickness, :screen_wire_diameter,
            :screen_wire_radius, :second_x, :semicon_tape_thickness,
            :sheath_thickness, :strand_diameter, :water_blocking_thickness,
        ),
        "m²" => (:core_cross_section,),
        "Ω·m" => (
            :aluminum_rho, :core_rho, :earth_rho, :lead_rho, :pe_rho, :pp_rho,
            :semicon_rho, :steel_rho, :xlpe_rho,
        ),
        "Ω/km" => (:core_r20,),
        "kV" => (:voltage,),
        "°C" => (:temperature,),
        "Hz" => (:frequencies,),
        "dimensionless" => (
            :armor_lay_ratio, :armor_packing_clearance_ratio,
            :copper_tape_lay_ratio, :core_fillet_factor, :core_lay_ratio,
            :core_mu_r, :core_ring_lay_ratio_1, :core_ring_lay_ratio_2,
            :core_ring_lay_ratio_3, :core_ring_lay_ratio_4, :earth_eps_r,
            :formation_clearance_ratio, :lead_mu_r, :pe_eps_r, :pp_eps_r,
            :pp_mu_r, :screen_lay_ratio, :screen_wire_lay_ratio, :semicon_eps_r,
            :steel_mu_r, :xlpe_eps_r, :xlpe_mu_r, :xlpe_tan_delta,
        ),
        "dimensionless (count)" => (
            :armor_wires, :core_courses, :core_layers, :core_ring_wire_counts,
            :core_sectors, :ring_counts, :screen_wires, :strand_layers,
            :strands_per_layer,
        ),
    ) for id in ids]
)

function case_parameter_unit(id::Symbol)
    haskey(CASE_PARAMETER_UNITS, id) ||
        error("document the nominal input unit for case parameter :$id")
    return CASE_PARAMETER_UNITS[id]
end
