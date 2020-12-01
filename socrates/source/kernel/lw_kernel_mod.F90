!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to Socrates for longwave (thermal) fluxes

module lw_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_INTEGER,      &
                              GH_READ, GH_WRITE,         &
                              GH_READWRITE, CELLS,       &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3
use fs_continuity_mod, only:  W3, Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

public :: lw_kernel_type
public :: lw_code

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: lw_kernel_type
  private
  type(arg_type) :: meta_args(32) = (/                             &
    arg_type(GH_FIELD,   GH_WRITE,     Wtheta),                    & ! lw_heating_rate
    arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! lw_down_surf
    arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_2), & ! lw_up_tile
    arg_type(GH_FIELD,   GH_READWRITE, Wtheta),                    & ! lw_heating_rate_rts
    arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! lw_down_surf_rts
    arg_type(GH_FIELD,   GH_READWRITE, ANY_DISCONTINUOUS_SPACE_2), & ! lw_up_tile_rts
    arg_type(GH_FIELD,   GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1), & ! cloud_cover_rts
    arg_type(GH_FIELD,   GH_WRITE,     Wtheta),                    & ! cloud_fraction_rts
    arg_type(GH_FIELD,   GH_WRITE,     Wtheta),                    & ! cloud_droplet_re_rts
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! theta
    arg_type(GH_FIELD,   GH_READ,      W3),                        & ! theta_in_w3
    arg_type(GH_FIELD,   GH_READ,      W3),                        & ! exner
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! exner_in_wth
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! rho_in_wth
    arg_type(GH_FIELD,   GH_READ,      W3),                        & ! height_w3
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! height_wth
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! ozone
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! mv
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! mcl
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! mci
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! area_fraction
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! liquid_fraction
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! ice_fraction
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! cca
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! ccw
    arg_type(GH_FIELD,   GH_READ,      Wtheta),                    & ! cloud_drop_no_conc
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! tile_temperature
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_3), & ! tile_lw_albedo
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! latitude
    arg_type(GH_FIELD,   GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! longitude
    arg_type(GH_INTEGER, GH_READ                                )  & ! timestep
    /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: lw_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

! @param[in]     nlayers                 Number of layers
! @param[out]    lw_heating_rate         LW heating rate
! @param[out]    lw_down_surf            LW downward flux at the surface
! @param[out]    lw_up_tile              LW upward tiled surface flux
! @param[in,out] lw_heating_rate_rts     LW heating rate
! @param[in,out] lw_down_surf_rts        LW downward flux at the surface
! @param[in,out] lw_up_tile_rts          LW upward tiled surface flux
! @param[out]    cloud_cover_rts         Total cloud cover 2D field
! @param[out]    cloud_fraction_rts      Total cloud fraction 3D field
! @param[out]    cloud_droplet_re_rts    Cloud droplet effective radius 3D field
! @param[in]     theta                   Potential temperature
! @param[in]     theta_in_w3             Potential temperature in density space
! @param[in]     exner                   Exner pressure in density space
! @param[in]     exner_in_wth            Exner pressure in wth space
! @param[in]     rho_in_wth              Density in potential temperature space
! @param[in]     height_w3               Height of w3 levels above surface
! @param[in]     height_wth              Height of wth levels above surface
! @param[in]     ozone                   Ozone field
! @param[in]     mv                      Water vapour field
! @param[in]     mcl                     Cloud liquid field
! @param[in]     mci                     Cloud ice field
! @param[in]     area_fraction           Total cloud area fraction field
! @param[in]     liquid_fraction         Liquid cloud fraction field
! @param[in]     ice_fraction            Ice cloud fraction field
! @param[in]     cca                     Convective cloud amount (fraction)
! @param[in]     ccw                     Convective cloud water (kg/kg) (can be ice or liquid)
! @param[in]     cloud_drop_no_conc      Cloud Droplet Number Concentration
! @param[in]     tile_fraction           Surface tile fractions
! @param[in]     tile_temperature        Surface tile temperature
! @param[in]     tile_lw_albedo          LW tile albedos
! @param[in]     latitude                Latitude field
! @param[in]     longitude               Longitude field
! @param[in]     timestep                Timestep number
! @param[in]     ndf_wth                 No. DOFs per cell for wth space
! @param[in]     undf_wth                No. unique of DOFs for wth space
! @param[in]     map_wth                 Dofmap for wth space column base cell
! @param[in]     ndf_2d                  No. of DOFs per cell for 2D space
! @param[in]     undf_2d                 No. unique of DOFs for 2D space
! @param[in]     map_2d                  Dofmap for 2D space column base cell
! @param[in]     ndf_tile                Number of DOFs per cell for tiles
! @param[in]     undf_tile               Number of total DOFs for tiles
! @param[in]     map_tile                Dofmap for tile space column base cell
! @param[in]     ndf_w3                  No. of DOFs per cell for w3 space
! @param[in]     undf_w3                 No. unique of DOFs for w3 space
! @param[in]     map_w3                  Dofmap for w3 space column base cell
! @param[in]     ndf_rtile               No. of DOFs per cell for rtile space
! @param[in]     undf_rtile              No. unique of DOFs for rtile space
! @param[in]     map_rtile               Dofmap for rtile space column base cell
subroutine lw_code(nlayers,                          &
                   lw_heating_rate,                  &
                   lw_down_surf,                     &
                   lw_up_tile,                       &
                   lw_heating_rate_rts,              &
                   lw_down_surf_rts,                 &
                   lw_up_tile_rts,                   &
                   cloud_cover_rts,                  &
                   cloud_fraction_rts,               &
                   cloud_droplet_re_rts,             &
                   theta,                            &
                   theta_in_w3,                      &
                   exner,                            &
                   exner_in_wth,                     &
                   rho_in_wth,                       &
                   height_w3,                        &
                   height_wth,                       &
                   ozone,                            &
                   mv,                               &
                   mcl,                              &
                   mci,                              &
                   area_fraction,                    &
                   liquid_fraction,                  &
                   ice_fraction,                     &
                   cca, ccw,                         &
                   cloud_drop_no_conc,               &
                   tile_fraction,                    &
                   tile_temperature,                 &
                   tile_lw_albedo,                   &
                   latitude, longitude,              &
                   timestep,                         &
                   ndf_wth, undf_wth, map_wth,       &
                   ndf_2d, undf_2d, map_2d,          &
                   ndf_tile, undf_tile, map_tile,    &
                   ndf_w3, undf_w3, map_w3,          &
                   ndf_rtile, undf_rtile, map_rtile)

  use well_mixed_gases_config_mod, only:         &
    co2_mix_ratio, n2o_mix_ratio, ch4_mix_ratio, &
    cfc11_mix_ratio, cfc12_mix_ratio,            &
    cfc113_mix_ratio, hcfc22_mix_ratio, hfc134a_mix_ratio
  use radiation_config_mod, only: n_radstep,  &
    l_planet_grey_surface, planet_emissivity, &
    i_cloud_ice_type_lw, i_cloud_liq_type_lw, &
    cloud_horizontal_rsd, cloud_vertical_decorr
  use set_thermodynamic_mod, only: set_thermodynamic
  use set_cloud_field_mod, only: set_cloud_field
  use jules_control_init_mod, only: n_surf_tile, first_sea_tile
  use socrates_init_mod, only: n_lw_band
  use socrates_runes, only: runes, ip_source_thermal
  use socrates_bones, only: bones

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, timestep
  integer(i_def), intent(in) :: ndf_wth, ndf_w3, ndf_2d
  integer(i_def), intent(in) :: ndf_tile, ndf_rtile
  integer(i_def), intent(in) :: undf_wth, undf_w3, undf_2d
  integer(i_def), intent(in) :: undf_tile, undf_rtile

  integer(i_def), dimension(ndf_wth),   intent(in) :: map_wth
  integer(i_def), dimension(ndf_w3),    intent(in) :: map_w3
  integer(i_def), dimension(ndf_2d),    intent(in) :: map_2d
  integer(i_def), dimension(ndf_tile),  intent(in) :: map_tile
  integer(i_def), dimension(ndf_rtile), intent(in) :: map_rtile

  real(r_def), dimension(undf_wth),  intent(out) :: lw_heating_rate
  real(r_def), dimension(undf_2d),   intent(out) :: lw_down_surf
  real(r_def), dimension(undf_tile), intent(out) :: lw_up_tile

  real(r_def), dimension(undf_wth),  intent(inout) :: lw_heating_rate_rts
  real(r_def), dimension(undf_2d),   intent(inout) :: lw_down_surf_rts
  real(r_def), dimension(undf_tile), intent(inout) :: lw_up_tile_rts

  real(r_def), dimension(undf_2d),  intent(out) :: cloud_cover_rts
  real(r_def), dimension(undf_wth), intent(out) :: cloud_fraction_rts
  real(r_def), dimension(undf_wth), intent(out) :: cloud_droplet_re_rts

  real(r_def), dimension(undf_w3),  intent(in) :: theta_in_w3, exner, height_w3
  real(r_def), dimension(undf_wth), intent(in) :: theta, exner_in_wth, &
    rho_in_wth, height_wth, ozone, mv, mcl, mci, &
    area_fraction, liquid_fraction, ice_fraction, cca, ccw, cloud_drop_no_conc
  real(r_def), dimension(undf_tile),  intent(in) :: tile_fraction
  real(r_def), dimension(undf_tile),  intent(in) :: tile_temperature
  real(r_def), dimension(undf_rtile), intent(in) :: tile_lw_albedo
  real(r_def), dimension(undf_2d),    intent(in) :: latitude, longitude

  ! Local variables for the kernel
  integer(i_def), parameter :: n_profile = 1
  integer(i_def) :: i_cloud_representation, i_overlap, i_inhom, i_drop_re
  integer(i_def) :: rand_seed(n_profile)
  integer(i_def) :: k, n_cloud_layer
  integer(i_def) :: df_rtile, i_tile, i_band
  integer(i_def) :: wth_1, wth_nlayers
  real(r_def), dimension(nlayers) :: layer_heat_capacity
    ! Heat capacity for each layer
  real(r_def), dimension(nlayers) :: p_layer, t_layer
    ! Layer pressure and temperature
  real(r_def), dimension(nlayers) :: d_mass
    ! Mass of layer per square metre
  real(r_def), dimension(nlayers) :: cloud_frac, liq_dim
    ! Large-scale cloud fields
  real(r_def), dimension(nlayers) :: conv_frac, &
    liq_conv_frac, ice_conv_frac, liq_conv_mmr, ice_conv_mmr, liq_conv_dim
    ! Convective cloud fields
  real(r_def), dimension(0:nlayers) :: t_layer_boundaries
  real(r_def), dimension(0:nlayers) :: lw_down, lw_up

  ! Tiled surface fields
  real(r_def) :: frac_tile(n_profile, n_surf_tile)
  real(r_def) :: t_tile(n_profile, n_surf_tile)
  real(r_def) :: albedo_diff_tile(n_profile, n_surf_tile, n_lw_band)
  real(r_def) :: flux_up_tile_rts(n_surf_tile)
  real(r_def) :: flux_up_tile(n_surf_tile)


  ! Set wth indexing
  wth_1 = map_wth(1)+1
  wth_nlayers = map_wth(1)+nlayers

  ! Tile fractions & temperatures
  do i_tile = 1, n_surf_tile
    frac_tile(1, i_tile) = tile_fraction(map_tile(i_tile))
    t_tile(1, i_tile) = tile_temperature(map_tile(i_tile))
  end do

  ! Tile albedos
  df_rtile = 0
  do i_band = 1, n_lw_band
    do i_tile = 1, n_surf_tile
      df_rtile = df_rtile + 1
      albedo_diff_tile(1, i_tile, i_band) &
        = tile_lw_albedo(map_rtile(1)+df_rtile-1)
    end do
  end do

  if (mod(timestep-1_i_def, n_radstep) == 0) then
    ! Radiation time-step: full calculation of radiative fluxes

    ! Set up pressures, temperatures, masses and heat capacities
    call set_thermodynamic(nlayers,                &
      exner(map_w3(1):map_w3(1)+nlayers-1),        &
      exner_in_wth(map_wth(1):map_wth(1)+nlayers), &
      theta(map_wth(1):map_wth(1)+nlayers),        &
      rho_in_wth(map_wth(1):map_wth(1)+nlayers),   &
      height_w3(map_w3(1):map_w3(1)+nlayers-1),    &
      height_wth(map_wth(1):map_wth(1)+nlayers),   &
      p_layer, t_layer, d_mass, layer_heat_capacity)

    ! Calculate temperature at layer boundaries
    t_layer_boundaries(0) = theta(map_wth(1)) * exner_in_wth(map_wth(1))
    do k=1,nlayers-1
      t_layer_boundaries(k) = theta_in_w3(map_w3(1)+k) * exner(map_w3(1)+k)
    end do
    t_layer_boundaries(nlayers) = t_layer(nlayers)

    ! Set up cloud fields for radiation
    call set_cloud_field(nlayers, n_profile, &
      area_fraction(wth_1:wth_nlayers), &
      cca(wth_1:wth_nlayers), ccw(wth_1:wth_nlayers), t_layer, &
      latitude(map_2d(1):map_2d(1)), longitude(map_2d(1):map_2d(1)), &
      i_cloud_representation, i_overlap, i_inhom, i_drop_re, &
      rand_seed, n_cloud_layer, cloud_frac, liq_dim, conv_frac, &
      liq_conv_frac, ice_conv_frac, liq_conv_mmr, ice_conv_mmr, liq_conv_dim)

    ! Calculate the LW fluxes (RUN the Edwards-Slingo two-stream solver)
    call runes(n_profile, nlayers,                                             &
      spectrum_name           = 'lw',                                          &
      i_source                = ip_source_thermal,                             &
      n_cloud_layer           = n_cloud_layer,                                 &
      p_layer_1d              = p_layer,                                       &
      t_layer_1d              = t_layer,                                       &
      mass_1d                 = d_mass,                                        &
      density_1d              = rho_in_wth(wth_1:wth_nlayers),                 &
      t_level_1d              = t_layer_boundaries,                            &
      h2o_1d                  = mv(wth_1:wth_nlayers),                         &
      o3_1d                   = ozone(wth_1:wth_nlayers),                      &
      co2_mix_ratio           = co2_mix_ratio,                                 &
      n2o_mix_ratio           = n2o_mix_ratio,                                 &
      ch4_mix_ratio           = ch4_mix_ratio,                                 &
      cfc11_mix_ratio         = cfc11_mix_ratio,                               &
      cfc12_mix_ratio         = cfc12_mix_ratio,                               &
      cfc113_mix_ratio        = cfc113_mix_ratio,                              &
      hcfc22_mix_ratio        = hcfc22_mix_ratio,                              &
      hfc134a_mix_ratio       = hfc134a_mix_ratio,                             &
      t_ground                = t_tile(1, first_sea_tile),                     &
      l_grey_albedo           = l_planet_grey_surface,                         &
      grey_albedo             = 1.0_r_def - planet_emissivity,                 &
      l_tile                  = .not.l_planet_grey_surface,                    &
      n_tile                  = n_surf_tile,                                   &
      frac_tile               = frac_tile,                                     &
      t_tile                  = t_tile,                                        &
      albedo_diff_tile        = albedo_diff_tile,                              &
      cloud_frac_1d           = cloud_frac,                                    &
      liq_frac_1d             = liquid_fraction(wth_1:wth_nlayers),            &
      ice_frac_1d             = ice_fraction(wth_1:wth_nlayers),               &
      liq_mmr_1d              = mcl(wth_1:wth_nlayers),                        &
      ice_mmr_1d              = mci(wth_1:wth_nlayers),                        &
      liq_dim_1d              = liq_dim,                                       &
      liq_nc_1d               = cloud_drop_no_conc(wth_1:wth_nlayers),         &
      conv_frac_1d            = conv_frac,                                     &
      liq_conv_frac_1d        = liq_conv_frac,                                 &
      ice_conv_frac_1d        = ice_conv_frac,                                 &
      liq_conv_mmr_1d         = liq_conv_mmr,                                  &
      ice_conv_mmr_1d         = ice_conv_mmr,                                  &
      liq_conv_dim_1d         = liq_conv_dim,                                  &
      liq_conv_nc_1d          = cloud_drop_no_conc(wth_1:wth_nlayers),         &
      cloud_horizontal_rsd    = cloud_horizontal_rsd,                          &
      cloud_vertical_decorr   = cloud_vertical_decorr,                         &
      conv_vertical_decorr    = cloud_vertical_decorr,                         &
      rand_seed               = rand_seed,                                     &
      layer_heat_capacity_1d  = layer_heat_capacity,                           &
      l_mixing_ratio          = .true.,                                        &
      i_cloud_representation  = i_cloud_representation,                        &
      i_overlap               = i_overlap,                                     &
      i_inhom                 = i_inhom,                                       &
      i_drop_re               = i_drop_re,                                     &
      i_st_water              = i_cloud_liq_type_lw,                           &
      i_st_ice                = i_cloud_ice_type_lw,                           &
      i_cnv_water             = i_cloud_liq_type_lw,                           &
      i_cnv_ice               = i_cloud_ice_type_lw,                           &
      l_invert                = .true.,                                        &
      total_cloud_cover       = cloud_cover_rts(map_2d(1):map_2d(1)),          &
      total_cloud_fraction_1d = cloud_fraction_rts(wth_1:wth_nlayers),         &
      liq_dim_diag_1d         = cloud_droplet_re_rts(wth_1:wth_nlayers),       &
      flux_down_1d            = lw_down,                                       &
      flux_up_1d              = lw_up,                                         &
      flux_up_tile_1d         = flux_up_tile_rts,                              &
      heating_rate_1d         = lw_heating_rate_rts(wth_1:wth_nlayers))

    ! Set level 0 increment such that theta increment will equal level 1
    lw_heating_rate_rts(map_wth(1)) = lw_heating_rate_rts(map_wth(1) + 1) &
                                    * exner_in_wth(map_wth(1))            &
                                    / exner_in_wth(map_wth(1) + 1)

    ! Set surface flux
    lw_down_surf_rts(map_2d(1)) = lw_down(0)

    ! Set tiled fluxes
    do i_tile = 1, n_surf_tile
      lw_up_tile_rts(map_tile(i_tile)) = flux_up_tile_rts(i_tile)
    end do

    ! No corrections needed for model timestep
    lw_heating_rate(map_wth(1):wth_nlayers) &
      = lw_heating_rate_rts(map_wth(1):wth_nlayers)

    lw_down_surf(map_2d(1)) = lw_down_surf_rts(map_2d(1))

    do i_tile = 1, n_surf_tile
      lw_up_tile(map_tile(i_tile)) = lw_up_tile_rts(map_tile(i_tile))
    end do

  else
    ! Not a radiation time-step: apply corrections to the model timestep

    ! The bare "bones" of a simple radiative transfer calculation.
    ! Update the upward fluxes for the change in surface temperature.
    call bones(n_profile, nlayers,                                             &
      n_tile                 = n_surf_tile,                                    &
      l_grey_emis_correction = .true.,                                         &
      grey_albedo_tile       = albedo_diff_tile(:, :, 1),                      &
      t_tile                 = t_tile,                                         &
      heating_rate_1d_rts    = lw_heating_rate_rts(wth_1:wth_nlayers),         &
      flux_down_surf_rts     = lw_down_surf_rts(map_2d(1):map_2d(1)),          &
      heating_rate_1d_mts    = lw_heating_rate(wth_1:wth_nlayers),             &
      flux_up_tile_1d_mts    = flux_up_tile,                                   &
      flux_down_surf_mts     = lw_down_surf(map_2d(1):map_2d(1)))

    ! Set level 0 increment such that theta increment will equal level 1
    lw_heating_rate(map_wth(1)) = lw_heating_rate(map_wth(1) + 1) &
                                * exner_in_wth(map_wth(1))        &
                                / exner_in_wth(map_wth(1) + 1)

    ! Set tiled fluxes
    do i_tile = 1, n_surf_tile
      lw_up_tile(map_tile(i_tile)) = flux_up_tile(i_tile)
    end do

  end if

end subroutine lw_code

end module lw_kernel_mod
