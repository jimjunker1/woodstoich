extensions [sr palette]

breed [cluster-makers cluster-maker]
breed [resource-makers resource-maker]

;;===========================================================================================================================================================================================;;
;;=========================================================== VARIABLE/PARAMETER DEFINITIONS (Netlogo code formality/requirement) ===========================================================;;
;;===========================================================================================================================================================================================;;
globals [
  ;;===== World globals
  world-area                                  ; Total number of patches  [# patches]
  world-avg-resource-CN-ratio                 ; Average carbon:nitrogen ratio for resources in entire world  [#]
  world-avg-resource-NC-ratio                 ; Average nitrogen:carbon ratio for resources in entire world  [#]
  world-total-initial-resource-N              ; Total amount of nitrogen in the world, constrained in closed system  [mol N]
  world-resource-stoichmismatch-severity      ; Factor to increase the carbon:nitrogen ratio in the world by, based on user input  [#]
  realized-avg-world-resource-NC-ratio        ; Calculated average nitrogen:carbon ratio in the landscape (post-setup)  [#]
  realized-avg-world-resource-CN-ratio        ; Calculated average carbon:nitrogen ratio in the landscape (post-setup)  [#]
  avg-resource-quantity                       ; Calculated average of resource amount in landscape  [mol Food]
  sum-resource-N-total                        ; Total amount of N in resources variables across landscape  [mol N]
  sum-world-N-total                           ; Total N in the world  [mol N]
  sum-consumer-N-total                        ; Total N in all consumer biomass  [mol N]
  sum-consumer-N-total_adult                  ;
  sum-consumer-N-total_juv                    ;
  sum-soil-N-total                            ; Total N in all soil pools  [mol N]
  total-world-N                               ; All biomass pools  [mol N]
  fraction-world-resources                    ; Proportion of patches that can contain resources  [#: 0 to 1]
  resource-patches                            ; Designation for patches that can grow resources [NA]
  no-resource-patches                         ; Designation for patches that can not grow resources  [NA]
  cluster-coords                              ; Coordinates for the start of a high-quality resource cluster  [#]
  num-coords                                  ; Number of coordinates for the start of a high-quality resource cluster  [#]

  ;;===== Patch related globals
  ;;======== Resource patch cluster globals
  resource-cluster-center-coords              ; The xy-coordinates of center of resource clusters  [#]
  current-per-resource-cluster-area           ; Calculated actual area of a resource cluster  [# patches]
  total-resource-cluster-patches              ; Number of patches that are part of a resource cluster  [# patches]
  absolute-max-cluster-resource-area          ; Calculated absolute maximum area of world that can be clusters  [# patches]
  max-cluster-resource-area                   ; Possible maximum number of patches dedicated to resource clusters, estimated from how to maintain world average carbon:nitrogen ratio  [# patches]
  noncluster-resource-area                    ; Minimum possible total patches that are allocated for non-clustered resources (low quality)  [# patches]
  max-per-cluster-area                        ; Maximum per-cluster area, based on max-cluster-resource-area and num-resource-clusters (number of total clusters, set by user interface)  [# patches]
  max-cluster-radius                          ; Maximum radius for each cluster  [# patches]
  realized-cluster-avg-resource-NC-ratio      ; Calculated average resource nitrogen:carbon ratio for cluster patches after setup is complete  [#]
  realized-cluster-avg-resource-CN-ratio      ; Calculated average resource carbon:nitrogen ratio for cluster patches after setup is complete  [#]
  target-cluster-avg-resource-CN-ratio        ; Goal average resource carbon:nitrogen ratio, used to setup patches  [#]
  target-cluster-avg-resource-NC-ratio        ; Goal average resource nitrogen:carbon ratio, used to setup patches  [#]
  target-noncluster-avg-resource-CN-ratio     ; Goal average resource carbon:nitrogen ratio, used to setup patches  [#]
  target-noncluster-avg-resource-NC-ratio     ; Goal average resource nitrogen:carbon ratio, used to setup patches  [#]
  min-cn-color-value                          ; Patch coloring  [#]
  max-cn-color-value                          ; Patch coloring  [#]

  ;;======== Resource patch NON cluster globals
  total-noncluster-patches                    ; Number of patches that are NOT part of a resource cluster  [#]
  noncluster-avg-resource-NC-ratio            ; Target average resource nitrogen:carbon ratio for NON cluster patches after setup is complete  [#]
  noncluster-avg-resource-CN-ratio            ; Target average resource carbon:nitrogen ratio for NON cluster patches after setup is complete  [#]
  realized-noncluster-avg-resource-NC-ratio   ; Calculated average resource nitrogen:carbon ratio for NONcluster patches after setup is complete  [#]
  realized-noncluster-avg-resource-CN-ratio   ; Calculated average resource carbon:nitrogen ratio for NONcluster patches after setup is complete  [#]
  resource-quantity-min                       ; A minimum value to set the resource quantity to  [#]
  total-resource-quantity                     ; Sum of resource quality  [mol Food]

  ;;===== Consumer related globals
  consumer-body-NC-ratio                      ; Nitrogen:carbon ratio of consumer's body, inverse of the carbon:nitrogen ratio provided in user interface  [#]
  total-consumers                             ; Sum of all turtles  [#]
  total-adult-consumers                       ; Sum of adult turtles  [#]
  total-juvenille-consumers                   ; Sum of turtles that aren't adults (babies/juveniles)  [#]
  avg-Creserve-scaleddensity                  ; Average scaled reserve density of turtles' C-reserve. For plotting  [#: 0 to 1]
  avg-Nreserve-scaleddensity                  ; Average scaled reserve density of turtles' N-reserve. For plotting  [#: 0 to 1]
  avg-M_V                                     ; Average turtle structure biomass (size). For plotting  [mol C]
  avg-M_R                                     ; Average turtle reproductive buffer biomass. For plotting  [mol C]
  avg-M_H                                     ; Average mass that has been invested into turtle maturity. For plotting  [mol C]
  Nmol_ODE_check_sum                          ; Simulation check variable to watch for any N-mol conservation issues arising out of the DEB model
  movement_response
  eatmore_response
]

patches-own [
  ;;===== Resource (primary production) patch variables
  resource-CN-ratio                           ; Carbon:nitrogen ratio of resource on patch  [#]
  resource-NC-ratio                           ; Nitrogen:carbon ratio of resource on patch  [#]
  resource-quantity                           ; Total amount of resource per patch, set by the user interface  [mol Food]
  resource-quantity-Cmol                      ; Total amount of resource carbon on a patch  [mol C]
  resource-quantity-Nmol                      ; Total amount of resource nitrogen on a patch  [mol N]
  new-resource-grown                          ; Increase in resource quantity per tick  [#]
  new-resource-quantity-Cmol                  ; Increase in resource mol C per tick  [mol C]
  new-resource-quantity-Nmol                  ; Increase in resource mol N per tick  [mol N]
  resource-growth-max                         ; Used as a ceiling to calculate how many resources can grow this tick  [#]
  resource-CN-ratio-test                      ; Used to test if projected C:N ratio will be above or below set thresholds

  ;;===== Patch binary info
  resource-patch?                             ; Binary assignment of if patch can grow resources  [0 or 1]
  cluster-center?                             ; Binary assignment of if patch is a center of a patch  [0 or 1]
  cluster-patch?                              ; Binary assignment of if patch is part of a patch  [0 or 1]

  ;;===== Soil pool
  soil-N-quantity                             ; Soil N pool size  [mol N]
  turtle-alterations-patch-pooled-N           ; How much N in moles is left by turtles on patch per turn  [mol N]
]

turtles-own [
  ;;===== Turtle binary info
  is-adult?                                   ; Binary assignment to indicate consumer lifestage  [0 or 1]
  did-baby?                                   ; Binary assignment for if consumer reproduced at end of tick  [0 or 1]
  starved-to-death?                           ; Binary assignment for if consumer starved, calculated by their reserve densities  [0 or 1]
  shrank-too-much?                            ; Binary assignment for if consumer starved, based on how much somatic structure they have broken down  [0 or 1]
  demographic-stochasticity-death?            ; Binary assignment for if consumer was randomly assigned to die at the end of the tick  [0 or 1]
  should-die?                                 ; Binary assignment for if consumer should die at the end of tick, based on any death criteria   [0 or 1]

  ;;===== Nutrient-related parameters
  turtle-input-N-this-tick                    ; Total amount of N in moles that is released into the patch at the end of turtle actions per tick  [mol N]
  consumer-total-Nmol                         ; The current amount of N moles in the turtle state variables that contribute towards its total composition  [mol N]
  corpse-C                                    ; Amount of C released to patch from consumer's corpse  [mol C]
  corpse-N                                    ; Amount of N released to patch from consumer's corpse  [mol N]

  ;;===== State Variables
  M_V                                         ; Somatic body biomass  [mol C]
  V                                           ; Somatic body volume  [cm^3]
  L                                           ; Length of body  [cm]
  largest-Mv-to-date                          ; Largest biomass reached to this tick for turtle  [mol C]
  M_Ec                                        ; Biomass of carbon-rich reserve  [mol C]
  M_En                                        ; Biomass of nitrogen-rich reserve  [mol C]
  Mdens_Ec                                    ; Carbon-rich reserve density  [mol C / mol C]
  scaled_Mdens_Ec                             ; Scaled carbon-rich reserve density  [# 0 - 1]
  Mdens_En                                    ; Nitrogen-rich reserve density  [mol C / mol C]
  scaled_Mdens_En                             ; Scaled nitrogen-rich reserve density  [#: 0 to 1]
  M_R                                         ; Biomass of reserves invested in reproduction  [mol C]
  M_H                                         ; Biomass of reserves invested in maturity/development  [mol C]
  Xc                                          ; Mol C of combined turtle excretion  [mol C]
  Xn                                          ; Mol N of combined turtle excretion  [mol N]
  Wc                                          ; Mol C of turtle egestion  [mol C]
  Wn                                          ; Mol N of turtle egestion  [mol N]

  ;;===== Mass Balance Differential Equations
  dFood-dt                                    ; Change in food density per tick calculation  [mol Food]
  food_ingested                               ; The amount that was changed above, but used for patch resource quantity calculations  [mol Food]
  dWc-dt                                      ; Change of carbon egestion per tick calculation  [mol C]
  dWn-dt                                      ; Change of nitrogen egestion per tick calculation  [mol N]
  dM_Ec-dt                                    ; Change of carbon-rich reserve biomass per tick calculation  [mol C]
  dM_En-dt                                    ; Change of nitrogen-rich reserve biomass per tick calculation  [mol C]
  dMdens_Ec-dt                                ; Change of carbon-rich reserve density per tick calculation  [mol C / mol C]
  dMdens_En-dt                                ; Change of nitrogen-rich reserve density per tick calculation  [mol C / mol C]
  dM_V-dt                                     ; Change of somatic structual biomass per tick calculation  [mol C]
  dM_H-dt                                     ; Change of development biomass per tick calculation  [mol C]
  dM_R-dt                                     ; Change of reproduction biomass per tick calculation  [mol C]
  dXc-dt                                      ; Change of excreted carbon biomass per tick calculation  [mol C]
  dXn-dt                                      ; Change of excreted nitrogen biomass per tick calculation  [mol N]
  ODE_Nconservation_test

  ;;===== Initial state variable conditions for juvenilles at birth
  M_V_baby                                    ; Starting amount of structural biomass for juvenille  [mol C]
  M_Ec_baby                                   ; Starting amount of carbon-rich reserve for juvenille  [mol C]
  M_En_baby                                   ; Starting amount of nitrogen-rich reserve for juvenille  [mol C]
  scaled_mdens_Ec_baby                        ; Starting scaled density of carbon-rich reserve for juvenille  [mol C / mol C]
  scaled_mdens_En_baby                        ; Starting scaled density of nitrogen-rich reserve for juvenille  [mol C / mol C]
  M_H_baby                                    ; Amount of maturity juvenille is born with  [mol C]

  ;;===== Maintenance parameters
  J_S                                         ; The cost of somatic maintenance this tick  [mol C]
  J_D                                         ; The cost of maturity maintenance this tick  [mol C]
  sigma                                       ; Factor increase of somatic maintenance requirement based on movement  [#]
  k_M                                         ; Somatic maintenance coefficient  [1 / day]
  k_J                                         ; Maturity maintenance coefficient  [1 / day]

  ;;===== Somatic structure/growth parameters
  r_dot                                       ; Specific growth rate of organism [1 / day]
  M_V_SqBrack                                 ; [M_V]: Volume-specific structural mass (M_V / V) [mol C / cm^3]

  ;;===== Reserve parameters
  kappa                                       ; Fraction of reserves sent to somatic vs reproductive processes, parameter given on Add-My-Pet DEB Portal  [#: 0 to 1]
  kappa_Ec                                    ; Fraction of C-rich reserves sent to somatic vs reproductive processes, user calculated  [#: 0 to 1]
  kappa_En                                    ; Fraction of N-rich reserves sent to somatic vs reproductive processes, user calculated  [#: 0 to 1]
  kappaG_Ec                                   ; Fraction of C-rich reserves rejected from the growth Syntheizing Unit that are recycled to the reserveto somatic vs reproductive processes, user calculated   [# 0 to 1]
  kappaG_En                                   ; Fraction of N-rich reserves rejected from the growth Syntheizing Unit that are recycled to the reserveto somatic vs reproductive processes, user calculated   [# 0 to 1]
  kappaR_Ec                                   ; Fraction of C-rich reserves rejected from the reproduction Syntheizing Unit that are recycled to the reserveto somatic vs reproductive processes, user calculated   [# 0 to 1]
  kappaR_En                                   ; Fraction of N-rich reserves rejected from the reproduction Syntheizing Unit that are recycled to the reserveto somatic vs reproductive processes, user calculated   [# 0 to 1]
  nu                                          ; Reserve energy conductance  [cm / day]
  mdens_tot_max                               ; Maximum total reserve density, parameter given on Add-My-Pet DEB Portal  [mol C / mol C]
  mdens_Ec_max                                ; Maximum carbon-rich reserve density, user calculated  [mol C / mol C]
  mdens_En_max                                ; Maximum nitrogen-rich reserve density, user calculated  [mol C / mol C]

  ;;===== Stoichiometry parameters
  q-F_C                                       ; Molar ratio of carbon per carbon in food  [#: 1]
  q-F_N                                       ; Molar ratio of nitrogen in food, normalized to carbon  [#: 0 to 1]
  q-Ec_C                                      ; Molar ratio of carbon per carbon in the carbon-rich reserve  [#: 1]
  q-Ec_N                                      ; Molar ratio of carbon in the nitrogen-rich reserve, normalized to carbon  [#: 0 to 1]
  q-En_C                                      ; Molar ratio of carbon in the nitrogen-rich reserve, normalized to carbon [#: 1]
  q-En_N                                      ; Molar ratio of nitrogen in the nitrogen-rich reserve, normalized to carbon[#: 0 to 1]
  q-V_C                                       ; Molar ratio of carbon per carbon in somatic structure  [#: 1]
  q-V_N                                       ; Molar ratio of nitrogen in somatic structure, normalized to carbon  [#: 0 to 1]
  q-R_C                                       ; Molar ratio of carbon per carbon in reproduction buffer  [#: 1]
  q-R_N                                       ; Molar ratio of nitrogen in reproduction buffer, normalized to carbon  [#: 0 to 1]


  ;;===== Yield coefficients and related parameters
  n-V_Ec                                      ; Number of mol carbon-rich needed to make one mol structure  [#]
  n-V_En                                      ; Number of mol nitrogen-rich needed to make one mol structure  [#]
  n-R_Ec                                      ; Number of mol carbon-rich needed to make one mol reproduction/maturity buffer  [#]
  n-R_En                                      ; Number of mol nitrogen-rich needed to make one mol reproduction/maturity buffer  [#]
  Y_EcF                                       ; Yield of food assimilated onto the carbon-rich reserve  [#]
  Y_EnF                                       ; Yield of food assimilated onto the nitrogen-rich reserve  [#]
  Y_EcV                                       ; Yield of carbon-rich reserve onto somatic growth  [#]
  Y_EnV                                       ; Yield of nitrogen-rich reserve onto somatic growth  [#]
  Y_EcR                                       ; Yield of carbon-rich reserve onto reproduction/maturity buffer  [#]
  Y_EnR                                       ; Yield of nitrogen-rich reserve onto reproduction/maturity buffer  [#]

  ;;===== Maturity and reproduction parameters
  M_H_p                                       ; Biomass needed in development to mature to reproductively active  [mol C]
  M_R_bufferbaby                              ; Amount of biomass needed in reproductive buffer to produce juvenille consumer  [mol C]
  M_R_bufferlitter                            ; Amount of biomass needed in reproductive buffer to produce a clutch of juvenille consumers  [mol C]
  avg_litter_size                             ; Number of individuals in a clutch  [#]

  ;;===== Ingestion/assimilation parameters
  J_XA_max_CrlBrack                           ; Surface-area specific maximum ingestion rate [mol Food / (cm^(2/3) * day)]
  F_h                                         ; Type II functional response half-saturation constant [mol Food]
  Food                                        ; Amount of resources available on patch  [mol Food]
  f                                           ; Type II functional response: (Food / (F_h + Food))  [#: 0 to 1]
  max-i-can-eat                               ; Parameter calculating the maximum food an individual could eat, if possible (for if more than one turtle is on the same patch)
  what-i-will-eat                             ; Parameter for the amount of food consumer will eat if it has to share food
  remaining-food                              ; Parameter to count how much food is left on a patch shared by turtles

  ;;===== Synthesizing Units parameters
  rho_S_N                                     ; Preference parameter for using nitrogen-rich reserve in meeting somatic maintenance costs, in comparison to using carbon-rich reserve  [#: 0 to 1]
  rho_S_V                                     ; Preference parameter for using structure in meeting somatic maintenance costs, in comparison to using carbon-rich reserve  [#: 0 to 1]
  rho_D_N                                     ; Preference parameter for using nitrogen-rich reserve in meeting maturity maintenance costs, in comparison to using carbon-rich reserve  [#: 0 to 1]
  rho_D_V                                     ; Preference parameter for using structure in meeting maturity maintenance costs, in comparison to using carbon-rich reserve  [#: 0 to 1]
  J_EcC_lastroot                              ; Storage of the last tick's root solution for the mobilization flux of carbon-rich reserve, from the numerical solving tool  [mol C / day]
  J_EcC_root                                  ; Storage of the current tick's root solution for the mobilization flux of carbon-rich reserve, from the numerical solving tool  [mol C / day]
  J_EnC_lastroot                              ; Storage of the last tick's root solution for the mobilization flux of nitrogen-rich reserve, from the numerical solving tool  [mol C / day]
  J_EnC_root                                  ; Storage of the current tick's root solution for the mobilization flux of nitrogen-rich reserve, from the numerical solving tool  [mol C / day]
  kS_lastroot                                 ; Storage of the last tick's root solution for the somatic maintenance disassociation parameter, from the numerical solving tool  [mol C / day]
  kS_root                                     ; Storage of the current tick's root solution for the somatic maintenance disassociation parameter, from the numerical solving tool  [mol C / day]
  kD_lastroot                                 ; Storage of the last tick's root solution for the maturity maintenance disassociation parameter, from the numerical solving tool  [mol C / day]
  kD_root                                     ; Storage of the current tick's root solution for the maturity maintenance disassociation parameter, from the numerical solving tool  [mol C / day]
  rdot_lastroot                               ; Storage of the last tick's root solution for the specific growth rate, from the numerical solving tool  [1 / day]
  rdot_root                                   ; Storage of the current tick's root solution for the specific growth rate, from the numerical solving tool  [1 / day]
]


;;==============================================================================================================================================;;
;;=========================================================== USER INTERFACE BUTTONS ===========================================================;;
;;==============================================================================================================================================;;


;;========================== SETUP: INITIALLIZE SIMULATIONS ==========================;;
;;====================================================================================;;

to setup
  ca
  clear-output

  ;;= Initialize the SimpleR Extension and load your R file into ABM
  sr:setup
  sr:run "source('C:/Users/jj0895/Documents/Projects/UNT_projects/woodstoich/code/DEBmodel_ODEsolver_for_NetLogo.R')"
;  sr:run "source('C:/Users/lynda/Documents/Work/Conferences/2024/Woodstoich/Model/NetLogo/DEBmodel_ODEsolver_for_NetLogo.R')"

  ;;= Procedures to set parameters/variables
  set world-area count patches  ; Create global that has total number of patches
  setup-parameters-from-interface

  ;;= Procedures to initialize patches
  setup-patches

  ;;= Initialize consumer population
  setup-turtles

  reset-ticks
end

;;=========================== GO: RUN PROCEDURES PER TICK ===========================;;
;;===================================================================================;;

to go
  ;;= Turtle procedures
  ask turtles [
    reset-ODE-and-fluxes
  ]
  new-tick-turtle-values                 ; Reset relevant turtle variables
  move                                   ; Movement for consumers
  hungry-hungry-hippos
  calculate-DEB-ODE                      ; Run bioenergetic individual-consumer model (Dynamic Energy Budget theory model)
  update-turtle-variables                ; Update the state variables and important parameters of consumer
  have-die?                              ; Determine if a consumer should die
  patch-N-added-by-turtles               ; Calculate the mol N that consumers will leave on the patch
  have-baby?                             ; Allow juvenilles to be born
  go-die                                 ; Kill consumers that were indicated to die
;  update-consumerN                       ; Re-update the total N of consumers after procedures

  ;;= Patch procedures
  update-patches                         ; Update patch changes caused by turtles, then determine new quantity and quality of patches
  update-patch-misc                      ; Update patch reporter/plot data and re-color patches

  ;;= World Updates
  update-N

  ;;= Plot procedure
  make-plots                             ; Update plots

  ;;= Increase timestep in simulation
  tick

  ;;= Stopping simulation conditions (Comment this out if using BehaviorSpace, and instead specify same condition)
;  if (count turtles = 0) and (initial-population > 0)[
;    stop
;  ]
;
; if any? turtles with [Xn < 0] [
;  stop
;]
;
;  if any? turtles with [Wn < 0] [
;  stop
;]

end

;;==============================================================================================================================================;;
;;=============================================== CALCULATE PARAMETERS FROM USER INTERFACE INPUT ===============================================;;
;;==============================================================================================================================================;;

;;========================== GENERAL PARAMETERS ==========================;;
;;========================================================================;;
to setup-parameters-from-interface
  ;;= Consumer-related general parameters
  set consumer-body-NC-ratio (1 / consumer-body-CN-ratio)  ; Convert carbon:nitrogen ratio given by user to N:C, for ease in future calculations (e.g., increased N:C ratio scales more intuitively with increased resource quality)
  set world-avg-resource-CN-ratio (target-world-avg-resource-CN-ratio)  ; Define the overall average carbon:nitrogen ratio of the entire landscape
  set world-avg-resource-NC-ratio (1 / world-avg-resource-CN-ratio)     ; Define the overall average nitrogen:carbon ratio of the entire landscape

  ;;= Resource cluster calculations
  ;; Determine what the quality of resources are in the increased quality resource patches
  set target-cluster-avg-resource-CN-ratio (target-min-CN-resource-ratio)  ; Set the intended resource quality for patches based on user input
  set target-cluster-avg-resource-NC-ratio (1 / target-cluster-avg-resource-CN-ratio)
  set target-noncluster-avg-resource-CN-ratio target-max-CN-resource-ratio
  set target-noncluster-avg-resource-NC-ratio (1 / target-noncluster-avg-resource-CN-ratio)

  ;;= Determine proportion of world-area that has no resources on it
  set fraction-world-resources (1 - (percent-world-no-resources / 100))

  ;;= Determine size of resource clusters from user input
  ;;=== Calculate the absolute maximum area that can be dedicated to resource clusters, based on the minimum, maximum, and average world nutrient ratios provided by user
  ifelse (num-resource-clusters = 0) [
    set max-cluster-resource-area 0
    set max-per-cluster-area 0
  ] [
    set max-cluster-resource-area floor ((fraction-world-resources * world-area * (world-avg-resource-CN-ratio - target-noncluster-avg-resource-CN-ratio)) / (target-cluster-avg-resource-CN-ratio - target-noncluster-avg-resource-CN-ratio))
    set max-per-cluster-area floor (max-cluster-resource-area / num-resource-clusters)
  ]

  (ifelse imbalanceddiet_response = "seek_higherquality_patch"[
    set movement_response 1
    set eatmore_response 0
  ]imbalanceddiet_response = "increase_ingestion" [
    set movement_response 0
    set eatmore_response 1
  ][
    print "Did you select a behavioral response to imbalanced nutrition?"
  ])
end


;;========================================================================================================================================================;;
;;=========================================================== SETUP/INITIALIZATION PROCEDURES: ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             PATCHES              ===========================================================;;
;;========================================================================================================================================================;;

;;========================== INITIALIZE RESOURCE PROCEDURES (Primary Production) ==========================;;
;;=========================================================================================================;;
to setup-patches

  ask patches [
  set resource-patch? 0 ; initialize all patches as not having resources
  ]

  ifelse (percent-world-no-resources > 0)[
    setup-resource-patches][
    ask patches [
      set resource-patch? 1 ; initialize all patches as having resources
    ]
  ]

  set resource-patches patches with [resource-patch? = 1]
  set no-resource-patches patches with [resource-patch? = 0]

  ;;= Resource Quality Calculations
  ;;=== Set up patches if there are no clusters, as set by user in interface
  ifelse num-resource-clusters = 0  [
    ask resource-patches [
      set resource-CN-ratio ( world-avg-resource-CN-ratio * (1 + (random-float (2 * overall-resource-patchiness) - overall-resource-patchiness)))
  ]
  ]

  ;;=== If there are clusters, perform the following procedures and commands to set up patches
  [
    setup-resource-cluster-center-coords
    setup-resource-clusters   ; Procedure to create clusters of high quality resources
    setup-resouce-nonclusters ; Procedure to create the baseline resource quality of noncluster patches
  ]

  ;;=== Set nitrogen:carbon ratios for future calculations
  ask resource-patches [
    set resource-NC-ratio (1 / resource-CN-ratio)
  ]

  ;;= Resource Quality Calculations
  ;;=== Set the initial quantity of resources per patch
  ask resource-patches[
    set resource-quantity (resource-carrying-capacity * initial-resource-quantity-as-percent-of-carrying-capacity) ; Sets the total quantity of resource to be what was defined in the user interface

  ; Set the intial amount of C and N per patch
    set resource-quantity-Cmol (resource-quantity)
    set resource-quantity-Nmol (resource-quantity-Cmol * resource-NC-ratio)
  ]

  ;;= Soil Nitrogen Calculations
  setup-soil

  ;;== Patch parameter update and aesthetics
  ;;=== Calculate the actual average world, cluster, and non-cluster patch resource quality (C:N and N:C) to ensure it matches the target
  ifelse (num-resource-clusters = 0) [
    set realized-avg-world-resource-CN-ratio mean [resource-CN-ratio] of resource-patches
    set realized-avg-world-resource-NC-ratio (1 / realized-avg-world-resource-CN-ratio)
  ]
  [
    set realized-avg-world-resource-CN-ratio mean [resource-CN-ratio] of resource-patches
    set realized-avg-world-resource-NC-ratio (1 / realized-avg-world-resource-CN-ratio)

    set realized-noncluster-avg-resource-CN-ratio mean [resource-CN-ratio] of resource-patches with [cluster-patch? = 0]
    set realized-noncluster-avg-resource-NC-ratio (1 / realized-noncluster-avg-resource-CN-ratio)

    set realized-cluster-avg-resource-CN-ratio mean [resource-CN-ratio] of resource-patches with [cluster-patch? = 1]
    set realized-cluster-avg-resource-NC-ratio (1 / realized-cluster-avg-resource-CN-ratio)
    ]

  ;;=== Color patches
  set min-cn-color-value (target-max-CN-resource-ratio + target-min-CN-resource-ratio)
  set max-cn-color-value (target-min-CN-resource-ratio - target-min-CN-resource-ratio)
  color-patches   ; Procedure that calculates pcolors
end

to color-patches
  ;;= Calculate the pcolor of patches
  ask resource-patches [
    set pcolor palette:scale-gradient [
  [140 127 90]
  [133 136 89]
  [125 145 83]
  [110 155 75]
  [90 130 65]
] resource-CN-ratio min-cn-color-value max-cn-color-value

    if resource-quantity <= 0 [
      set pcolor black
    ]
  ]

  ask no-resource-patches [
    set pcolor 8
  ]
end

;;========================== SETUP HIGH-QUALITY RESOURCE CLUSTERS ==========================;;
;;==========================================================================================;;
to setup-resource-patches
  ;;= Create temporary turtles that will walk around and "designate" which patches will be able to grow resources. It is possible to instead have random coordinates and radii to calculate patches of resources, but this creates more "round" shapes
  ;;==== Calculate how close together/homogenous resource patches will be
  let num-resource-turtles resource-patch-aggregation
  let max-resource-patches-per-turtle floor ((fraction-world-resources * world-area) / num-resource-turtles)

  ;;===== Create the temporary turtles
  create-resource-makers num-resource-turtles [
    setxy random-xcor random-ycor
    ask patch-here [
      set resource-patch? 1
    ]
  ]

  ;;===== Ask temporary turtles to walk around and label patches as resource-patches, then remove them from simulation
  ask resource-makers [
    repeat max-resource-patches-per-turtle [
      ;; Find the closest unmarked patch and update it
      ask min-one-of patches with [resource-patch? = 0] [distance myself] [
        set resource-patch? 1
      ]
      ;; Random walk
      rt random 90
      fd 1
    ]
    die  ; Remove the resource-maker turtle after completing its job
  ]

end

to setup-resource-cluster-center-coords
  ;;= Determine a list of coordinates for the starting patches of high-quality clusters of resources on resource-patches
  set resource-cluster-center-coords []  ; Initialize resource cluster coords list
  while [length resource-cluster-center-coords < 1] [
    let first-cluster-coords list random-pxcor random-pycor
    let first-cluster-coords-x first first-cluster-coords
    let first-cluster-coords-y last first-cluster-coords

    if [resource-patch?] of patch first-cluster-coords-x first-cluster-coords-y = 1 [
      set resource-cluster-center-coords lput first-cluster-coords resource-cluster-center-coords
    ]
  ]

  ;;=== Generate the rest of the coordinate sets (ensure each new pair is not within the radius of a previous coordinate pair)
  while [length resource-cluster-center-coords < num-resource-clusters] [
    let new-cluster-xcor random-pxcor
    let new-cluster-ycor random-pycor
    let new-cluster-coord list new-cluster-xcor new-cluster-ycor

    if [resource-patch?] of patch new-cluster-xcor new-cluster-ycor = 1 [
      set resource-cluster-center-coords lput new-cluster-coord resource-cluster-center-coords
    ]
  ]

  ;;=== Tag the patch for each center
  ask resource-patches [
    ; Start all patches with a binary assignment of not being in a cluster
    set cluster-patch? 0

  ;;=== Below, determine the center of clusters and create higher quality patches
    let current-coord list pxcor pycor ; Patch identifies its location

    if member? current-coord resource-cluster-center-coords [ ; If patch coordinate matches a coordinate of a cluster center, it will perform the follow commands to create clusters
      set cluster-patch? 1   ; Binary assignment for if a patch is part of a cluster at all. All centers automatically are in a cluster
      set cluster-center? 1
      set resource-CN-ratio target-cluster-avg-resource-CN-ratio
      set pcolor blue ; Coloring only for troubleshooting, can be removed
    ]
  ]
end

to setup-resource-clusters
  ;;= Procedure to create clusters of high-quality resources on resource patches

  set cluster-coords resource-cluster-center-coords
  set num-coords length cluster-coords            ; Number of cluster-makers to create
  let coord-index 0                               ; Initialize a counter for the list indices

  ;;=== Create temporary turtles to walk on patches and label patches as high-quality
  create-cluster-makers num-coords [
    let coord item coord-index cluster-coords     ; Assign each cluster-maker to a coordinate from the list
    setxy first coord last coord                  ; Set the cluster-maker's position based on the coordinates generated in the previous procedure
    set coord-index coord-index + 1              ; Increment the counter
  ]

  ;;=== Now have each cluster-maker perform the random walk and modify patches
  ask cluster-makers [
    repeat (max-per-cluster-area) [                                  ; Size of one cluster, normalized so that the patch the turtle is currently on won't overshoot the max cluster area
      ask min-one-of patches with [resource-patch? = 1 and cluster-patch? = 0 ] [ distance myself ] [
        set cluster-patch? 1
        set resource-CN-ratio (target-cluster-avg-resource-CN-ratio * (1 + (random-float (2 * overall-resource-patchiness) - overall-resource-patchiness)))
      ]
      rt random 360
      fd (0.5 + resource-cluster-scatter)
    ]
    die                                           ; Remove the cluster-maker after it finishes
  ]

  ;;=== It is possible that the clusters may not be exactly the resource cluster area. For simplicity and computational efficiency, this is corrected after the cluster formation:
  set total-resource-cluster-patches count resource-patches with [cluster-patch? = 1]                                                     ; Sum all patches in a resource cluster
  let discrepancy-resource-cluster-area (max-cluster-resource-area - total-resource-cluster-patches)                             ; Calculate the difference between the actual resource cluster patches and intended

  let resource-cluster-patches resource-patches with [cluster-patch? = 1]                                                                 ; Temp assignment to any cluster patches
  let resource-cluster-border-patches resource-cluster-patches with [count neighbors with [cluster-patch? = 1] < 8]              ; Find cluster patches with fewer than 8 neighbors that are also in a cluster
  let resource-cluster-edge-patches resource-cluster-border-patches with [count neighbors with [cluster-patch? = 0] >= 3]        ; Find the cluster patches with at least one neighbor that isn't in a cluster, and find patches that have at least 3 neighbors that are in a cluster to mimic morerealistic border

  let resource-noncluster-patches resource-patches with [cluster-patch? = 0]                                                              ; Temp assignment to any noncluster patches
  let resource-noncluster-border-patches resource-noncluster-patches with [count neighbors with [cluster-patch? = 1] > 0]        ; Finds noncluster patches with at least one neighbor that is in a cluster
  let resource-noncluster-edge-patches resource-noncluster-border-patches with [count neighbors with [cluster-patch? = 1] >= 3]  ; Find the cluster patches with at least one neighbor that isn't in a cluster, and find patches that have at least 3 neighbors that are in a cluster to mimic morerealistic border

  if total-resource-cluster-patches > max-cluster-resource-area [
    ask n-of (abs discrepancy-resource-cluster-area) resource-cluster-edge-patches [
      set resource-CN-ratio target-noncluster-avg-resource-CN-ratio
      set cluster-patch? 0
      set total-resource-cluster-patches count resource-patches with [cluster-patch? = 1]
    ]
  ]

  if total-resource-cluster-patches < max-cluster-resource-area[
    ask n-of (abs discrepancy-resource-cluster-area) resource-noncluster-edge-patches [
      set resource-CN-ratio target-cluster-avg-resource-CN-ratio
      set cluster-patch? 1
      set total-resource-cluster-patches count resource-patches with [cluster-patch? = 1]
    ]
  ]
end


;;========================== SETUP LOWER QUALITY NON-CLUSTER RESOURCES ==========================;;
;;===============================================================================================;;
to setup-resouce-nonclusters
 ;;= Set of the quality of patches that are lower-quality
  let noncluster-patches resource-patches with [cluster-patch? = 0]

  ;;=== Set the resource-NC-ratio to be the target noncluster average, with randomness
  ask noncluster-patches [
    set resource-CN-ratio (target-noncluster-avg-resource-CN-ratio * (1 + (random-float (2 * overall-resource-patchiness) -  overall-resource-patchiness)))
  ]
end


;;================================= SETUP SOILS RESOURCES =======================================;;
;;===============================================================================================;;
to setup-soil
  ;= Set the initial pool of N mol in soil, based on user interface
  ask resource-patches [
    set soil-N-quantity (resource-quantity-Nmol * patch-soil-N-multiplier)
  ]
end


;;========================================================================================================================================================;;
;;=========================================================== SETUP/INITIALIZATION PROCEDURES: ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             TURTLES              ===========================================================;;
;;========================================================================================================================================================;;

;;================================= INITIALIZE CONSUMER POPULATION ==================================;;
;;===================================================================================================;;
to setup-turtles
  ;;= Create consumers
  let turtle-cluster-center-x 0  ; X-coordinate of the cluster center
  let turtle-cluster-center-y 0  ; Y-coordinate of the cluster center
  let turtle-cluster-radius (max-pxcor * turtle-initial-location-spread) ; Radius within which turtles will be distributed
  create-turtles initial-population
  [
    ;;=== Randomly place turtles, with their spread across the landscape defined by the user interface
    let angle random-float 360                         ; Random angle in degrees
    let radius random-float turtle-cluster-radius      ; Random distance from the center
    let x turtle-cluster-center-x + radius * cos angle ; X-coordinate within the circle
    let y turtle-cluster-center-y + radius * sin angle ; Y-coordinate within the circle

    setxy x y
  ]

  ;;= Set consumers to be adults (others not set will remain juveniles)
  ; Determine the stage structure of population
  ask n-of (initial-population * initial-fraction-adults) turtles [
    set is-adult? 1
  ]

  ;;= Call procedures to setup detailed consumer states and parameters for DEB model
  ask turtles[
    initialize-DEB-parameters
    initialize-DEB-variables
  ]

  ;;= Set turtle aesthetics
  ask turtles[
    set shape "rabbit"
    ifelse (is-adult? = 0)[
      set color 127
    ][
      set color 106
    ]
  ]
end

;;=============================== SPECIFY INITIAL CONSUMER PARAMETERS ===============================;;
;;===================================================================================================;;
to initialize-DEB-parameters
  ;;= Stoichiometry parameters
  set q-F_C 1
  set q-F_N [resource-NC-ratio] of patch-here
  set q-Ec_C 1
  set q-Ec_N 0
  set q-En_C 1
  set q-En_N 0.33
  set q-V_C 1
  set q-V_N 0.2
  set q-R_C q-V_C
  set q-R_N q-V_N

  ;;= Dynamic Energy Budget parameters from Add-My-Pet portal
  ;;=== For Lepus timidus
  let p_AM_CrlBrack 12717.6
  let F_m 6.5
  set nu 0.028585
  let p_M 3044.77
  set kappa 0.96113
  set k_J 0.002
  set M_V_SqBrack 0.0125523
  let mu_E 550000
  set mdens_tot_max	 64.44363
  set k_M 0.3879724
  set J_XA_max_CrlBrack 0.03027993
  set M_H_baby 0.01218364

  ;;= Calculate other standard DEB parameters
  set F_h (J_XA_max_CrlBrack / F_m)
  set Food [resource-quantity] of patch-here
  set f (Food / (Food + F_h))  ; Type-II functional response

  ;;= Imposed multi-reserve DEB parameters
  ;;=== Reserve parameters
  let fraction_mdens_Ec 0.6
  let fraction_mdens_En (1 - fraction_mdens_Ec)
  set mdens_Ec_max (fraction_mdens_Ec * mdens_tot_max)
  set mdens_En_max (fraction_mdens_En * mdens_tot_max)
  set kappa_Ec kappa
  set kappa_En 0.8
  set kappaG_Ec 0.8
  set kappaG_En 0.9
  set kappaR_Ec 0.8
  set kappaR_En 0.9
  let J_EcA_M_CrlBrack (mdens_Ec_max * nu * M_V_SqBrack)

  ;;=== Reproduction parameters
  set M_H_p (10 * M_H_baby) ; Lower than DEB Add-My-Pet parameter
  let M_V_baby_AMPparm 0.0212256
  set M_V_baby 0.25 * M_V_baby_AMPparm  ; Lower than DEB Add-My-Pet parameter
  set M_R_bufferbaby 0.2779085 ; Calculation done for a pup to meet size, maturity, and reserve densities at birth
  set avg_litter_size 3
  set M_R_bufferlitter (M_R_bufferbaby * avg_litter_size)
  set scaled_mdens_Ec_baby 0.6
  set scaled_mdens_En_baby 0.9230769

  ;;=== Yield coefficients
  set Y_EcF (J_EcA_M_CrlBrack / J_XA_max_CrlBrack)
  set Y_EnF 0.85
  set Y_EnV ((q-V_N / q-V_C) * (q-En_C / q-En_N))
  set Y_EcV (1 - Y_EnV)
  set Y_EnR ((q-R_N / q-R_C) * (q-En_C / q-En_N))
  set Y_EcR (1 - Y_EnR)
  ;;=== Maintenance
  set sigma movement-somatic-maintenance-cost-multiplier

  ;;=== SU preference parameters
  set rho_S_N 0.01
  set rho_S_V 0.001
  set rho_D_N 0.01
  set rho_D_V 0.001

  ;;= Other important parameters to set
  set should-die? 0

  ;;= Set initial root solve guesses for numerical estimation
  set J_EcC_lastroot 0.1
  set J_EnC_lastroot 0.1
  set kS_lastroot 0
  set kD_lastroot 0
  set rdot_lastroot 0
  set ODE_Nconservation_test 0
end

;;============================ SPECIFY INITIAL CONSUMER STATE VARIABLES =============================;;
;;===================================================================================================;;
to initialize-DEB-variables
  ;;= Initialize turtle state variables before simulation begins.
  ;;=== Set up size, maturity, and reproduction buffer based on if juvenile or adult
  ifelse is-adult? = 0 [
    ; Set state variables for juveniles
    set M_V (M_V_baby * (1 + random-float 0.1))             ; Size (biomass)
    set M_H (M_H_baby + (random-float 1) * (M_H_p - M_H_baby))   ; Maturity
    set M_R 0                                                      ; Reproduction buffer
  ][
    set M_V (0.15 * (1 - random-float 0.1))                         ; Size (biomass)
    set M_H M_H_p                                                  ; Maturity
    set M_R random-float M_R_bufferlitter                          ; Reproduction buffer
  ]

  ;;=== Define other forms of somatic structure
  set V (M_V / M_V_SqBrack)                                        ; Size (volume)
  set L (V ^ (1 / 3))                                              ; Size (length)

  ;;=== Calculate reserves
  ;;===== Calculate scaled reserve densities
  let initial-mean-carbon-reserve-filled-fraction 0.7     ; Define starting value of carbon reserve density
  let initial-mean-nitrogen-reserve-filled-fraction 0.7   ; Define starting value of carbon reserve density
  let initial-error-carbon-reserve-filled-fraction 0.15   ; Define range of initial scaled reserve density
  let initial-error-nitrogen-reserve-filled-fraction 0.15 ; Define range of initial scaled reserve density
  set scaled_Mdens_Ec initial-mean-carbon-reserve-filled-fraction + random-float (2 * initial-error-carbon-reserve-filled-fraction) - initial-error-carbon-reserve-filled-fraction
  set scaled_Mdens_En initial-mean-nitrogen-reserve-filled-fraction + random-float (2 * initial-error-nitrogen-reserve-filled-fraction) - initial-error-nitrogen-reserve-filled-fraction

  ;;===== Calculate reserve densities
  if scaled_Mdens_Ec < 0.1 [
    set scaled_Mdens_Ec 0.1
  ]
  if scaled_Mdens_Ec > 1 [
    set scaled_Mdens_Ec 0.99
  ]

  if scaled_Mdens_En < 0.1 [
    set scaled_Mdens_En 0.1
  ]
  if scaled_Mdens_En > 1 [
    set scaled_Mdens_En 0.99
  ]

  set Mdens_Ec (scaled_Mdens_Ec * mdens_Ec_max)
  set Mdens_En (scaled_Mdens_En * mdens_En_max)

  ;;===== Calculate reserve biomass
  set M_Ec (Mdens_Ec * M_V)
  set M_En (Mdens_En * M_V)

  ;;=== Calculate pooled consumer excretion from all metabolic processes
  set Xc 0
  set Xn 0

  ;;=== Calculate consumer egestion
  set Wc 0
  set Wn 0

  ;;=== Calculate tick-specific consumer summaries
  update-consumerN
  set largest-Mv-to-date M_V       ; Largest somatic biomass achieved
end


;;========================================================================================================================================================;;
;;===========================================================   GO/RUN SIMULATION PROCEDURES:  ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             TURTLES              ===========================================================;;
;;========================================================================================================================================================;;

;;================ UPDATE TURTLE VALUES THAT MUST BE ZERO AT START OF TICK (Failsafe) ================;;
;;====================================================================================================;;
to new-tick-turtle-values
  ;;= Initialize important binary variables as 0 before rest of go procedures
  ask turtles[
    set should-die? 0
    set starved-to-death? 0
    set shrank-too-much? 0
    set demographic-stochasticity-death? 0
    set did-baby? 0
  ]
end

;;======================================== CONSUMER MOVEMENT ========================================;;
;;===================================================================================================;;
to move
;;= First, ask turtles to evaluate if their nitrogen-rich reserve is low (based on user-set threshold)
  ask turtles [
    ifelse movement_response = 1 [
    ;;=== If nitrogen-rich reserve is low, trigger movement to increase reserve
    ifelse scaled_Mdens_En <= N-reserve-threshold [
      ;;===== Seek high nitrogen food when body nitrogen (C:N) is above the threshold
      let best-patch min-one-of patches in-radius sight [resource-CN-ratio] ;; Look for best nitrogen patch nearby
      forward distance-consumer-walks
    ][
      if movement-type = "correlated-random-walk" [
        ;;=== Move in a correlated random walk
        rt random-normal 0 30
        forward distance-consumer-walks
      ]
      if movement-type = "brownian-walk" [
        rt random-float 360
        forward distance-consumer-walks
    ]
    ]
  ][
      if movement-type = "correlated-random-walk" [
        ;;=== Move in a correlated random walk
        rt random-normal 0 30
        forward distance-consumer-walks
      ]
      if movement-type = "brownian-walk" [
        rt random-float 360
        forward distance-consumer-walks
    ]
  ]
  ]
end

;;==================== CALCULATE HOW MUCH CONSUMERS ON EACH PATCH SHOULD CONSUME ====================;;
;;===================================================================================================;;
to hungry-hungry-hippos
  ask turtles [
    ;;=== Calculate the relevant simulation environemental information to send to R to compute the ODE
    set Food [resource-quantity] of patch-here
    set q-F_C 1
    set q-F_N [resource-NC-ratio] of patch-here

    set f (Food / (Food + F_h))  ; Type-II functional response

    ifelse eatmore_response = 1[
      set max-i-can-eat (f * (L ^ 2) * (J_XA_max_CrlBrack * (1 + ingestion_increase)))
    ][
      set max-i-can-eat (f * (L ^ 2) * J_XA_max_CrlBrack)
    ]

    ;;=== If there are more than a single turtle on this patch, calculate how much each can eat
    ifelse count turtles-here > 1 [
      let maximum-eating-potential sum [max-i-can-eat] of turtles-here
      ;;===== If there is less food than the turtles collectively plan to eat, create a system whereby the largest turtles get to eat first. If no food is left for smaller individuals, they will not eat
      if maximum-eating-potential > Food [
        set remaining-food Food
        foreach reverse sort-on [max-i-can-eat] turtles-here [
          the-turtle -> ask the-turtle [
            set what-i-will-eat min list [max-i-can-eat] of the-turtle remaining-food
            set remaining-food (remaining-food - what-i-will-eat)
          ]
        ]
      ]
    ][
     ;;=== If there is only one turtle on the patch, it can eat its full share
      set what-i-will-eat max-i-can-eat
    ]
  ]
end

;;=================== CALCULATE FLUXES AND ODE SOLUTIONS FOR INDIVIDUAL CONSUMERS ===================;;
;;===================================================================================================;;
to calculate-DEB-ODE
 ;;= Calculate change in consumer state variables this tick, performed in R using the SimpleR NetLogo Extension
  ask turtles [
    ;;=== Calculate the relevant simulation environemental information to send to R to compute the ODE
    set Food [resource-quantity] of patch-here
    if Food > 0 [
      set q-F_C 1
      set q-F_N [resource-NC-ratio] of patch-here
    ]

    let J_S_basal (M_V * k_M)
    set J_S (1 + (distance-consumer-walks * sigma)) * (M_V * k_M)
    set J_D (M_H * k_J)

    ;;=== Create a list of consumer state variables, patch variables, and parameters that R will need to calculate ODE
    let netlogo_input []  ; Create empty list
    ;;===== Send state variable information
    set netlogo_input lput what-i-will-eat netlogo_input
    set netlogo_input lput M_V netlogo_input
    set netlogo_input lput V netlogo_input
    set netlogo_input lput M_Ec netlogo_input
    set netlogo_input lput M_En netlogo_input
    set netlogo_input lput M_H netlogo_input
    set netlogo_input lput M_R netlogo_input
    set netlogo_input lput Wc netlogo_input
    set netlogo_input lput Wn netlogo_input
    set netlogo_input lput Xc netlogo_input
    set netlogo_input lput Xn netlogo_input
    ;;===== Send stoichiometry information
    set netlogo_input lput q-F_C netlogo_input
    set netlogo_input lput q-V_C netlogo_input
    set netlogo_input lput q-R_C netlogo_input
    set netlogo_input lput q-En_C netlogo_input
    set netlogo_input lput q-F_N netlogo_input
    set netlogo_input lput q-V_N netlogo_input
    set netlogo_input lput q-R_N netlogo_input
    set netlogo_input lput q-En_N netlogo_input
    ;;===== Send maintenance costs
    set netlogo_input lput J_S netlogo_input
    set netlogo_input lput J_D netlogo_input
    ;;===== Send relevant parameter values
    set netlogo_input lput kappa_Ec netlogo_input
    set netlogo_input lput kappa_En netlogo_input
    set netlogo_input lput kappaG_Ec netlogo_input
    set netlogo_input lput kappaG_En netlogo_input
    set netlogo_input lput kappaR_Ec netlogo_input
    set netlogo_input lput kappaR_En netlogo_input
    set netlogo_input lput F_h netlogo_input
    set netlogo_input lput J_XA_max_CrlBrack netlogo_input
    set netlogo_input lput nu netlogo_input
    set netlogo_input lput mdens_Ec_max netlogo_input
    set netlogo_input lput mdens_En_max netlogo_input
    set netlogo_input lput Y_EcF netlogo_input
    set netlogo_input lput Y_EnF netlogo_input
    set netlogo_input lput rho_S_N netlogo_input
    set netlogo_input lput rho_S_V netlogo_input
    set netlogo_input lput Y_EcV netlogo_input
    set netlogo_input lput Y_EnV netlogo_input
    set netlogo_input lput rho_D_N netlogo_input
    set netlogo_input lput rho_D_V netlogo_input
    set netlogo_input lput Y_EcR netlogo_input
    set netlogo_input lput Y_EnR netlogo_input
    set netlogo_input lput M_H_p netlogo_input
    ;;===== Send the previous numerical root solutions
    set netlogo_input lput J_EcC_lastroot netlogo_input
    set netlogo_input lput J_EnC_lastroot netlogo_input
    set netlogo_input lput kS_lastroot netlogo_input
    set netlogo_input lput kD_lastroot netlogo_input
    set netlogo_input lput rdot_lastroot netlogo_input

    ;;=== Send list to R
    (sr:set "netlogo_input" netlogo_input)

    ;;=== Perform ODE calculations with function specified in R file
    sr:run "x <- DEB_ODE_calc(netlogo_input)"

    ;;=== Set ODE results from information sent back to NetLogo by R
    set dFood-dt sr:runresult "x[1]"
    set dM_V-dt sr:runresult "x[2]"
    set dM_Ec-dt sr:runresult "x[3]"
    set dM_En-dt sr:runresult "x[4]"
    set dM_H-dt sr:runresult "x[5]"
    set dM_R-dt sr:runresult "x[6]"
    set dXc-dt sr:runresult "x[7]"
    set dXn-dt sr:runresult "x[8]"
    set dWc-dt sr:runresult "x[9]"
    set dWn-dt sr:runresult "x[10]"
    set dWn-dt sr:runresult "x[10]"

    ;;=== Update numerical root solutions returned from R
    set J_EcC_root sr:runresult "x[11]"
    set J_EnC_root sr:runresult "x[12]"
    set kS_root sr:runresult "x[13]"
    set kD_root sr:runresult "x[14]"
    set rdot_root sr:runresult "x[15]"

    set ODE_Nconservation_test sr:runresult "x[16]"

    ;;=== Calculate the total amount of food eaten this tick, for convenience in later patch calculations
    set food_ingested (-1 * dFood-dt)
  ]
end

;;===================== UPDATE CONSUMER STATE VARIABLES AND RELEVANT PARAMETERS =====================;;
;;===================================================================================================;;
to update-turtle-variables
  ;;= Update DEB model variables
  ask turtles[
    set M_V (M_V + dM_V-dt)
    set M_H (M_H + dM_H-dt)
    set M_R (M_R + dM_R-dt)
    set M_Ec (M_Ec + dM_Ec-dt)
    set M_En (M_En + dM_En-dt)

    set Xc (Xc + dXc-dt)
    set Xn (Xn + dXn-dt)
    set Wc (Wc + dWc-dt)
    set Wn (Wn + dWn-dt)

    set V (M_V / M_V_SqBrack)
    set L (V ^ (1 / 3))
    set mdens_Ec (M_Ec / M_V)
    set mdens_En (M_En / M_V)

    set scaled_mdens_Ec (mdens_Ec / mdens_Ec_max)
    set scaled_mdens_En (mdens_En / mdens_En_max)

    if M_V > largest-Mv-to-date [
      set largest-Mv-to-date M_V
    ]

    ifelse M_H >= M_H_p [
      set is-adult? 1
      set color 106
    ][
      set is-adult? 0
      set color 127
    ]

    ;;= Update the roots of the numerical solving procedure, for the start of the next tick
    set J_EcC_lastroot J_EcC_root
    set J_EnC_lastroot J_EnC_root
    set kS_lastroot kS_root
    set kD_lastroot kD_root
    set rdot_lastroot rdot_root

    update-consumerN
  ]
end

;;====================================== REPRODUCTION PROCEDURE =====================================;;
;;===================================================================================================;;
to have-baby?
  ;;= For turtles whose reproduction buffer has accumulated enough biomass, allow them to reproduce, and re-calculate their reproduction buffer
  ask turtles [
    if M_R >= M_R_bufferlitter [
      ;;=== Initialize newborns' parameters and states
      hatch avg_litter_size [

        ;;===== First, clear the parent variables
        reset-ODE-and-fluxes

        ;;===== Set important parameters for juveniles
        set is-adult? 0
        set did-baby? 0
        set starved-to-death? 0
        set demographic-stochasticity-death? 0
        set should-die? 0

        initialize-DEB-parameters

        ;;===== Set state variables
        set M_V M_V_baby
        set V (M_V / M_V_SqBrack)
        set L (V ^ (1 / 3))
        set scaled_mdens_Ec scaled_mdens_Ec_baby
        set scaled_mdens_En scaled_mdens_En_baby
        set mdens_Ec (scaled_mdens_Ec * mdens_Ec_max)
        set mdens_En (scaled_mdens_En * mdens_En_max)
        set M_Ec (mdens_Ec * M_V)
        set M_En (mdens_En * M_V)
        set M_H M_H_baby
        set M_R 0
        set Xc 0
        set Xn 0
        set Wc 0
        set Wn 0
        set largest-Mv-to-date M_V

        update-consumerN

        set color 127
      ]
      ;;=== For new parent, reduce their reproduction buffer, indicate they reproduced this turn, and re-calculate their N composition
      set M_R (M_R - M_R_bufferlitter)
      set did-baby? 1
      write "Babies! :)"

      update-consumerN
    ]
  ]
end

;;================================ MORTALITY CALCULATIONS PROCEDURES ================================;;
;;===================================================================================================;;
to have-die?
  ;;= Define death conditions of turtles each tick
  ask turtles [
   ;;=== Starvation scenarios
    if (scaled_Mdens_Ec <= 0.01) and (scaled_Mdens_En <= 0.01) [
      set starved-to-death? 1
      set should-die? 1
    ]

    if (M_V <= 0) [
      set shrank-too-much? 1
      set should-die? 1
    ]

    if (1 - (M_V / largest-Mv-to-date)) > fraction-structure-shrinkage-allowed[
      set shrank-too-much? 1
      set should-die? 1
    ]

   ;;=== Demographic stochasiticity
    if random-float 1 < mortality-risk [
      set demographic-stochasticity-death? 1
      set should-die? 1
    ]
  ]
end

;;========================== NITROGEN DEPOSITION BY CONSUMERS ONTO PATCHES ==========================;;
;;===================================================================================================;;
to patch-N-added-by-turtles
  ;;= Reflect the amount of N available for resource growth that was left by consumer during tick procedures, calculated during patch-N-added-by-turtles
  ask turtles [
    ;;=== The amount of nitrogen (in moles) that is added to the environment is composed of its excretion (X_n) and egestion (W_n) calculated that timestep, plus its corpse if it is to die
    ifelse should-die? = 1 [
      set corpse-n consumer-total-Nmol
      set turtle-input-N-this-tick (Wn + Xn + corpse-N)
    ]
    [
      set turtle-input-N-this-tick (Wn + Xn)
    ]
    ;;=== Tell patch how much N was left on it by turtle
    ask patch-here [
      set turtle-alterations-patch-pooled-N (turtle-alterations-patch-pooled-N + ([turtle-input-N-this-tick] of myself))
    ]
  ]
end


;;=================================== DEATH/REMOVAL OF CONSUMERS ====================================;;
;;===================================================================================================;;
to go-die
 ;;= Remove turtles assigned to die from the simulation
  ask turtles[
    if should-die? = 1 [
    write "RIP </3"
    die
    ]
  ]
end

;;======================= RESET RELEVANT TURTLE STATE VARIABLES AND PARAMETERS ======================;;
;;===================================================================================================;;
to reset-ODE-and-fluxes
;;= Reset relevant turtle variables
  ;;=== Reset all ODEs to 0 for the next tick
  set dFood-dt 0
  set food_ingested 0
  set dWc-dt 0
  set dWn-dt 0
  set dM_Ec-dt 0
  set dM_En-dt 0
  set dMdens_Ec-dt 0
  set dMdens_En-dt 0
  set dM_V-dt 0
  set dM_H-dt 0
  set dM_R-dt 0
  set dXc-dt 0
  set dXn-dt 0

  ;;=== Reset environmental output to 0 again
  set Xc 0
  set Xn 0
  set Wc 0
  set Wn 0
  set turtle-input-N-this-tick 0

  set ODE_Nconservation_test 0
end

to update-consumerN
   set consumer-total-Nmol (((q-V_N / q-V_C) * M_V) + ((q-En_N / q-En_C) * M_En) + ((q-R_N / q-R_C) * (M_H + M_R)))
end

;;========================================================================================================================================================;;
;;===========================================================   GO/RUN SIMULATION PROCEDURES:  ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             PATCHES              ===========================================================;;
;;========================================================================================================================================================;;

;;========================= UPDATE PATCH PARAMETERS AFTER TURTLE INTERACTION ========================;;
;;===================================================================================================;;
to update-patches
  ;;=== Update resource quantity and the amount of C mol and N mol in resources, after turtles have eaten and deposited N
  ask patches [
    ; Remove eaten food and recalculate the C and N moles of remaining resources
    let food-eaten sum [food_ingested] of turtles-here
    set resource-quantity (resource-quantity - food-eaten)
    set resource-quantity-Cmol resource-quantity
    set resource-quantity-Nmol (resource-quantity-Cmol * resource-NC-ratio)

    ; Update the N in soil from consumer input (all patches have soil, regardless of if resources are there or not)
    set soil-N-quantity (soil-N-quantity + turtle-alterations-patch-pooled-N)
  ]

  ;;=== Update resource quality (through primary production) and quality (through carbon dilution and/or soil nitrogen uptake)
  ask resource-patches[
    ; Calculate how many resources would grow if they could infinitely grow
    let resource-growth-infinite ((resource-specific-growth-rate / 100) * resource-carrying-capacity)
    ; Correct for the resource carrying capacity
    ifelse (resource-growth-infinite + resource-quantity) > resource-carrying-capacity [
      set resource-growth-max (resource-carrying-capacity - resource-quantity)
    ][
      set resource-growth-max resource-growth-infinite
    ]

    ;;===== Test: How much C will be added if maximum growth occurs
    let new-resource-quantity-Cmol-max resource-growth-max
    let max-total-Cmol (new-resource-quantity-Cmol-max + resource-quantity-Cmol)

    ;;===== Test: How much N will be available for growth and quality improvement
    let new-resource-quantity-Nmol-max (soil-N-quantity * (percent-resource-quality-improve / 100))
    let max-total-Nmol (new-resource-quantity-Nmol-max + resource-quantity-Nmol)

    ;;===== Test: Calculate the new projected C:N resource ratio
    ifelse max-total-Nmol > 0 [
      set resource-CN-ratio-test (max-total-Cmol / max-total-Nmol)
    ][
      set resource-CN-ratio-test 0
    ]

    ;;===== Correct for edge cases where the new CN ratio is too much or too little
    (ifelse (resource-CN-ratio-test > target-max-CN-resource-ratio) [
      set new-resource-quantity-Cmol ((target-max-CN-resource-ratio * max-total-Nmol) - resource-quantity-Cmol)
      set resource-growth-max new-resource-quantity-Cmol
      set new-resource-quantity-Nmol new-resource-quantity-Nmol-max
    ]
    (resource-CN-ratio-test < target-min-CN-resource-ratio)[
      set new-resource-quantity-Nmol (((1 / target-min-CN-resource-ratio) * max-total-Cmol) - resource-quantity-Nmol)
      set new-resource-quantity-Cmol new-resource-quantity-Cmol-max
    ][
      set new-resource-quantity-Cmol new-resource-quantity-Cmol-max
      set new-resource-quantity-Nmol new-resource-quantity-Nmol-max
    ])

    ;;===== Finally, add new growth, Cmol, and Nmol to resources
    set resource-quantity (resource-quantity + resource-growth-max)
    set resource-quantity-Cmol (resource-quantity-Cmol + new-resource-quantity-Cmol)
    set resource-quantity-Nmol (resource-quantity-Nmol + new-resource-quantity-Nmol)

    ;;=== Remove N from soil that was used for resources
    set soil-N-quantity (soil-N-quantity - new-resource-quantity-Nmol)

    ;;=== Re-calculate resource quality ratios
    ifelse resource-quantity > 0 and resource-quantity-Nmol > 0 [
      set resource-CN-ratio (resource-quantity-Cmol / resource-quantity-Nmol)
      set resource-NC-ratio (1 / resource-CN-ratio)
      ]
    [
      set resource-CN-ratio 0
      set resource-NC-ratio 0
    ]
  ]
end

;;======================= UPDATE PATCH REPORTER VALUES AND COLORING/AESTHETICS ======================;;
;;===================================================================================================;;
to update-patch-misc
  ;;= Update informational monitors
  ifelse (num-resource-clusters = 0) [
    set realized-avg-world-resource-CN-ratio mean [resource-CN-ratio] of resource-patches with [resource-quantity > 0]
    set realized-avg-world-resource-NC-ratio (1 / realized-avg-world-resource-CN-ratio)
  ]
  [
    ifelse (any? patches with [resource-quantity > 0])[
      set realized-avg-world-resource-CN-ratio mean [resource-CN-ratio] of resource-patches with [resource-quantity > 0]
      set realized-avg-world-resource-NC-ratio (1 / realized-avg-world-resource-CN-ratio)
    ][
      set realized-avg-world-resource-CN-ratio 0
      set realized-avg-world-resource-NC-ratio 0
    ]
  ]

  ;;= Re-color patches after tick
  color-patches

  ;;= Reset the turtle input of N this tick
    ;;= Ask patches to clear the amount of N mol deposited by consumers
  ask patches [
    set turtle-alterations-patch-pooled-N 0
  ]
end

;;========================================================================================================================================================;;
;;===========================================================   GO/RUN SIMULATION PROCEDURES:  ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================              WORLD               ===========================================================;;
;;========================================================================================================================================================;;

to update-N
  set sum-resource-N-total sum [resource-quantity-Nmol] of patches
  set sum-consumer-N-total sum [consumer-total-Nmol] of turtles
  set sum-consumer-N-total_adult sum [consumer-total-Nmol] of turtles with [is-adult? = 1]
  set sum-consumer-N-total_juv sum [consumer-total-Nmol] of turtles with [is-adult? = 0]
  set sum-soil-N-total sum [soil-N-quantity] of patches
  set sum-world-N-total (sum-resource-N-total + sum-consumer-N-total + sum-soil-N-total)
end

;;========================================================================================================================================================;;
;;===========================================================   GO/RUN SIMULATION PROCEDURES:  ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================              PLOTS               ===========================================================;;
;;========================================================================================================================================================;;

to make-plots
  plot-total-N
  plot-resource-densities
  plot-consumer-densities
  plot-consumer-reserve
  plot-consumer-size
  plot-consumer-repro
  plot-resource-mean-quality
  plot-consumer-Nmol-waste
  plot-consumer-Nmol-excretion
  plot-worldN-alone
end

to plot-total-N
  ;;= Calculate the sum of the patch variable
  set-current-plot "Total World N Dynamics"
end

to plot-resource-densities
  ;= Plot total resources in simulation
  set total-resource-quantity sum [resource-quantity] of resource-patches
  set-current-plot "Total Resources"
end

to plot-consumer-densities
  ;= Plot consumer abundances in simulation
  set total-consumers count turtles
  set total-adult-consumers count turtles with [is-adult? = 1]
  set total-juvenille-consumers count turtles with [is-adult? = 0]
  set-current-plot "Total Consumers"
end

to plot-consumer-reserve
  ;= Plot the average consumer scaled reserve densities
  ifelse count turtles > 0 [
    set avg-Creserve-scaleddensity mean [scaled_Mdens_Ec] of turtles
    set avg-Nreserve-scaleddensity mean [scaled_Mdens_En] of turtles
  ][
    set avg-Creserve-scaleddensity 0
    set avg-Nreserve-scaleddensity 0
  ]
  set-current-plot "Consumer Average Scaled Reserve Densities"
end

to plot-consumer-size
  ;= Plot the average consumer size
  ifelse count turtles > 0 [
    set avg-M_V mean [M_V] of turtles
  ][
    set avg-M_V 0
  ]
  set-current-plot "Consumer Avg Size"
end

to plot-consumer-repro
  ;= Plot the average consumer reproduction buffer of adults and maturity level of juveniles
  ifelse count turtles > 0 [
    ifelse count turtles with [is-adult? = 1] > 0 [
      set avg-M_R mean [M_R] of turtles with [is-adult? = 1]
    ][
      set avg-M_R 0
    ]
    ifelse count turtles with [is-adult? = 0] > 0 [
      set avg-M_H mean [M_H] of turtles with [is-adult? = 0]
    ][
      set avg-M_H 0
    ]
  ][
    set avg-M_R 0
    set avg-M_H 0
  ]
  set-current-plot "Consumer Avg Repro and Maturity"
end

to plot-resource-mean-quality
  set-current-plot "Avg World C:N Ratio"
end

to plot-consumer-Nmol-waste
  set-current-plot "Total Consumer Waste N"
end

to plot-consumer-Nmol-excretion
  set-current-plot "Total Consumer Excreted N"
end

to plot-worldN-alone
  set Nmol_ODE_check_sum sum[ODE_Nconservation_test] of turtles
  set-current-plot "Nmol Conservation ODE Check"
end
@#$#@#$#@
GRAPHICS-WINDOW
5
15
451
462
-1
-1
10.7
1
10
1
1
1
0
1
1
1
-20
20
-20
20
1
1
1
ticks
30.0

SLIDER
841
191
1048
224
initial-population
initial-population
0
100
14.0
1
1
NIL
HORIZONTAL

BUTTON
591
58
655
92
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
885
35
994
95
consumer-body-CN-ratio
5.0
1
0
Number

SLIDER
841
263
1048
296
turtle-initial-location-spread
turtle-initial-location-spread
0
1
1.0
0.1
1
NIL
HORIZONTAL

BUTTON
657
58
721
92
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
577
328
829
361
overall-resource-patchiness
overall-resource-patchiness
0
0.5
0.1
0.05
1
NIL
HORIZONTAL

SLIDER
840
449
1051
482
N-reserve-threshold
N-reserve-threshold
0
0.75
0.5
0.05
1
NIL
HORIZONTAL

CHOOSER
840
402
1050
447
movement-type
movement-type
"brownian-walk" "correlated-random-walk"
1

SLIDER
840
485
1051
518
movement-somatic-maintenance-cost-multiplier
movement-somatic-maintenance-cost-multiplier
0
1
0.0
0.1
1
NIL
HORIZONTAL

SLIDER
841
227
1049
260
initial-fraction-adults
initial-fraction-adults
0
1
0.6
0.1
1
NIL
HORIZONTAL

SLIDER
577
258
832
291
resource-specific-growth-rate
resource-specific-growth-rate
0
5
1.0
0.5
1
NIL
HORIZONTAL

SLIDER
577
224
831
257
resource-carrying-capacity
resource-carrying-capacity
0
5
1.0
0.25
1
NIL
HORIZONTAL

SLIDER
577
190
832
223
initial-resource-quantity-as-percent-of-carrying-capacity
initial-resource-quantity-as-percent-of-carrying-capacity
0.1
1
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
841
298
1048
331
mortality-risk
mortality-risk
0.0001
0.1
6.0E-4
0.0005
1
NIL
HORIZONTAL

PLOT
1087
43
1458
163
Total World N Dynamics
NIL
NIL
0.0
10.0
0.0
0.05
true
true
"" ""
PENS
"World Total N" 1.0 0 -16777216 true "" "plot sum-world-N-total"

INPUTBOX
870
100
997
160
target-min-CN-resource-ratio
30.0
1
0
Number

INPUTBOX
739
100
869
160
target-max-CN-resource-ratio
60.0
1
0
Number

INPUTBOX
588
100
739
160
target-world-avg-resource-CN-ratio
45.0
1
0
Number

SLIDER
575
470
826
503
patch-soil-N-multiplier
patch-soil-N-multiplier
0
1
0.25
0.25
1
NIL
HORIZONTAL

SLIDER
577
293
830
326
num-resource-clusters
num-resource-clusters
0
30
11.0
1
1
NIL
HORIZONTAL

MONITOR
1219
646
1459
695
% Resources Initially High-Quality
count patches with [cluster-patch? = 1] / (fraction-world-resources * world-area)
2
1
12

SLIDER
577
363
828
396
resource-cluster-scatter
resource-cluster-scatter
0
10
1.0
0.5
1
NIL
HORIZONTAL

PLOT
1458
43
1829
171
Total Consumers
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -8630108 true "" "plot total-consumers"
"pen-1" 1.0 0 -13345367 true "" "plot total-adult-consumers"
"pen-2" 1.0 0 -4699768 true "" "plot total-juvenille-consumers"

PLOT
1087
528
1272
648
Total Resources
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -14439633 true "" "plot total-resource-quantity"

TEXTBOX
658
10
953
45
User Input/Simulation Controls
18
0.0
1

TEXTBOX
648
165
767
195
Patch Controls
16
0.0
1

TEXTBOX
878
167
1028
187
Consumer Controls
16
0.0
1

PLOT
1458
170
1829
290
Consumer Average Scaled Reserve Densities
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -8431303 true "" "plot avg-Creserve-scaleddensity"
"pen-1" 1.0 0 -2674135 true "" "plot avg-Nreserve-scaleddensity"

SLIDER
574
506
826
539
percent-resource-quality-improve
percent-resource-quality-improve
0
10
5.0
0.5
1
NIL
HORIZONTAL

PLOT
1271
528
1458
648
Avg World C:N Ratio
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [resource-CN-ratio] of resource-patches with [resource-quantity > 0]"

MONITOR
1087
646
1220
695
Avg Food Ingested
mean [food_ingested] of turtles
7
1
12

SLIDER
577
398
827
431
percent-world-no-resources
percent-world-no-resources
0
99
10.0
5
1
NIL
HORIZONTAL

SLIDER
575
434
827
467
resource-patch-aggregation
resource-patch-aggregation
1
10
10.0
1
1
NIL
HORIZONTAL

PLOT
1458
420
1641
544
Total Consumer Waste N
NIL
NIL
0.0
10.0
0.0
0.01
true
false
"" ""
PENS
"pen-2" 1.0 0 -955883 true "" "plot sum [Wn] of turtles"

PLOT
1643
420
1829
544
Total Consumer Excreted N
NIL
NIL
0.0
10.0
0.0
0.01
true
false
"" ""
PENS
"default" 1.0 0 -13791810 true "" "plot sum [Xn] of turtles"

SLIDER
838
332
1049
365
distance-consumer-walks
distance-consumer-walks
0
3
1.0
0.5
1
NIL
HORIZONTAL

SLIDER
840
368
1050
401
fraction-structure-shrinkage-allowed
fraction-structure-shrinkage-allowed
0
0.8
0.5
0.1
1
NIL
HORIZONTAL

PLOT
1458
290
1652
421
Consumer Avg Size
NIL
NIL
0.0
10.0
0.0
0.1
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot avg-M_V"

PLOT
1652
290
1829
420
Consumer Avg Repro and Maturity
NIL
NIL
0.0
10.0
0.0
0.8
true
false
"" ""
PENS
"default" 1.0 0 -8630108 true "" "plot avg-M_R"
"pen-1" 1.0 0 -7858858 true "" "plot avg-M_H"

TEXTBOX
1382
12
1532
34
Data Output/Plots
18
0.0
1

PLOT
1458
546
1831
696
Nmol Conservation ODE Check
NIL
NIL
0.0
10.0
0.0
1.0E-14
true
true
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot Nmol_ODE_check_sum"

PLOT
1087
164
1458
284
Total Consumer N Dynamics
NIL
NIL
0.0
10.0
0.0
0.1
true
false
"" ""
PENS
"default" 1.0 0 -8630108 true "" "plot sum-consumer-N-total"
"pen-1" 1.0 0 -13345367 true "" "plot sum-consumer-N-total_adult"
"pen-2" 1.0 0 -5825686 true "" "plot sum-consumer-N-total_juv"

PLOT
1088
285
1459
405
Total Soil N Dynamics
NIL
NIL
0.0
10.0
0.0
0.1
true
false
"" ""
PENS
"default" 1.0 0 -8431303 true "" "plot sum-soil-N-total"

PLOT
1087
408
1455
528
Total Resource N Dynamics
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -14439633 true "" "plot sum-resource-N-total"

SLIDER
839
520
1050
553
sight
sight
1
15
10.0
1
1
NIL
HORIZONTAL

CHOOSER
361
482
518
527
imbalanceddiet_response
imbalanceddiet_response
"seek_higherquality_patch" "increase_ingestion"
1

SLIDER
838
554
1051
587
ingestion_increase
ingestion_increase
0
5
3.0
1
1
NIL
HORIZONTAL

@#$#@#$#@
## WHAT IS IT?

The model creates a landscape of resources of varying quality that agents can ingest. The resource quality is set up in a variety of quality patchiness such that scenario modeling can be done for the "strength" of spatial patchiness.

## HOW IT WORKS

Agents operate on a 

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

rabbit
false
0
Polygon -7500403 true true 61 150 76 180 91 195 103 214 91 240 76 255 61 270 76 270 106 255 132 209 151 210 181 210 211 240 196 255 181 255 166 247 151 255 166 270 211 270 241 255 240 210 270 225 285 165 256 135 226 105 166 90 91 105
Polygon -7500403 true true 75 164 94 104 70 82 45 89 19 104 4 149 19 164 37 162 59 153
Polygon -7500403 true true 64 98 96 87 138 26 130 15 97 36 54 86
Polygon -7500403 true true 49 89 57 47 78 4 89 20 70 88
Circle -16777216 true false 37 103 16
Line -16777216 false 44 150 104 150
Line -16777216 false 39 158 84 175
Line -16777216 false 29 159 57 195
Polygon -5825686 true false 0 150 15 165 15 150
Polygon -5825686 true false 76 90 97 47 130 32
Line -16777216 false 180 210 165 180
Line -16777216 false 165 180 180 165
Line -16777216 false 180 165 225 165
Line -16777216 false 180 210 210 240

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="exp_redo_movement" repetitions="5" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1000"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles with [is-adult? = 0]</metric>
    <metric>count turtles with [is-adult? = 1]</metric>
    <metric>mean [M_V] of turtles</metric>
    <metric>mean [M_Ec] of turtles</metric>
    <metric>variance [M_Ec] of turtles</metric>
    <metric>mean [M_En] of turtles</metric>
    <metric>variance [M_En] of turtles</metric>
    <metric>mean [M_H] of turtles</metric>
    <metric>mean [M_R] of turtles</metric>
    <metric>mean [Wc] of turtles</metric>
    <metric>variance [Wc] of turtles</metric>
    <metric>mean [Wn] of turtles</metric>
    <metric>variance [Wn] of turtles</metric>
    <metric>mean [Xc] of turtles</metric>
    <metric>variance [Xc] of turtles</metric>
    <metric>mean [Xn] of turtles</metric>
    <metric>variance [Xn] of turtles</metric>
    <metric>mean [resource-quantity-Nmol] of patches with [resource-patch? = 1]</metric>
    <metric>variance [resource-quantity-Nmol] of patches with [resource-patch? = 1]</metric>
    <metric>mean [resource-CN-ratio] of patches</metric>
    <metric>variance [resource-CN-ratio] of patches with [resource-patch? = 1]</metric>
    <metric>mean [resource-quantity] of patches</metric>
    <metric>variance [resource-quantity] of patches</metric>
    <metric>sum [resource-quantity] of patches</metric>
    <metric>sum [resource-quantity-Nmol] of patches</metric>
    <metric>sum [soil-N-quantity] of patches with [resource-patch? = 0]</metric>
    <metric>sum [soil-N-quantity] of patches with [resource-patch? = 1]</metric>
    <metric>sum [consumer-total-Nmol] of turtles</metric>
    <metric>sum-world-N-total</metric>
    <enumeratedValueSet variable="imbalanceddiet_response">
      <value value="&quot;seek_higherquality_patch&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="target-max-CN-resource-ratio">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cluster-scatter">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patch-soil-N-multiplier">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-patch-aggregation">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-specific-growth-rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-type">
      <value value="&quot;correlated-random-walk&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-resource-quality-improve">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="turtle-initial-location-spread">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-risk">
      <value value="6.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sight">
      <value value="10"/>
    </enumeratedValueSet>
    <steppedValueSet variable="movement-somatic-maintenance-cost-multiplier" first="0" step="0.25" last="1"/>
    <enumeratedValueSet variable="percent-world-no-resources">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-carrying-capacity">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consumer-body-CN-ratio">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="target-world-avg-resource-CN-ratio">
      <value value="45"/>
    </enumeratedValueSet>
    <steppedValueSet variable="N-reserve-threshold" first="0" step="0.25" last="0.75"/>
    <steppedValueSet variable="num-resource-clusters" first="0" step="10" last="10"/>
    <enumeratedValueSet variable="distance-consumer-walks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-structure-shrinkage-allowed">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-fraction-adults">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-resource-quantity-as-percent-of-carrying-capacity">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overall-resource-patchiness">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="target-min-CN-resource-ratio">
      <value value="30"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="exp_redo_compensatoryingest" repetitions="5" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <exitCondition>not any? turtles</exitCondition>
    <metric>count turtles with [is-adult? = 0]</metric>
    <metric>count turtles with [is-adult? = 1]</metric>
    <metric>mean [M_V] of turtles</metric>
    <metric>mean [M_Ec] of turtles</metric>
    <metric>variance [M_Ec] of turtles</metric>
    <metric>mean [M_En] of turtles</metric>
    <metric>variance [M_En] of turtles</metric>
    <metric>mean [M_H] of turtles</metric>
    <metric>mean [M_R] of turtles</metric>
    <metric>mean [Wc] of turtles</metric>
    <metric>variance [Wc] of turtles</metric>
    <metric>mean [Wn] of turtles</metric>
    <metric>variance [Wn] of turtles</metric>
    <metric>mean [Xc] of turtles</metric>
    <metric>variance [Xc] of turtles</metric>
    <metric>mean [Xn] of turtles</metric>
    <metric>variance [Xn] of turtles</metric>
    <metric>mean [resource-quantity-Nmol] of patches with [resource-patch? = 1]</metric>
    <metric>variance [resource-quantity-Nmol] of patches with [resource-patch? = 1]</metric>
    <metric>mean [resource-CN-ratio] of patches</metric>
    <metric>variance [resource-CN-ratio] of patches with [resource-patch? = 1]</metric>
    <metric>mean [resource-quantity] of patches</metric>
    <metric>variance [resource-quantity] of patches</metric>
    <metric>sum [resource-quantity] of patches</metric>
    <metric>sum [resource-quantity-Nmol] of patches</metric>
    <metric>sum [soil-N-quantity] of patches with [resource-patch? = 0]</metric>
    <metric>sum [soil-N-quantity] of patches with [resource-patch? = 1]</metric>
    <metric>sum [consumer-total-Nmol] of turtles</metric>
    <metric>sum-world-N-total</metric>
    <enumeratedValueSet variable="imbalanceddiet_response">
      <value value="&quot;increase_ingestion&quot;"/>
    </enumeratedValueSet>
    <steppedValueSet variable="ingestion_increase" first="0" step="1" last="3"/>
    <enumeratedValueSet variable="target-max-CN-resource-ratio">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-cluster-scatter">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patch-soil-N-multiplier">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-patch-aggregation">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-specific-growth-rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="movement-type">
      <value value="&quot;correlated-random-walk&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="percent-resource-quality-improve">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="turtle-initial-location-spread">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-population">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality-risk">
      <value value="6.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sight">
      <value value="10"/>
    </enumeratedValueSet>
    <steppedValueSet variable="movement-somatic-maintenance-cost-multiplier" first="0" step="0.25" last="1"/>
    <enumeratedValueSet variable="percent-world-no-resources">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="resource-carrying-capacity">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="consumer-body-CN-ratio">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="target-world-avg-resource-CN-ratio">
      <value value="45"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N-reserve-threshold">
      <value value="0.5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="num-resource-clusters" first="0" step="10" last="10"/>
    <enumeratedValueSet variable="distance-consumer-walks">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fraction-structure-shrinkage-allowed">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-fraction-adults">
      <value value="0.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-resource-quantity-as-percent-of-carrying-capacity">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overall-resource-patchiness">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="target-min-CN-resource-ratio">
      <value value="30"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
