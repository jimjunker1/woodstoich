extensions [sr]

;;===========================================================================================================================================================================================;;
;;=========================================================== VARIABLE/PARAMETER DEFINITIONS (Netlogo code formality/requirement) ===========================================================;;
;;===========================================================================================================================================================================================;;
globals [
  ;; World globals
  world-area                                  ; Total number of patches
  world-avg-resource-CN-ratio                 ; Average carbon:nitrogen ratio for resources in entire world
  world-avg-resource-NC-ratio                 ; Average nitrogen:carbon ratio for resources in entire world
  world-total-initial-resource-N              ; Total amount of nitrogen in the world, constrained in closed system
  world-resource-stoichmismatch-severity      ; Factor to increase the carbon:nitrogen ratio in the world by, based on user input
  realized-avg-world-resource-NC-ratio        ; Calculated average nitrogen:carbon ratio in the landscape (post-setup)
  realized-avg-world-resource-CN-ratio        ; Calculated average carbon:nitrogen ratio in the landscape (post-setup)
  avg-resource-quantity                       ; Calculated average of resource amount in landscape
  sum-resource-N-total                        ; Total amount of N in resources variables across landscape
  sum-world-N-total                           ; Total N in the world
  sum-consumer-N-total                        ; Total N in all consumer biomass
  sum-soil-N-total                            ; Total N in all soil pools
  total-world-N                               ; All biomass pools
  ;; Patch related globals
  ;;=== Resource patch cluster globals
  resource-cluster-center-coords              ; The xy-coordinates of center of resource clusters
  resource-patch-radius                       ; Radius of a resource cluster, selected partially by user in interface
  current-per-resource-cluster-area           ; Calculated actual area of a resource cluster
  total-cluster-patches                       ; Number of patches that are part of a resource cluster
  max-cluster-resource-area                   ; Possible maximum number of patches dedicated to resource clusters, estimated from how to maintain world average carbon:nitrogen ratio
  max-per-cluster-area                        ; Maximum per-cluster area, based on max-cluster-resource-area and num-resource-clusters (number of total clusters, set by user interface)
  max-cluster-radius                          ; Maximum radius for each cluster
  realized-cluster-avg-resource-NC-ratio      ; Calculated average resource nitrogen:carbon ratio for cluster patches after setup is complete
  realized-cluster-avg-resource-CN-ratio      ; Calculated average resource carbon:nitrogen ratio for cluster patches after setup is complete
  target-cluster-avg-resource-CN-ratio        ; Goal average resource carbon:nitrogen ratio, used to setup patches
  target-cluster-avg-resource-NC-ratio        ; Goal average resource nitrogen:carbon ratio, used to setup patches


  ;;=== Resource patch NON cluster globals
  total-noncluster-patches                    ; Number of patches that are NOT part of a resource cluster
  noncluster-avg-resource-NC-ratio            ; Target average resource nitrogen:carbon ratio for NON cluster patches after setup is complete
  noncluster-avg-resource-CN-ratio            ; Target average resource carbon:nitrogen ratio for NON cluster patches after setup is complete
  realized-noncluster-avg-resource-NC-ratio   ; Calculated average resource nitrogen:carbon ratio for NONcluster patches after setup is complete
  realized-noncluster-avg-resource-CN-ratio   ; Calculated average resource carbon:nitrogen ratio for NONcluster patches after setup is complete
;  background-soil-N-quantity                 ; how much initial soil N should a patch have? This is slider.
;  patch-soil-N-multiplier                    ; how much should patch-level soil N pool mirror local patch quantity
  resource-CN-optimum                         ; Optimum Carbon:Nitrogen ratio of resources
  resource-NC-optimum                         ; Optimum Nitrogen:Carbon ratio of resources
  resource-C-content                        ; The percent of resource-quantity that is Cmol as a proportion. when resource-C-content = 1, resource-quantity = resource-quantity-Cmol.
; This was set as a global, but could moved to patches-own to make more dynamic
  resource-quantity-min                       ; A minimum value to set the resource quantity to

  ;;=== Soil patch globals

  soil-N-crit                                 ; A critical miniumum soil pool

  ;; Consumer related globals
  consumer-body-NC-ratio                      ; Nitrogen:carbon ratio of consumer's body, inverse of the carbon:nitrogen ratio provided in user interface
]

patches-own [
  ;; Resource (primary production) patch variables
  resource-CN-ratio                           ; Carbon:nitrogen ratio of resource on patch
  resource-NC-ratio                           ; Nitrogen:carbon ratio of resource on patch
  resource-quantity                           ; Total amount of resource per patch, set by the user interface, but this term is reflected in the
  resource-quantity-Cmol                      ; Total amount of resource carbon on a patch
  resource-quantity-Nmol                      ; Total amount of resource nitrogen on a patch
  new-resource-grown                         ; Increase in resource quantity as calculated by the discretized logistic growth model
  pre-growth-resource-quantity                ; Accounting of how much resource was there pre-
  resource-N-uptake                          ; uptake of N from soil pool
  resource-grown-opt-N                      ; resource N demand based on new resources grown at optimal NC ratio
  resource-N-deficit                          ; resource biomass N deficit compared to optimal NC

  ;; Resource cluster-specific patch variables
  cluster-center?                             ; Binary assignment of if patch is a center of a patch
  cluster-patch?                              ; Binary assignment of if patch is part of a patch
  cluster-ids                                 ; Value assigned for each unique cluster. Is a list, each patch can technically allow for overlapping clusters (though it is programmed to avoid them)
  cluster-center-id                           ; Value assigned for each unique cluster center patch
  cluster-quality-level                       ; How "good" resource quality is of patch, based on which ring patch is in around cluster center. Lower number -> lower quality

  ;; Overall patch
  total-patch-quantity-Nmol                   ; The total N moles in patch between resource and soil
  turtle-alterations-patch-pooled-N           ; How much N in moles is left by turtles on patch per turn

  ;; Soil pool
  soil-N-quantity                             ; Soil N pool size
  soil-N-add-quantity                         ; soil pool N addition same as turtle-alterations-patch-pooled-N currently

]

turtles-own [
  ;; Dynamic Energy Budget (DEB) theory terms
  ;;=== State Variables
  M_V                                         ; Somatic body biomass (mol C)
  V                                           ; Somatic body volume (m^3)
  M_Ec                                        ; Biomass of carbon-rich reserve (mol C)
  M_En                                        ; Biomass of nitrogen-rich reserve (mol C)
  Mdens_Ec                                    ; Carbon-rich reserve density (M_EC / M_V)(mol C / mol C)
  Mdens_En                                    ; Nitrogen-rich reserve density (M_EC / M_V)(mol C / mol C)
  M_R                                         ; Biomass of reserves invested in reproduction (mol C)
  M_D                                         ; Biomass of reserves invested in maturity/development (mol C)

  ;;=== Mass Balance Differential Equations
  dF-dt                                       ; Change of food density per tick calculation
  dWc-dt                                      ; Change of carbon egestion per tick calculation
  dWn-dt                                      ; Change of nitrogen egestion per tick calculation
  dM_Ec-dt                                    ; Change of carbon-rich reserve biomass per tick calculation
  dM_En-dt                                    ; Change of nitrogen-rich reserve biomass per tick calculation
  dMdens_Ec-dt                                ; Change of carbon-rich reserve density per tick calculation
  dMdens_En-dt                                ; Change of nitrogen-rich reserve density per tick calculation
  dM_V-dt                                     ; Change of somatic structual biomass per tick calculation
  dM_D-dt                                     ; Change of development biomass per tick calculation
  dM_R-dt                                     ; Change of reproduction biomass per tick calculation
  dXc-dt                                      ; Change of excreted carbon biomass per tick calculation
  dXn-dt                                      ; Change of excreted nitrogen biomass per tick calculation

  ;;=== Flux quantities (Absolute fluxes, per DEB definitions: all units of (mol C / t ) unless otherwise stated)
  ;;======= Ingestion/Assimilation/Waste
  J_FA                                        ; Ingestion
  J_AEc                                       ; Assimilation of food into carbon-rich reserve
  J_AEn                                       ; Assimilation of food into nitrogen-rich reserve
  J_Wc                                        ; Egestion of unassimilated food, carbon flux
  J_Wn                                        ; Egestion of unassimilated food, nitrogen flux (mol N / t)
  ;;======= Reserve dynamics
  J_EcC                                       ; Mobilization of carbon-rich reserve
  J_EnC                                       ; Mobilization of nitrogen-rich reserve
  ;;======= Somatic maintenance
  J_EcS                                       ; Mobilized carbon reserve sent to the somatic maintenance SU
  J_EnS                                       ; Mobilized nitrogen reserve sent to the somatic maintenance SU
  ;;======= Developmental maintenance
  J_EcD                                       ; Mobilized carbon reserve sent to the developmental maintenance SU
  J_EnD                                       ; Mobilized nitrogen reserve sent to the developmental maintenance SU
  ;;======= Somatic growth
  J_EcG                                       ; Mobilized carbon reserve sent to the growth SU
  J_EnG                                       ; Mobilized nitrogen reserve sent to the growth SU
  J-pos_EcG                                   ; Mobilized carbon reserve accepted to bind to the growth SU
  J-pos_EnG                                   ; Mobilized nitrogen reserve accepted to bind to the growth SU
  J-neg_EcG                                   ; Mobilized carbon reserve rejected by the growth SU
  J-neg_EnG                                   ; Mobilized nitrogen reserve rejected by the growth SU
  J_GV                                        ; Structural biomass growth
  ;;======= Development/Reproduction
  J_EcR                                       ; Mobilized carbon reserve sent to the reproduction SU
  J_EnR                                       ; Mobilized nitrogen reserve sent to the reproduction SU
  J-pos_EcR                                   ; Mobilized carbon reserve accepted to bind to the growth SU
  J-pos_EnR                                   ; Mobilized nitrogen reserve accepted to bind to the growth SU
  J-neg_EcR                                   ; Mobilized carbon reserve rejected by the growth SU
  J-neg_EnR                                   ; Mobilized nitrogen reserve rejected by the growth SU
  J_RR                                        ; Development/reproduction buffer biomass growth
  ;;======= Excretion dynamics
  J_SX_Ec                                     ; Excreted carbon-rich reserve from the somatic maintenance SU
  J_SX_En                                     ; Excreted nitrogen-rich reserve from the somatic maintenance SU
  J_SX_C                                      ; Excreted total carbon from the somatic maintenance SU
  J_SX_N                                      ; Excreted total nitrogen from the somatic maintenance SU (mol N / t)
  J_DX_Ec                                     ; Excreted carbon-rich reserve from the developmental maintenance SU
  J_DX_En                                     ; Excreted nitrogen-rich reserve from the developmental maintenance SU
  J_DX_C                                      ; Excreted total carbon from the developmental maintenance SU
  J_DX_N                                      ; Excreted total nitrogen from the developmental maintenance SU (mol N / t)
  J_GX_Ec                                     ; Excreted carbon-rich reserve from the growth SU
  J_GX_En                                     ; Excreted nitrogen-rich reserve from the growth SU
  J_GX_C                                      ; Excreted total carbon from the growth SU
  J_GX_N                                      ; Excreted total nitrogen from the growth SU (mol N / t)
  J_RX_Ec                                     ; Excreted carbon-rich reserve from the reproduction SU
  J_RX_En                                     ; Excreted nitrogen-rich reserve from the reproduction SU
  J_RX_C                                      ; Excreted total carbon from the reproduction SU
  J_RX_N                                      ; Excreted total nitrogen from the reproduction SU (mol N / t)

  ;;=== Parameters: constraining/converting state variables
  M_V_SqBrack                                 ; [M_V]: Volume-specific structural mass (M_V / V) (mol C / m^3)
  M_V_max                                     ; Absolute maximum biomass consumer can reach (mol C)
  M_V_max-to-this-point                       ; The maximum biomass an individual consumer has actually reached (mol C)
  M_V_min-to-die                              ; Minimum somatic biomass compared to M_V_M-to-this-point (mol C)
  M_V_min-to-die-scalar                       ; Fraction of M_V_M-to-this-point that determines (must be between 0 and 1)
  Mdens_Ec_max                                ; Maximum carbon-rich reserve density (mol C)
  Mdens_En_max                                ; Maximum nitrogen-rich reserve density (mol C)
  M_DR_threshold                              ; Biomass needed in development to mature to reproductively active (mol C)
  M_R_bufferbaby                              ; Amount of biomass needed in reproductive buffer to produce juvenille consumer (mol C)
  M_V_baby                                    ; Amount of structural biomass juvenille starts with (mol C)
  M_Ec_baby                                   ; Amount of carbon reserve juvenille starts with (mol C)
  M_En_baby                                   ; Amount of nitrogen reserve juvenille starts with (mol C)
  Mdens_Ec_baby                               ; Density of carbon reserve juvenille starts with (mol C / mol C)
  Mdens_En_baby                               ; Density of nitrogen reserve juvenille starts with (mol C / mol C)
  M_D_baby                                    ; Amount of development juvenille is born with (mol C)

  sigma                                       ; Factor increase of somatic maintenance requirement based on movement
  r-dot                                       ; Specific growth rate of organism (1/t)
  kappa                                       ; Fraction of reserves sent to somatic vs reproductive processes
  rho_Ec                                      ; Fraction of SU-rejected reserve flux recycled back to reserves vs excreted
  rho_En                                      ; Fraction of SU-rejected reserve flux recycled back to reserves vs excreted
  q-F_C                                       ; Molar ratio of carbon in food, normalized to where nitrogen content = 1 (mol)
  q-F_N                                       ; Molar ratio of nitrogen in food, normalized to 1 (mol)
  q-Ec_C                                      ; Molar ratio of carbon in the carbon-rich reserve normalized to itself (1 C : 0 N) (mol)
  q-En_C                                      ; Molar ratio of carbon in the nitrogen-rich reserve normalized to the nitrogen content = 1 (mol)
  q-En_N                                      ; Molar ratio of nitrogen in the nitrogen-rich reserve, normalized to 1 (mol)
  q-V_C                                       ; Molar ratio of carbon in somatic structure, normalized to where nitrogen content = 1 (mol)
  q-V_N                                       ; Molar ratio of nitrogen in somatic structure, normalized to 1 (mol)
  q-R_C                                       ; Molar ratio of carbon in reproduction buffer, normalized to where nitrogen content = 1 (mol)
  q-R_N                                       ; Molar ratio of nitrogen in reproduction buffer, normalized to 1 (mol)

  ;;=== Parameters: relating to ingestion/assimilation
  i_max_CrlBrack                              ; Surface-area specific maximum ingestion rate (mol Food / (V^(2/3) * t))
  F_a_CrlBrack                                ; Type II functional response surface-area specific search rate (1 / (V^(2/3) * t))
  F_h                                         ; Type II functional response half-saturation constant ({i_max} / {F_a}) (mol Food)
  Food                                        ; Amount of resources available on patch
  f                                           ; Type II functional response: (Food / (F_h + Food))
  aEc_max                                     ; Maximum assimilation rate for carbon from food, based on max ingestion rate and yield of reserves from food
  aEn_max                                     ; Maximum assimilation rate for nitrogen from food, based on max ingestion rate and yield of reserves from food
  J_AEn_max_CrlBrack                          ; Surface-area specific maximum assimilation flux to the nitrogen-rich reserve (mol C / (V^(2/3) * t))
  J_AEc_max_CrlBrack                          ; Surface-area specific maximum assimilation flux to the carbon-rich reserve (mol C / (V^(2/3) * t))

  ;;=== Parameters: yield coefficients
  Y_EcF                                       ; Yield of food assimilated to the carbon reserve
  Y_EnF                                       ; Yield of food assimilated to the nitrogen reserve
  Y_SEc                                       ; Yield of somatic maintenance from the carbon reserve
  Y_SEn                                       ; Yield of somatic maintenance from the nitrogen reserve
  Y_DEc                                       ; Yield of developmental maintenance from the carbon reserve
  Y_DEn                                       ; Yield of developmental maintenance from the nitrogen reserve
  Y_VEc                                       ; Yield of structure from the carbon reserve
  Y_VEn                                       ; Yield of structure from the nitrogen reserve
  Y_REc                                       ; Yield of reproduction/development buffer from the carbon reserve
  Y_REn                                       ; Yield of reproduction/development buffer from the nitrogen reserve

  ;;=== Parameters: relating to synthesizing units
  k_M                                         ; Rate of SU turnover for somatic maintenance SU (1/t)
  k_D                                         ; Rate of SU turnover for developmental maintenance SU (1/t)
  k_G                                         ; Rate of SU turnover for somatic growth SU (1/t)
  k_R                                         ; Rate of SU turnover for reproduction/maturity SU (1/t)
  beta-S                                      ; Binding preference for carbon-rich reserve to turnover somatic maintenance SU compared to nitrogen rich reserve (between 0 and 1)
  nu                                          ; enegetic conductance
 ;; Other turtle scenario commands/variables
  did-baby?                                   ; Binary assignment for if consumer reproduced at end of tick
  low-Nreserve-threshold                      ; Value of N reserve density that alerts consumer to find food
  is-Nreserve-low?                            ; Binary assignment for if consumer's N-rich reserve density is low, affects how it chooses to move
  turtle-input-N-this-tick                    ; Total amount of N in moles that is released into the patch at the end of turtle actions per tick

 ;;=== Mortality related parameters
  starved-to-death?                           ; Binary assignment for if consumer starved (M_V fell below M_V_min)
  demographic-stochasticity-death?            ; Binary assignment for if consumer was randomly assigned to die at the end of the tick
  should-die?                                 ; Binary assignment for if consumer should die at the end of tick, based on any death criteria
  corpse-C                                    ; Amount of C released to patch from consumer's corpse
  corpse-N                                    ; Amount of N released to patch from consumer's corpse

  Nmols-eaten                                 ; Temporary variable for test eating mass balance sans DEB model
]


;;==============================================================================================================================================;;
;;=========================================================== USER INTERFACE BUTTONS ===========================================================;;
;;==============================================================================================================================================;;


;;========================== SETUP: INITIALLIZE SIMULATIONS ==========================;;
;;====================================================================================;;

to setup
  init-R-functions
  ca
  clear-output
  set world-area count patches  ; Create global that has total number of patches

  ; Procedures to calculate parameters/variables from user input on the interface
  setup-world-parameters-from-interface  ; General parameters
  setup-resource-parameters-from-interface  ; Parameters specific to resource quality on patches (clusters, etc)

  ; Procedures to initialize model
  setup-patches          ; Initialize patch landscape

  ; Initialize consumer population
  setup-turtles

  set soil-N-crit 1
  reset-ticks
end

to init-R-functions

  sr:setup
  sr:run "source('C:/Users/jrjunker/Documents/Projects/UNT_projects/woodstoich/code/calc_roots.R')"

end
;;=========================== GO: RUN PROCEDURES PER TICK ===========================;;
;;===================================================================================;;

to go
 ; Consumers first move in the simulation, then eat at the patch they arrive at. Then the DEB equations are calculated for the day, and other consumer updates are done
  new-tick-turtle-values                 ; Procedure to ensure all binary variables for turtles are correct
;print "new-tick"
  move                                   ; Procedure to move consumers
;print "move"
  go-DEB-fluxes                          ; Calculates the fluxes for each time step
;print "go-DEB-fluxes"
  go-DEB-ODE                             ; Calculate the change in state variables at each timestep
;print "go-DEB-ODE"
  update-turtle-variables                ; Update the state variables and important parameters of consumer
;print "update-turtles"
  have-baby                              ; Procedure to reproduce
;print "have-baby"
  have-die?                              ; Procedure to determine if turtle dies
;print "have-die"
  patch-N-added-by-turtles
;print "patch-N-add"
  go-die                                 ; Procedure that kills turtles that were indicated to die, after it deposits the N it would leave on the patch from its future corpse
;print "go-die"

  ; Patch commands
  resource-grow                            ; Procedure to increase the quantity of resources per patch
;print "resource-grow"
  update-resource-N-quantity               ; Procedure to calculate how much
;print "update-resrouce-N"
  update-resource-quality                  ; Procedure to calculate how much N should be in the resource that was just grown
;print "update-resource-quality"


  ; Set the color of patches scaled to the nitrogen:carbon ratio (increasing green with increasing resource quality proxy)
  ask patches [
   set pcolor scale-color green resource-NC-ratio 0.5 0
    ]

  make-plots

  tick
end

;;==============================================================================================================================================;;
;;=============================================== CALCULATE PARAMETERS FROM USER INTERFACE INPUT ===============================================;;
;;==============================================================================================================================================;;

;;========================== GENERAL PARAMETERS ==========================;;
;;========================================================================;;

to setup-world-parameters-from-interface
  ; Consumer-related general parameters
  set consumer-body-NC-ratio (1 / consumer-body-CN-ratio)  ; Convert carbon:nitrogen ratio given by user to N:C, for ease in future calculations (e.g., increased N:C ratio scales more intuitively with increased resource quality)

  set world-avg-resource-CN-ratio (target-world-avg-resource-CN-ratio)  ; Define the overall average carbon:nitrogen ratio of the entire landscape
  set world-avg-resource-NC-ratio (1 / world-avg-resource-CN-ratio)     ; Define the overall average nitrogen:carbon ratio of the entire landscape
end


;;========================== RESOURCE-SPECIFIC PARAMETERS ==========================;;
;;==================================================================================;;

to setup-resource-parameters-from-interface
 ; Resource cluster calculations
  ;; Determine what the quality of resources are in the increased quality resource patches
  set target-cluster-avg-resource-CN-ratio (target-min-CN-resource-ratio)  ; Set the intended resource quality for patches based on user input
  set target-cluster-avg-resource-NC-ratio (1 / target-cluster-avg-resource-CN-ratio)
  set resource-NC-optimum target-cluster-avg-resource-NC-ratio

  ;; Calculate the maximum area for each cluster based on the number of clusters set by user. If user sets it to 0, defaults to setting the area to 0.
  ifelse num-resource-clusters = 0[
  ][
    set max-per-cluster-area floor (world-area  / (3 * num-resource-clusters))
    set resource-patch-radius floor sqrt (max-per-cluster-area / pi)  ; Determine the max radius, and take the lower round to avoid overestimation
    set current-per-resource-cluster-area ((resource-patch-radius ^ 2) * pi)
  ]  ; How many patches can actually be within an individual cluster, based on the number of patches selected by user

end



;;========================================================================================================================================================;;
;;=========================================================== SETUP/INITIALIZATION PROCEDURES: ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             PATCHES              ===========================================================;;
;;========================================================================================================================================================;;


;;========================== INITIALIZE RESOURCE PROCEDURES (Primary Production) ==========================;;
;;=========================================================================================================;;

to setup-patches
  ; Set up patches if there are no clusters, as set by user in interface
  ifelse num-resource-clusters = 0  [ask patches [
   set resource-NC-ratio (world-avg-resource-NC-ratio * (1 + (random-float 1 - 0.5) * 2 * noncluster-resource-patchiness))
   set resource-CN-ratio (1 / resource-NC-ratio)
   set pcolor scale-color green resource-NC-ratio 0.5 0
    ]
  ]

  ; If there are clusters, perform the following procedures and commands to set up patches
  [
  setup-resource-clusters   ; Procedure to create clusters of high quality resources
  setup-resouce-nonclusters ; Procedure to create the baseline resource quality of noncluster patches

  ; Set the color of patches scaled to the nitrogen:carbon ratio (increasing green with increasing resource quality proxy)
  ask patches [
   set pcolor scale-color green resource-NC-ratio 0.5 0
    ]
  ]

  ; Calculate the actual average world resource quality (C:N and N:C) to ensure it matches the target
  ask patches[
    set realized-avg-world-resource-NC-ratio mean [resource-NC-ratio] of patches
    set realized-avg-world-resource-CN-ratio (1 / realized-avg-world-resource-NC-ratio)
  ]

  ; Set the initial quantity of resources per patch
  ask patches[
    set resource-quantity (resource-carrying-capacity * initial-resource-quantity-as-percent-of-carrying-capacity) ; Sets the total quantity of resource to be what was defined in the user interface
  ]

  ; Set the intial amount of C and N per patch
  ask patches[
    set resource-quantity-Cmol (resource-quantity)
    set resource-quantity-Nmol (resource-quantity-Cmol * resource-NC-ratio)
  ]

  ask patches[
    setup-soil
  ]

end

;=========================== SETUP SOILS RESOURCES =================================
; ==================================================================================
to setup-soil
  ;====workflow====
  ; setup initial soil N pool from bar
  ; add % proportion of total resource N patch
  set soil-N-quantity background-soil-N-quantity + (patch-soil-N-multiplier * resource-quantity-Nmol)
end

;;========================== SETUP HIGH-QUALITY RESOURCE CLUSTERS ==========================;;
;;==========================================================================================;;

to setup-resource-cluster-center-coords
  ; Calculate the area in which clusters can occur, based on the resource-cluster radius
  let cluster-min-xcor (min-pxcor + resource-patch-radius)
  let cluster-max-xcor (max-pxcor - resource-patch-radius)
  let cluster-min-ycor (min-pycor + resource-patch-radius)
  let cluster-max-ycor (max-pycor - resource-patch-radius)
  set resource-cluster-center-coords []  ; Initialize resource cluster coords list

  ask patches [
    set cluster-center? 0
  ]

  ; Generate a random set of coordinates within the new world area range for the first cluster center
  let first-cluster-xcor cluster-min-xcor + random (cluster-max-xcor - cluster-min-xcor + 1)
  let first-cluster-ycor cluster-min-ycor + random (cluster-max-ycor - cluster-min-ycor + 1)
  let first-cluster-coords list first-cluster-xcor first-cluster-ycor
  set resource-cluster-center-coords lput first-cluster-coords resource-cluster-center-coords

  ask patch first-cluster-xcor first-cluster-ycor [
    set cluster-center? 1
    set pcolor red
  ]

  ; Generate the rest of the coordinate sets (ensure each new pair is not within the radius of a previous coordinate pair)
  while [length resource-cluster-center-coords < num-resource-clusters] [
    let new-cluster-xcor cluster-min-xcor + random (cluster-max-xcor - cluster-min-xcor + 1)
    let new-cluster-ycor cluster-min-ycor + random (cluster-max-ycor - cluster-min-ycor + 1)
    let new-cluster-coord list new-cluster-xcor new-cluster-ycor

    ; Check the distance from the new coordinates to all existing cluster centers
    let too-close? false
    foreach resource-cluster-center-coords [
      existing-coord ->
      let existing-xcor item 0 existing-coord
      let existing-ycor item 1 existing-coord
      ; Calculate distance manually
      let distance-calc sqrt ((new-cluster-xcor - existing-xcor) ^ 2 + (new-cluster-ycor - existing-ycor) ^ 2)
      if distance-calc < (2 * resource-patch-radius) [
        set too-close? true
      ]
    ]

    ; If the new coordinates are valid, add them to the list and mark the patch
    if not too-close? [
      set resource-cluster-center-coords lput new-cluster-coord resource-cluster-center-coords
      ask patch new-cluster-xcor new-cluster-ycor [
        set cluster-center? 1
        set pcolor blue
      ]
    ]
  ]
end


to setup-resource-clusters
setup-resource-cluster-center-coords ; Procedure to determine the locations of clusters without overlapping cluster or being off the edge of the world

 ask patches [
    ; Start all patches with a binary assignment of not being in a cluster
    set cluster-patch? 0

  ; Below, determine the center of clusters and create higher quality patches
    let current-coord list pxcor pycor ; Patch identifies its location

    if member? current-coord resource-cluster-center-coords [ ; If patch coordinate matches a coordinate of a cluster center, it will perform the follow commands to create clusters
      set cluster-patch? 1   ; Binary assignment for if a patch is part of a cluster at all. All centers automatically are in a cluster
      let random-adjusted-cluster-area (current-per-resource-cluster-area * (1.25 + (random-float 0.2 - 0.1)))  ; adjust by +/-10% randomness and add 25% more area
      let random-adjusted-radius (resource-patch-radius * (1 + (random-float 0.2 - 0.2)))  ; adjust by +/-20%
       ask up-to-n-of random-adjusted-cluster-area patches in-radius random-adjusted-radius [
        set resource-NC-ratio (target-cluster-avg-resource-NC-ratio * (1 + (random-float 1 - 0.5) * 2 ))
        set resource-CN-ratio (1 / resource-NC-ratio)
        set cluster-patch? 1
      ]
    ]
  ]

  ; Check mean of cluster patches to see if it is on target average, and calculate a correction
 let cluster-patches patches with [cluster-patch? = 1]

 ask cluster-patches [
   set realized-cluster-avg-resource-NC-ratio mean [resource-NC-ratio] of cluster-patches
   set realized-cluster-avg-resource-CN-ratio (1 / realized-cluster-avg-resource-NC-ratio)
  ]

end


;;========================== SETUP LOWER QUALITY NON-CLUSTER RESOURCES ==========================;;
;;===============================================================================================;;

to setup-resouce-nonclusters
  ; Calculate the total number of noncluster patches available
  set total-cluster-patches count patches with [cluster-patch? = 1] ; Calculate the actual number of patches that are part of a cluster
  set total-noncluster-patches world-area - total-cluster-patches ; Calculate the remaining number of patches, alloted for nonclusters

  ; Calculate what the average N:C ratio should be
  set noncluster-avg-resource-NC-ratio (((world-avg-resource-NC-ratio * world-area) - (total-cluster-patches * realized-cluster-avg-resource-NC-ratio)) / total-noncluster-patches)

  let noncluster-patches patches with [cluster-patch? = 0]

  ; Set the resource-NC-ratio to be the target noncluster average, with randomness
  ask noncluster-patches [
   set resource-CN-ratio 20;(1 / resource-NC-ratio)
  ]

  ; Calculate the actual mean of the noncluster patches
  set realized-noncluster-avg-resource-NC-ratio mean [resource-NC-ratio] of noncluster-patches
  set realized-noncluster-avg-resource-CN-ratio (1 / realized-noncluster-avg-resource-NC-ratio)
end

;;========================================================================================================================================================;;
;;=========================================================== SETUP/INITIALIZATION PROCEDURES: ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             TURTLES              ===========================================================;;
;;========================================================================================================================================================;;

to setup-turtles
  let turtle-cluster-center-x 0  ; X-coordinate of the cluster center
  let turtle-cluster-center-y 0  ; Y-coordinate of the cluster center
  let turtle-cluster-radius (max-pxcor * turtle-initial-location-spread) ; Radius within which turtles will be distributed
  let num-turtles initial-population

  create-turtles num-turtles
  [
    set color blue

    ; Randomly place turtles within the cluster radius
    let angle random-float 360                         ; Random angle in degrees
    let radius random-float turtle-cluster-radius      ; Random distance from the center
    let x turtle-cluster-center-x + radius * cos angle ; X-coordinate within the circle
    let y turtle-cluster-center-y + radius * sin angle ; Y-coordinate within the circle

    setxy x y
    ; Assign turtles' stage
    let num-initial-adults (num-turtles * initial-fraction-adults)


    ; Assign values for DEB equations
    initialize-DEB-parameters
    initialize-DEB-variables
  ]
end


to initialize-DEB-parameters

  ; Yield coefficients
  set Y_EcF 0.6
  set Y_EnF 0.8

  ; Ingestion/Assimilation parameters
  set i_max_CrlBrack 5                           ; Arbitrarily set
  set F_a_CrlBrack 8                               ; Arbitrarily set
  set F_h (i_max_CrlBrack / F_a_CrlBrack)
  set Food [resource-quantity] of patch-here
  set f (Food / (Food + F_h))
  set aEc_max (Y_EcF * i_max_CrlBrack)
  set aEn_max (Y_EnF * i_max_CrlBrack)

  ; Movement parameters
  set sigma movement-somatic-maintenance-cost-multiplier
  ; Reserve fractionation parameters
  set kappa 0.7
  set rho_Ec 0.6
  set rho_En rho_Ec

  ; Stoichiometry parameters
  set q-F_C [resource-CN-ratio] of patch-here
  set q-F_N 1
  set q-Ec_C 1
  set q-En_C 4
  set q-En_N 1
  set q-V_C consumer-body-CN-ratio
  set q-V_N 1
  set q-R_C q-V_C
  set q-R_N 1

  ; SU parameters
  set beta-S 0.2
  set nu 0.4

  ; State variable misc
  set M_DR_threshold  30
  set M_V_SqBrack 4.1

  set k_G 4  ;
  set Y_VEc 0.8
  set Y_VEn Y_VEc
  ; Maintenance parameters
  set k_M 0.45
  set Y_SEc 0.8
  set Mdens_Ec_max 70   ; Maximum carbon-rich reserve density (mol C)
  set Mdens_En_max 50
  ; Development parameters
  set k_D 0.002

  ; Reproduction parameters
  set k_R 0.95
  set Y_REc 0.5
  set Y_REn 0.5
  set M_R_bufferbaby 5
  set M_D_baby 0
end


to initialize-DEB-variables
    set M_V 50
    set V ((M_V / M_V_SqBrack) ^ 3)
    set M_D 30
    set M_R 5
    set M_Ec 10  ; Maximum carbon-rich reserve density (mol C)
    set M_En 10  ; Maximum nitrogen-rich reserve density (mol C)
    set Mdens_Ec (M_Ec / M_V)
    set Mdens_En (M_En / M_V)
end


;;========================================================================================================================================================;;
;;===========================================================   GO/RUN SIMULATION PROCEDURES:  ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             TURTLES              ===========================================================;;
;;========================================================================================================================================================;;


;; ==================== UPDATE TURTLE VALUES THAT MUST BE RESET AT START OF TICK ====================;;
;;===================================================================================================;;

to new-tick-turtle-values
; Initialize important binary variables as 0 before rest of go procedures

ask turtles[
  set starved-to-death? 0
  set demographic-stochasticity-death? 0
  set should-die? 0
  set did-baby? 0

  ; Reset all variable parameters and fluxes to 0 for purposes
  set q-F_C 0
  set Food 0
  set f 0
  set r-dot 0

  set J_FA 0
  set J_AEn 0
  set J_AEc 0
  set J_Wc 0
  set J_Wn 0
  set J_EcC 0
  set J_EnC 0
  set J_EcS 0
  set J_EnS 0
  set J_EcD 0
  set J_EnD 0
  set J_EcG 0
  set J_EnG 0
  set J-pos_EcG 0
  set J-pos_EnG 0
  set J-neg_EcG 0
  set J-neg_EnG 0
  set J_GV 0
  set J_EcR 0
  set J_EnR 0
  set J-pos_EcR 0
  set J-pos_EnR 0
  set J-neg_EcR 0
  set J-neg_EnR 0
  set J_RR 0
  set J_SX_Ec 0
  set J_SX_En 0
  set J_SX_C 0
  set J_SX_N 0
  set J_DX_Ec 0
  set J_DX_En 0
  set J_DX_C 0
  set J_DX_N 0
  set J_GX_Ec 0
  set J_GX_En 0
  set J_GX_C 0
  set J_GX_N 0
  set J_RX_Ec 0
  set J_RX_En 0
  set J_RX_C 0
  set J_RX_N 0

  ; Set all differential equations to 0 for new day
  set dWc-dt 0
  set dWn-dt 0
  set dM_Ec-dt 0
  set dM_En-dt 0
  set dMdens_Ec-dt 0
  set dMdens_En-dt 0
  set dM_V-dt 0
  set dM_D-dt 0
  set dM_R-dt 0
  set dXc-dt 0
  set dXn-dt 0
  ]
end


;; ==================== TURTLE MOVEMENT PROCEDURES ====================;;
;;=====================================================================;;

to move
  ask turtles[
    if consumer-body-CN-ratio > mismatch-threshold
    [
      ;; Seek high nitrogen food when body nitrogen (C:N) is above the threshold
      let best-patch max-one-of patches in-radius 10 [resource-NC-ratio]  ;; Look for best nitrogen patch nearby
      face best-patch  ;; Face the patch with the highest nitrogen
      forward 1 ;; move forward
    ]
    if consumer-body-CN-ratio <= mismatch-threshold
    [
      if movement-type = "brownian-walk" ;; consumer moves in random direction at each time interval
                                [ rt random-float 360
                                  forward 1
      ]
      if movement-type = "correlated-random-walk" ;; consumer moves in a more directional motion at each time interval
                                [ rt random-normal 0 30
                                  forward 1]
    ]
  ]
end

;;======================== DEB MODEL PROCEDURES ======================;;
;;====================================================================;;

;; ========== Calculate fluxes for tick ==========;;
;;================================================;;

to go-DEB-fluxes
ask turtles [
  ; Ingestion
  set q-F_C [resource-CN-ratio] of patch-here
  set q-F_N 1
  set Food [resource-quantity] of patch-here
  set f (Food / (F_H + Food))
  set J_FA ((f * i_max_CrlBrack * ((M_V / M_V_SqBrack) ^ (2 / 3))) * q-F_C)

;   Ensure J_FA does not exceed available resource-quantity
  if J_FA > [resource-quantity] of patch-here [
    set J_FA ([resource-quantity] of patch-here - resource-quantity-min)
  ]

  ;== Remove eaten food from resource quantity on patch and remove the N and C total resource moles thus far eaten
  ask patch-here [
    let max-J_FA [J_FA] of myself  ; Capture J_FA value from the turtle
    ; Ensure resource-quantity never goes below zero
    if max-J_FA > resource-quantity [
      set max-J_FA (resource-quantity - resource-quantity-min)
    ]

    ; Update the resource-quantity and mole values based on the actual amount eaten
    set resource-quantity (resource-quantity - max-J_FA)
      let Cmol-eaten-this-tick  max-J_FA
      let Nmol-eaten-this-tick (max-J_FA );* resource-NC-ratio) ;Jim made an adjustment here, not sure if this is correct

    set resource-quantity-Cmol (resource-quantity-Cmol - Cmol-eaten-this-tick)
    set resource-quantity-Nmol (resource-quantity-Nmol - Nmol-eaten-this-tick)
  ]


    ; Assimilation and egestion
    set J_AEn (Y_EnF  * J_FA * q-En_C * (q-F_N / q-F_C))
    set J_AEc ((Y_EcF  * J_FA) - J_AEn)
    set J_Wc ((1 - Y_EcF) * J_FA)
    set J_Wn ((1 - Y_EnF) * J_FA * (q-F_N / q-F_C))

;    let C-eatbalance-check (J_FA - J_AEc - J_Wc - J_AEn)
;    ifelse C-eatbalance-check = 0 [
;    ][
;      write "C ingestion/assimilation/waste isn't balanced"]
;
;    let N-eatbalance-check (J_FA * ((q-F_N / q-F_C))) - (J_AEn * (q-En_N / q-En_C)) - J_Wn
;    ifelse N-eatbalance-check = 0 [
;    ][
;      write "N ingestion/assimilation/waste isn't balanced"
;    ]

    ; TEMP DYNAMICS PROCEDURE FOR TESTING RESOURCE FEEDBACK
;    temp-growth-dynamics-test
    root-find-reservemobilization-and-rdot
    ; Reserve dynamics
    set J_AEc_max_CrlBrack (aEc_max / (V ^ (2 / 3)))
    set J_AEn_max_CrlBrack (aEn_max / (V ^ (2 / 3)))
    set J_EcC (Mdens_Ec * ((J_AEc_max_CrlBrack / Mdens_Ec_max) - (M_V * r-dot)))
    set J_EnC (Mdens_En * ((J_AEn_max_CrlBrack / Mdens_En_max) - (M_V * r-dot)))

    ; Somatic maintenance
    set J_EcS (1 + sigma) * k_M * M_V
    set J_EnS (1 + sigma) * k_M * M_V

    ; Developmental maintenance
    set J_EcD k_D * M_D
    set J_EnD k_D * M_D

    ; Somatic growth
    set J_EcG ((kappa * J_EcC) - J_EcS)
    set J_EnG ((kappa * J_EnC) - J_EnS)
    set J-pos_EcG (J_GV / Y_VEc)
    set J-pos_EnG (J_GV / Y_VEn)
    set J-neg_EcG (J_EcG - J-pos_EcG)
    set J-neg_EnG (J_EnG - J-pos_EnG)
    set J_GV (1 / ((1 / k_G) + (1 / (J_EcG * Y_VEc)) + (1 / (J_EnG * Y_VEn) ) + (1 / ((J_EcG * Y_VEc) + (J_EnG * Y_VEn))) ))

    ; Development/Reproduction
    set J_EcR (((1 - kappa) * J_EcC) - J_EcD)
    set J_EnR (((1 - kappa) * J_EnC) - J_EnD)
    set J-pos_EcR (J_RR / Y_REc)
    set J-pos_EnR (J_RR / Y_REn)
    set J-neg_EcR (J_EcR - J-pos_EcR)
    set J-neg_EnR (J_EnR - J-pos_EnR)
    set J_RR (1 / ((1 / k_R) + (1 / (J_EcR * Y_REc)) + (1 / (J_EnR * Y_REn) ) + (1 / ((J_EcR * Y_REc) + (J_EnR * Y_REn))) ))

    ; Excretion dynamics
    set J_SX_Ec (k_M * M_V * (1 - Y_SEc))
    set J_SX_En (k_M * M_V * (1 - Y_SEn))
    set J_SX_C (J_SX_Ec + J_SX_En)
    set J_DX_Ec (k_D * M_D * (1 - Y_DEc))
    set J_DX_En (k_D * M_D * (1 - Y_DEn))
    set J_DX_C (J_DX_Ec + J_DX_En)
    set J_DX_N (((q-En_N / q-En_C) * J_DX_En))
    set J_GX_Ec (J-pos_EcG * (1 - Y_VEc))
    set J_GX_En (J-pos_EnG * (1 - Y_VEn))
    set J_GX_C (J_GX_Ec + J_GX_En)
    set J_GX_N ((q-En_N / q-En_C) * J_GX_En)
    set J_RX_Ec (J-pos_EcR * (1 - Y_REc))
    set J_RX_En (J-pos_EnR * (1 - Y_REn))
    set J_RX_C (J_RX_Ec + J_RX_En)
    set J_RX_N ((q-En_N / q-En_C) * J_RX_En)
  ]

end

;; ========== USE ROOT FINDER TO CALCULATE FLUXES ==========;;
;;==========================================================;;
to root-find-reservemobilization-and-rdot
    ;;set parameters for this example
  let parList []
  set parList lput k_G parList
  set parList lput M_V parList
  set parList lput M_V_SqBrack parList
  set parList lput k_M parList
  set parList lput Y_SEc parList
  set parList lput Y_VEc parList
  set parList lput kappa parList
  set parList lput sigma parList
  set parList lput Mdens_Ec_max parList
  set parList lput Mdens_En_max parList
  set parList lput M_Ec parList
  set parList lput M_En parList
  set parList lput nu  parList

;  print parList

  (sr:set "parList" parList)

  sr:run "x <- calc_roots(parList)"
;  print sr:runresult "x"
  set r-dot sr:runresult "x[1]"
  set J_EcC sr:runresult "x[2]"
  set J_EnC sr:runresult "x[3]"
end
;; ========== TEST PROCEDURE TEMPORARY GROWTH DYN ==========;;
;;==========================================================;;
;to temp-growth-dynamics-test
;      ; Calculate the amount of nutrients released into the patch by corpse, if you were to die
;   let body-N_test (M_V * consumer-body-NC-ratio)
;;   set Nmols-eaten (Nmols-eaten + (J_AEn * (q-En_N / q-En_C)))
;;  set corpse-n (body-N_test + Nmols-eaten)
;
;end

;; ========== Calculate change in state variables ==========;;
;;==========================================================;;

to  go-DEB-ODE
  ask turtles[
    ; Differential Equations/ Mass balances
    set dWc-dt (J_Wc)
    set dWn-dt (J_Wn)
    set dM_Ec-dt (J_AEc + (rho_Ec * (J-neg_EcG + J-neg_EcR)) - J_EcC)
    set dM_En-dt (J_AEn + (rho_En * (J-neg_EnG + J-neg_EnR)) - J_EnC)
    set dMdens_Ec-dt (((1 / M_V) * (J_AEc + (rho_Ec * (J-neg_EcG + J-neg_EcR)) - J_EcC)) - (Mdens_Ec * r-dot))
    set dMdens_En-dt (((1 / M_V) * (J_AEn + (rho_En * (J-neg_EnG + J-neg_EnR)) - J_EnC)) - (Mdens_En * r-dot))
    set dM_V-dt (J_GV)
    ifelse (M_D < M_DR_threshold) [
      set dM_D-dt J_RR][
      set dM_D-dt 0]
    ifelse (M_D >= M_DR_threshold) [
      set dM_R-dt J_RR][
      set dM_R-dt 0]
    set dXc-dt (J_GX_C + J_RX_C + (((1 - rho_Ec) + (1 - rho_En)) *(J-neg_EcG + J-neg_EcR)))
    set dXn-dt (J_GX_N + J_RX_N + ((q-En_N / q-En_C) * ((1 - rho_En) *(J-neg_EnG + J-neg_EnR))))

    ; Update other terms that vary
    set r-dot ((1 / M_V) * dM_V-dt)
  ]
end

;; ========== Update turtle model values  ==========;;
;;===============================================;;

to update-turtle-variables

  ;; Update DEB model variables
  ask turtles[
    set M_V (M_V + dM_V-dt)
    set V ((M_V / M_V_SqBrack) ^ 3)
    ifelse ((M_D + dM_D-dt) >= M_DR_threshold)[
      set M_D M_DR_threshold][
      set M_D (M_D + dM_D-dt)
    ]
    if ((M_D + dM_D-dt) >= M_DR_threshold)[
      set M_R (M_R + dM_R-dt)]
    set M_Ec (M_Ec + dM_Ec-dt)
    set M_En (M_En + dM_En-dt)
    set Mdens_Ec (Mdens_Ec + dMdens_Ec-dt)
    set Mdens_En (Mdens_En + dMdens_En-dt)
  ]
end

;; ========== Reproduction submodel  ==========;;
;;=============================================;;

to have-baby
  ask turtles [
    if M_R >= M_R_bufferbaby [
      hatch 1 [
        set M_V M_V_baby
        set V ((M_V / M_V_SqBrack) ^ 3)
        set M_Ec M_Ec_baby
        set M_En M_En_baby
        set Mdens_Ec Mdens_Ec_baby
        set Mdens_En Mdens_En_baby
        set M_D M_D_baby
        set M_R 0


        ; Set important parameters
        set did-baby? 0
        set low-Nreserve-threshold 0.2
        set is-Nreserve-low? 0
        set starved-to-death? 0
        set demographic-stochasticity-death? 0
        set should-die? 0

        ; Below, same parameters that original consumer turtles intialized with
        set sigma movement-somatic-maintenance-cost-multiplier
        set kappa 0.7
        set rho_Ec 0.6
        set rho_En rho_Ec
        set q-F_C [resource-CN-ratio] of patch-here
        set q-F_N 1
        set q-Ec_C 1
        set q-En_C 4
        set q-En_N 1
        set q-V_C consumer-body-CN-ratio
        set q-V_N 1
        set q-R_C q-V_C
        set q-R_N 1
        set beta-S 0.2
      ]
    ]
  ]
end

;; ========== Mortality submodel  ==========;;
;;==========================================;;

to have-die?
  ask turtles [
   ; Starvation scenarios
    if (Mdens_Ec <= 0) [
      set starved-to-death? 1
      set should-die? 1
    ]
    if (Mdens_En <= 0) [
      set starved-to-death? 1
      set should-die? 1
    ]

    if (k_D * M_D < ((1 - kappa) * J_EcC) + ((1 - kappa) * J_EnC))[
      set starved-to-death? 1
      set should-die? 1
    ]

    if (J_EcS < ((1 + sigma) * k_M * M_V ))[
      set starved-to-death? 1
      set should-die? 1
    ]

   ; Demographic stochasiticity
    if random-float 1 < mortality-risk [
      set demographic-stochasticity-death? 1
      set should-die? 1
    ]

  ]
end

;; ========== Nutrient cycling submodel  ==========;;
;;=================================================;;

; Reflect the amount of N available for resource growth that was left by consumer during tick procedures, calculated during patch-N-added-by-turtles

to patch-N-added-by-turtles
  ask turtles [
    ;; The amount of nitrogen (in moles) that is added to the environment is composed of its excretion (X_n) and egestion (W_n) calculated that timestep, plus its corpse if it is to die
    ifelse should-die? = 1 [
      set corpse-n ((q-V_N / q-V_C) * M_V) + ((q-En_N / q-En_C) * M_En) + ((q-R_N / q-R_C) * M_R)
      set turtle-input-N-this-tick ((dXn-dt + dWn-dt) + corpse-N)
      ask patch-here [set turtle-alterations-patch-pooled-N [turtle-input-N-this-tick] of myself]
    ]
    [
      set turtle-input-N-this-tick (dXn-dt + dWn-dt)
      ask patch-here [set turtle-alterations-patch-pooled-N [turtle-input-N-this-tick] of myself]
    ]
  ]

  ;Diffuse the amount of turtle deposited N
; diffuse turtle-alterations-patch-pooled-N 0.5   ; The value to the right is the % of the quantity that gets spread to the patches
end


;; ============== ...Death submodel  ==============;;
;;=================================================;;

to go-die
  ask turtles[
    if should-die? = 1 [
    write "RIP </3"
    die
    ]
  ]
end


;;========================================================================================================================================================;;
;;===========================================================   GO/RUN SIMULATION PROCEDURES:  ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             PATCHES              ===========================================================;;
;;========================================================================================================================================================;;

;;======================================================================================================================;;
;; ============================================ PATCH RESOURCE CALCULATIONS ============================================;;
;;======================================================================================================================;;

;; ================= CALCULATE RESOURCE QUANTITY PER TICK =================;;
;;=========================================================================;;

;to resource-grow
;  ask patches [
;
;    ;; Ensure resource quantity doesn't start below zero
;    if resource-quantity <= 0 [
;      set resource-quantity 0
;    ]
;
;    ;; Grow grass at a constant rate, set by a fraction (arbitrary) of the resource carrying capacity
;    set resource-quantity (resource-quantity + (resource-carrying-capacity * 0.01))
;    ;; Set to carrying capacity if already fully grown
;    if resource-quantity >= resource-carrying-capacity [
;      set resource-quantity resource-carrying-capacity
;    ]
;
;
;  ]
;end
to resource-grow
   set resource-quantity-min 1
  ask patches [
    ;; Check if resource-quantity is 0, if so, apply minimal growth to prevent being stuck at 0
    if resource-quantity = 0 [
      set resource-quantity resource-quantity-min  ;; Small positive value to kickstart growth
    ]

    ;; Logistic growth equation with adjustable growth rate and carrying capacity
    ifelse soil-N-quantity < soil-N-crit [
      set new-resource-grown 0
    ] [
      set new-resource-grown (resource-quantity * resource-specific-growth-rate *
      (resource-carrying-capacity - resource-quantity) / resource-carrying-capacity)
    ]
;    print resource-specific-growth-rate
;    set pre-growth-resource-quantity resource-quantity
    set resource-quantity (resource-quantity + new-resource-grown)

    ;; Set to carrying capacity if already fully grown
;    if resource-quantity >= resource-carrying-capacity [
;      set resource-quantity resource-carrying-capacity
;    ]

;    ;; Ensure resource quantity doesn't go below zero
;    if resource-quantity <= 0 [
;      set resource-quantity resource-quantity-min
;    ]
  ]
;  print [resource-quantity] of patches
end

to update-resource-N-quantity
  ask patches[
    ;; Set minimum N quantity
    if soil-N-quantity <= soil-N-crit [
      set soil-N-quantity soil-N-crit
    ]
    set resource-grown-opt-N (new-resource-grown * resource-NC-optimum) ; what is the growth N demand based on resources growing at optimal NC ratio
    set resource-N-deficit (resource-quantity * resource-NC-optimum) - (resource-quantity * resource-NC-ratio)

    ifelse (resource-grown-opt-N) > soil-N-quantity [    ; If resource N demand is greater than soil N pool, set to available soil N pool - N minimum
    set resource-N-uptake (soil-N-quantity - soil-N-crit)
    ][
     ifelse resource-NC-ratio <= resource-NC-optimum  and (resource-grown-opt-N + resource-N-deficit) > soil-N-quantity[
       set resource-N-uptake (resource-grown-opt-N + (soil-N-quantity - soil-N-crit))
      ][
       set resource-N-uptake new-resource-grown + resource-N-deficit
  ]
    ]
  ]
end

to update-soil-N-quantity
  ask patches[
;    set soil-N-add-quantity
    set soil-N-quantity (soil-N-quantity - resource-N-uptake + turtle-alterations-patch-pooled-N)
  ]
end


to update-resource-quality
 ask patches [
    set resource-quantity-Cmol (resource-quantity)
    set resource-quantity-Nmol (resource-quantity-Nmol + resource-N-uptake)
   if resource-quantity-Nmol = 0 [
    print "ERROR: total amount of N in patch went to 0. Try increasing initial amounts of C and N to avoid this and start simulation again"
    ]
    set resource-NC-ratio (resource-quantity-Nmol / resource-quantity-Cmol)
    set resource-CN-ratio (resource-quantity-Cmol / resource-quantity-Nmol)

;    ; Impose constraints to not have unrealistic values on resource quality
;    if resource-CN-ratio > target-max-CN-resource-ratio [
;      set resource-quantity-Cmol (target-max-CN-resource-ratio * resource-quantity-Nmol)
;      set resource-CN-ratio (resource-quantity-Cmol / resource-quantity-Nmol)
;      set resource-NC-ratio (1 / resource-CN-ratio)
;    ]
;    if resource-CN-ratio < target-min-CN-resource-ratio [
;      set resource-quantity-Cmol (resource-quantity-Nmol * target-min-CN-resource-ratio)
;      set resource-CN-ratio (resource-quantity-Cmol / resource-quantity-Nmol)
;    ]
  ]

end



;;========================================================================================================================================================;;
;;===========================================================   GO/RUN SIMULATION PROCEDURES:  ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================              PLOTS               ===========================================================;;
;;========================================================================================================================================================;;

to make-plots
  plot-resource-quantity
  plot-total-N
  ; Call all plotting procedures
end

to plot-resource-quantity
  ;; Calculate and plot the average resource-quantity
  set avg-resource-quantity mean [resource-quantity] of patches
  set-current-plot "Average Resource Quantity"
  plot avg-resource-quantity
end

to plot-total-N
;; Calculate the sum of the patch variable
  set sum-resource-N-total sum [resource-quantity-Nmol] of patches
  set sum-consumer-N-total (sum [ ((q-V_N / q-V_C) * M_V) + ((q-En_N / q-En_C) * M_En)] of turtles)
  set sum-soil-N-total (sum [soil-N-quantity] of patches)
  set sum-world-N-total (sum-resource-N-total + sum-consumer-N-total + sum-soil-N-total)


  ;; Update the plot with the sum
  set-current-plot "Total N Dynamics"  ;; Replace with the name of your plot
;  plot sum-consumer-N-total
;  plot sum-soil-N-total
;  plot sum-resource-N-total
;;  print sum-world-N-total
;  plot sum-world-N-total

end
@#$#@#$#@
GRAPHICS-WINDOW
12
10
668
667
-1
-1
6.42
1
10
1
1
1
0
0
0
1
-50
50
-50
50
1
1
1
ticks
30.0

SLIDER
887
224
1059
257
num-resource-clusters
num-resource-clusters
0
15
6.0
1
1
NIL
HORIZONTAL

SLIDER
671
176
843
209
initial-population
initial-population
0
100
98.0
1
1
NIL
HORIZONTAL

BUTTON
670
10
734
44
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
670
70
825
130
consumer-body-CN-ratio
6.0
1
0
Number

SLIDER
674
367
869
400
turtle-initial-location-spread
turtle-initial-location-spread
0
1
0.9
0.1
1
NIL
HORIZONTAL

BUTTON
734
10
798
44
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

MONITOR
463
716
683
761
Target world mean resource N:C ratio
world-avg-resource-NC-ratio
5
1
11

MONITOR
462
819
681
860
Target high-quality cluster resource N:C ratio
target-cluster-avg-resource-NC-ratio
5
1
10

MONITOR
462
863
681
904
True high-quality cluster resource N:C ratio
realized-cluster-avg-resource-NC-ratio
5
1
10

MONITOR
462
918
682
959
Target non-cluster resource N:C ratio
noncluster-avg-resource-NC-ratio
5
1
10

MONITOR
236
718
454
763
Target world mean resource C:N ratio
world-avg-resource-CN-ratio
5
1
11

MONITOR
462
963
683
1004
True non-cluster resource N:C ratio
realized-noncluster-avg-resource-NC-ratio
5
1
10

SLIDER
887
324
1072
357
noncluster-resource-patchiness
noncluster-resource-patchiness
0
0.8
0.8
0.1
1
NIL
HORIZONTAL

MONITOR
235
965
456
1006
True non-cluster resource C:N ratio
realized-noncluster-avg-resource-CN-ratio
5
1
10

MONITOR
462
763
682
808
True world mean resource N:C ratio
realized-avg-world-resource-NC-ratio
5
1
11

MONITOR
236
765
455
810
True world mean resource C:N ratio
realized-avg-world-resource-CN-ratio
5
1
11

SLIDER
672
296
844
329
mismatch-threshold
mismatch-threshold
1
50
50.0
1
1
NIL
HORIZONTAL

CHOOSER
672
248
847
293
movement-type
movement-type
"brownian-walk" "correlated-random-walk"
1

SLIDER
673
403
925
436
movement-somatic-maintenance-cost-multiplier
movement-somatic-maintenance-cost-multiplier
0
10
10.0
1
1
NIL
HORIZONTAL

SLIDER
672
210
844
243
initial-fraction-adults
initial-fraction-adults
0
1
1.0
0.1
1
NIL
HORIZONTAL

SLIDER
835
101
1012
134
resource-specific-growth-rate
resource-specific-growth-rate
0.001
2
0.011
0.005
1
NIL
HORIZONTAL

SLIDER
836
66
1008
99
resource-carrying-capacity
resource-carrying-capacity
20
200
200.0
5
1
NIL
HORIZONTAL

TEXTBOX
675
51
825
69
Turtle Controls
14
0.0
1

MONITOR
671
130
809
175
NIL
consumer-body-NC-ratio
5
1
11

TEXTBOX
903
10
1115
50
Resource Quantity Controls
14
0.0
1

TEXTBOX
888
200
1076
223
Resource Quality Controls
14
0.0
1

MONITOR
233
823
452
864
Target high-quality cluster resource C:N ratio
target-cluster-avg-resource-CN-ratio
5
1
10

MONITOR
232
868
452
909
True high-quality cluster resource C:N ratio
realized-cluster-avg-resource-CN-ratio
5
1
10

MONITOR
233
921
455
962
Target non-cluster resource C:N ratio
(1 / noncluster-avg-resource-NC-ratio)
5
1
10

TEXTBOX
80
687
637
707
Resource Moniter for Troubleshooting(Will be deleted for final product)
14
0.0
1

SLIDER
809
30
1134
63
initial-resource-quantity-as-percent-of-carrying-capacity
initial-resource-quantity-as-percent-of-carrying-capacity
0.01
1
0.51
0.1
1
NIL
HORIZONTAL

PLOT
24
728
224
878
Average Resource Quantity
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
"default" 1.0 0 -12087248 true "" "plot avg-resource-quantity"

SLIDER
673
331
846
364
mortality-risk
mortality-risk
0
0.01
0.0
0.0001
1
NIL
HORIZONTAL

PLOT
1100
230
1502
586
Total N Dynamics
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"World Total N" 1.0 0 -955883 true "" "plot sum-world-N-total"
"Resource N" 1.0 0 -13840069 true "" "plot sum-resource-N-total"
"Consumer N" 1.0 0 -13345367 true "" "plot sum-consumer-N-total"
"Soil N" 1.0 0 -10146808 true "" "plot sum-soil-N-total"

INPUTBOX
1149
135
1307
195
target-min-CN-resource-ratio
20.0
1
0
Number

INPUTBOX
1148
72
1307
132
target-max-CN-resource-ratio
40.0
1
0
Number

INPUTBOX
1145
10
1300
70
target-world-avg-resource-CN-ratio
30.0
1
0
Number

PLOT
1296
421
1496
571
Total World Turtles
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
"default" 1.0 0 -16777216 true "" "plot count turtles"

SLIDER
886
291
1075
324
background-soil-N-quantity
background-soil-N-quantity
1
100
10.0
1
1
NIL
HORIZONTAL

SLIDER
886
258
1083
291
patch-soil-N-multiplier
patch-soil-N-multiplier
0
20
6.0
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
