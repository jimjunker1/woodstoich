extensions [sr palette]

breed [cluster-makers cluster-maker]   ; Define a special breed for cluster-making turtles

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
  current-per-resource-cluster-area           ; Calculated actual area of a resource cluster
  total-resource-cluster-patches              ; Number of patches that are part of a resource cluster
  absolute-max-cluster-resource-area          ; Calculated absolute maximum area of world that can be clusters
  max-cluster-resource-area                   ; Possible maximum number of patches dedicated to resource clusters, estimated from how to maintain world average carbon:nitrogen ratio
  noncluster-resource-area                    ; Minimum possible total patches that are allocated for non-clustered resources (low quality)
  max-per-cluster-area                        ; Maximum per-cluster area, based on max-cluster-resource-area and num-resource-clusters (number of total clusters, set by user interface)
  max-cluster-radius                          ; Maximum radius for each cluster
  realized-cluster-avg-resource-NC-ratio      ; Calculated average resource nitrogen:carbon ratio for cluster patches after setup is complete
  realized-cluster-avg-resource-CN-ratio      ; Calculated average resource carbon:nitrogen ratio for cluster patches after setup is complete
  target-cluster-avg-resource-CN-ratio        ; Goal average resource carbon:nitrogen ratio, used to setup patches
  target-cluster-avg-resource-NC-ratio        ; Goal average resource nitrogen:carbon ratio, used to setup patches
  target-noncluster-avg-resource-CN-ratio     ; Goal average resource carbon:nitrogen ratio, used to setup patches
  target-noncluster-avg-resource-NC-ratio     ; Goal average resource nitrogen:carbon ratio, used to setup patches
  min-cn-color-value
  max-cn-color-value


  ;;=== Resource patch NON cluster globals
  total-noncluster-patches                    ; Number of patches that are NOT part of a resource cluster
  noncluster-avg-resource-NC-ratio            ; Target average resource nitrogen:carbon ratio for NON cluster patches after setup is complete
  noncluster-avg-resource-CN-ratio            ; Target average resource carbon:nitrogen ratio for NON cluster patches after setup is complete
  realized-noncluster-avg-resource-NC-ratio   ; Calculated average resource nitrogen:carbon ratio for NONcluster patches after setup is complete
  realized-noncluster-avg-resource-CN-ratio   ; Calculated average resource carbon:nitrogen ratio for NONcluster patches after setup is complete
  resource-quantity-min                       ; A minimum value to set the resource quantity to
  total-resource-quantity

  ;;=== Soil patch globals

  soil-N-crit                                 ; A critical miniumum soil pool

  ;; Consumer related globals
  consumer-body-NC-ratio                      ; Nitrogen:carbon ratio of consumer's body, inverse of the carbon:nitrogen ratio provided in user interface
  total-consumers
  total-adult-consumers
  total-juvenille-consumers
]

patches-own [
  ;; Resource (primary production) patch variables
  resource-CN-ratio                           ; Carbon:nitrogen ratio of resource on patch
  resource-NC-ratio                           ; Nitrogen:carbon ratio of resource on patch
  resource-quantity                           ; Total amount of resource per patch, set by the user interface, but this term is reflected in the
  resource-quantity-Cmol                      ; Total amount of resource carbon on a patch
  resource-quantity-Nmol                      ; Total amount of resource nitrogen on a patch
  max-resource-grown
  new-resource-grown                          ; Increase in resource quantity as calculated by the discretized logistic growth model
  new-resource-quantity-Cmol
  new-resource-quantity-Nmol
  max-new-resource-Nmol
  max-new-resource-Cmol

  pre-growth-resource-quantity                ; Accounting of how much resource was there pre-
  resource-N-uptake                           ; uptake of N from soil pool
  resource-grown-opt-N                        ; resource N demand based on new resources grown at optimal NC ratio
  resource-N-deficit                          ; resource biomass N deficit compared to optimal NC

  ;; Resource cluster-specific patch variables
  cluster-center?                             ; Binary assignment of if patch is a center of a patch
  cluster-patch?                              ; Binary assignment of if patch is part of a patch

  ;; Overall patch
  turtle-alterations-patch-pooled-N           ; How much N in moles is left by turtles on patch per turn

  ;; Soil pool
  soil-N-quantity                             ; Soil N pool size
  soil-N-add-quantity                         ; soil pool N addition same as turtle-alterations-patch-pooled-N currently

]

turtles-own [
  ;;=== State Variables
  M_V                                         ; Somatic body biomass (mol C)
  V                                           ; Somatic body volume (m^3)
  M_Ec                                        ; Biomass of carbon-rich reserve (mol C)
  M_En                                        ; Biomass of nitrogen-rich reserve (mol C)
  Mdens_Ec                                    ; Carbon-rich reserve density (M_EC / M_V)(mol C / mol C)
  scaled_Mdens_Ec                             ; Scaled carbon-rich reserve density (range 0-1)
  Mdens_En                                    ; Nitrogen-rich reserve density (M_EC / M_V)(mol C / mol C)
  scaled_Mdens_En                             ; Scaled nitrogen-rich reserve density (range 0-1)
  M_R                                         ; Biomass of reserves invested in reproduction (mol C)
  M_H                                         ; Biomass of reserves invested in maturity/development (mol C)

  ;;=== Mass Balance Differential Equations
  dWc-dt                                      ; Change of carbon egestion per tick calculation
  dWn-dt                                      ; Change of nitrogen egestion per tick calculation
  dM_Ec-dt                                    ; Change of carbon-rich reserve biomass per tick calculation
  dM_En-dt                                    ; Change of nitrogen-rich reserve biomass per tick calculation
  dMdens_Ec-dt                                ; Change of carbon-rich reserve density per tick calculation
  dMdens_En-dt                                ; Change of nitrogen-rich reserve density per tick calculation
  dM_V-dt                                     ; Change of somatic structual biomass per tick calculation
  dM_H-dt                                     ; Change of development biomass per tick calculation
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
  J_GV                                        ; Structural biomass growth
  J_VS
  J_EcG_plus
  J_EnG_plus
  J_EcG_neg
  J_EnG_neg

  ;;======= Development/Reproduction
  J_EcR                                       ; Mobilized carbon reserve sent to the reproduction SU
  J_EnR                                       ; Mobilized nitrogen reserve sent to the reproduction SU
  J_RR                                        ; Development/reproduction buffer biomass growth
  J_HD
  J_EcR_plus
  J_EnR_plus
  J_EcR_neg
  J_EnR_neg

  ;;======= Excretion dynamics
  J_GREc
  J_GREn
  J_GRX_Ec
  J_GRX_En

  ; Initial state variable conditions for juvenilles at birth
  M_V_baby                                    ; Amount of structural biomass juvenille starts with (mol C)
  M_Ec_baby                                   ; Amount of carbon reserve juvenille starts with (mol C)
  M_En_baby                                   ; Amount of nitrogen reserve juvenille starts with (mol C)
  Mdens_Ec_baby                               ; Density of carbon reserve juvenille starts with (mol C / mol C)
  Mdens_En_baby                               ; Density of nitrogen reserve juvenille starts with (mol C / mol C)
  M_H_baby                                    ; Amount of development juvenille is born with (mol C)


  ; Maintenance parameters
  sigma                                       ; Factor increase of somatic maintenance requirement based on movement
  k_M
  k_J

  ; Somatic structure/growth parameters
  r_dot                                       ; Specific growth rate of organism (1/t)
  M_V_SqBrack                                 ; [M_V]: Volume-specific structural mass (M_V / V) (mol C / m^3)
  M_V_max                                     ; Absolute maximum biomass consumer can reach (mol C)

  ; Reserve parameters
  kappa                                       ; Fraction of reserves sent to somatic vs reproductive processes
  kappaG_Ec
  kappaG_En
  kappaR_Ec
  kappaR_En
  nu
  Mdens_Ec_max                                ; Maximum carbon-rich reserve density (mol C)
  Mdens_En_max                                ; Maximum nitrogen-rich reserve density (mol C)
  mdens_tot_max
  mdens_tot_max_SqBrack
  mdens_Ec_max_SqBrack
  mdens_En_max_SqBrack

  ; Stoichiometry parameters
  q-F_C                                       ; Molar ratio of carbon in food, normalized to where nitrogen content = 1 (mol)
  q-F_N                                       ; Molar ratio of nitrogen in food, normalized to 1 (mol)
  q-Ec_C                                      ; Molar ratio of carbon in the carbon-rich reserve normalized to itself (1 C : 0 N) (mol)
  q-En_C                                      ; Molar ratio of carbon in the nitrogen-rich reserve normalized to the nitrogen content = 1 (mol)
  q-Ec_N                                      ; Molar ratio of carbon in the nitrogen-rich reserve, normalized to 1 (mol)
  q-En_N                                      ; Molar ratio of nitrogen in the nitrogen-rich reserve, normalized to 1 (mol)
  q-V_C                                       ; Molar ratio of carbon in somatic structure, normalized to where nitrogen content = 1 (mol)
  q-V_N                                       ; Molar ratio of nitrogen in somatic structure, normalized to 1 (mol)
  q-R_C                                       ; Molar ratio of carbon in reproduction buffer, normalized to where nitrogen content = 1 (mol)
  q-R_N                                       ; Molar ratio of nitrogen in reproduction buffer, normalized to 1 (mol)
  n_V
  n_R
  n-V_En
  n-R_En
  n-V_Ec
  n-R_Ec

  ; Yield coefficients parameters
  Y_EcF                                       ; Yield of food assimilated to the carbon reserve
  Y_EnF                                       ; Yield of food assimilated to the nitrogen reserve
  Y_VS
  Y_HRD
  Y_EcV                                       ; Yield of structure from the carbon reserve
  Y_EnV                                       ; Yield of structure from the nitrogen reserve
  Y_EcR                                       ; Yield of reproduction/development buffer from the carbon reserve
  Y_EnR                                       ; Yield of reproduction/development buffer from the nitrogen reserve

  ; Maturity and reproduction parameters
  M_H_p                              ; Biomass needed in development to mature to reproductively active (mol C)
  M_R_bufferbaby                              ; Amount of biomass needed in reproductive buffer to produce juvenille consumer (mol C)

  ; Ingestion/assimilation parameters
  i_max_CrlBrack                              ; Surface-area specific maximum ingestion rate (mol Food / (V^(2/3) * t))
  F_a_CrlBrack                                ; Type II functional response surface-area specific search rate (1 / (V^(2/3) * t))
  F_h                                         ; Type II functional response half-saturation constant ({i_max} / {F_a}) (mol Food)
  Food                                        ; Amount of resources available on patch
  f                                           ; Type II functional response: (Food / (F_h + Food))
  temp-N-eaten


  ; Synthesizing Units parameters
  k_S                                         ; Rate of SU turnover for somatic maintenance SU (1/t)
  k_D                                         ; Rate of SU turnover for developmental maintenance SU (1/t)
  k_G                                         ; Rate of SU turnover for somatic growth SU (1/t)
  k_R                                         ; Rate of SU turnover for reproduction/maturity SU (1/t)
  rho_S_N
  rho_S_V
  rho_D_N
  rho_D_H

 ; Other turtle scenario commands/variables
  is-adult?                                   ; Binary assignment to indicate consumer lifestage
  did-baby?                                   ; Binary assignment for if consumer reproduced at end of tick
  low-Nreserve-threshold                      ; Value of N reserve density that alerts consumer to find food
  is-Nreserve-low?                            ; Binary assignment for if consumer's N-rich reserve density is low, affects how it chooses to move
  turtle-input-N-this-tick                    ; Total amount of N in moles that is released into the patch at the end of turtle actions per tick

 ; Mortality related parameters
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

  ; Procedures to set parameters/variables
  setup-parameters-from-interface

  ; Procedures to initialize model
  setup-patches          ; Initialize patch landscape

  ; Initialize consumer population
  setup-turtles

  reset-ticks
end

to init-R-functions

  sr:setup
  sr:run "source('C:/Users/lynda/Documents/Work/Conferences/2024/Woodstoich/Model/NetLogo/rootfinder_woodstoich_netlogofile.R')"

end
;;=========================== GO: RUN PROCEDURES PER TICK ===========================;;
;;===================================================================================;;

to go
 ; Consumers first move in the simulation, then eat at the patch they arrive at. Then the DEB equations are calculated for the day, and other consumer updates are done
  reset-ODE-and-fluxes
  new-tick-turtle-values                 ; Procedure to ensure all binary variables for turtles are correct
  move                                   ; Procedure to move consumers
  go-DEB-flux-ingest-assimilate-egest
  go-DEB-fluxes                          ; Calculates the fluxes for each time step
  go-DEB-ODE                             ; Calculate the change in state variables at each timestep
  update-turtle-variables                ; Update the state variables and important parameters of consumer
;  have-baby                              ; Procedure to reproduce
;  have-die?                              ; Procedure to determine if turtle dies
;  patch-N-added-by-turtles
;  go-die                                 ; Procedure that kills turtles that were indicated to die, after it deposits the N it would leave on the patch from its future corpse


  ; Patch commands
  update-patches                           ; Procedure to increase the quantity of resources per patch
  update-patch-misc


  make-plots

  tick
end

;;==============================================================================================================================================;;
;;=============================================== CALCULATE PARAMETERS FROM USER INTERFACE INPUT ===============================================;;
;;==============================================================================================================================================;;

;;========================== GENERAL PARAMETERS ==========================;;
;;========================================================================;;

to setup-parameters-from-interface
  ; Consumer-related general parameters
  set consumer-body-NC-ratio (1 / consumer-body-CN-ratio)  ; Convert carbon:nitrogen ratio given by user to N:C, for ease in future calculations (e.g., increased N:C ratio scales more intuitively with increased resource quality)

  set world-avg-resource-CN-ratio (target-world-avg-resource-CN-ratio)  ; Define the overall average carbon:nitrogen ratio of the entire landscape
  set world-avg-resource-NC-ratio (1 / world-avg-resource-CN-ratio)     ; Define the overall average nitrogen:carbon ratio of the entire landscape

 ; Resource cluster calculations
  ;; Determine what the quality of resources are in the increased quality resource patches
  set target-cluster-avg-resource-CN-ratio (target-min-CN-resource-ratio)  ; Set the intended resource quality for patches based on user input
  set target-cluster-avg-resource-NC-ratio (1 / target-cluster-avg-resource-CN-ratio)
  set target-noncluster-avg-resource-CN-ratio target-max-CN-resource-ratio
  set target-noncluster-avg-resource-NC-ratio (1 / target-noncluster-avg-resource-CN-ratio)


 ; Determine size of resource clusters from user input
  ;; Calculate the absolute maximum area that can be dedicated to resource clusters, based on the minimum, maximum, and average world nutrient ratios provided by user
  ifelse (num-resource-clusters = 0) [
    set max-cluster-resource-area 0
    set max-per-cluster-area 0
  ] [
    set max-cluster-resource-area floor ((world-area * (world-avg-resource-CN-ratio - target-noncluster-avg-resource-CN-ratio)) / (target-cluster-avg-resource-CN-ratio - target-noncluster-avg-resource-CN-ratio))
    set max-per-cluster-area floor (max-cluster-resource-area / num-resource-clusters)
  ]
end



;;========================================================================================================================================================;;
;;=========================================================== SETUP/INITIALIZATION PROCEDURES: ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             PATCHES              ===========================================================;;
;;========================================================================================================================================================;;


;;========================== INITIALIZE RESOURCE PROCEDURES (Primary Production) ==========================;;
;;=========================================================================================================;;

to setup-patches

;; SETUP RESOURCE QUALITY

  ; Set up patches if there are no clusters, as set by user in interface
  ifelse num-resource-clusters = 0  [ask patches [
    set resource-CN-ratio ( world-avg-resource-CN-ratio * (1 + (random-float (2 * overall-resource-patchiness) - overall-resource-patchiness)))
    ]
  ]

  ; If there are clusters, perform the following procedures and commands to set up patches
  [
    setup-resource-cluster-center-coords
    setup-resource-clusters   ; Procedure to create clusters of high quality resources
    setup-resouce-nonclusters ; Procedure to create the baseline resource quality of noncluster patches
  ]

  ; Set nitrogen:carbon ratios for future calculations
  ask patches [
    set resource-NC-ratio (1 / resource-CN-ratio)
  ]

;; SETUP RESOURCE QUANTITY

  ; Set the initial quantity of resources per patch
  ask patches[
    set resource-quantity (resource-carrying-capacity * initial-resource-quantity-as-percent-of-carrying-capacity) ; Sets the total quantity of resource to be what was defined in the user interface

  ; Set the intial amount of C and N per patch
    set resource-quantity-Cmol (resource-quantity)
    set resource-quantity-Nmol (resource-quantity-Cmol * resource-NC-ratio)
  ]

 ;; SETUP SOIL NITROGEN POOL
  setup-soil

 ;; POST PATCH SETUP CALCULATIONS AND AESTHETIC SETTINGS

  ; Calculate the actual average world, cluster, and non-cluster patch resource quality (C:N and N:C) to ensure it matches the target
  ifelse (num-resource-clusters = 0) [
    set realized-avg-world-resource-CN-ratio mean [resource-CN-ratio] of patches
    set realized-avg-world-resource-NC-ratio (1 / realized-avg-world-resource-CN-ratio)
  ]
  [
    set realized-avg-world-resource-CN-ratio mean [resource-CN-ratio] of patches
    set realized-avg-world-resource-NC-ratio (1 / realized-avg-world-resource-CN-ratio)

    set realized-noncluster-avg-resource-CN-ratio mean [resource-CN-ratio] of patches with [cluster-patch? = 0]
    set realized-noncluster-avg-resource-NC-ratio (1 / realized-noncluster-avg-resource-CN-ratio)

    set realized-cluster-avg-resource-CN-ratio mean [resource-CN-ratio] of patches with [cluster-patch? = 1]
    set realized-cluster-avg-resource-NC-ratio (1 / realized-cluster-avg-resource-CN-ratio)
    ]

  set min-cn-color-value (target-max-CN-resource-ratio + target-min-CN-resource-ratio)
  set max-cn-color-value (target-min-CN-resource-ratio - target-min-CN-resource-ratio)

color-patches

end

to color-patches
  ask patches [
  set pcolor palette:scale-gradient [[84 83 60] [96 94 56] [109 106 59] [121 115 63] [133 125 68] [123 132 69] [112 140 74] [103 154 81] [94 167 88]] resource-CN-ratio min-cn-color-value max-cn-color-value
]
end

;;========================== SETUP HIGH-QUALITY RESOURCE CLUSTERS ==========================;;
;;==========================================================================================;;

to setup-resource-cluster-center-coords
  set resource-cluster-center-coords []  ; Initialize resource cluster coords list
  let first-cluster-coords list random-pxcor random-pycor
  set resource-cluster-center-coords lput first-cluster-coords resource-cluster-center-coords

    ; Generate the rest of the coordinate sets (ensure each new pair is not within the radius of a previous coordinate pair)
  while [length resource-cluster-center-coords < num-resource-clusters] [
    let new-cluster-xcor random-pxcor
    let new-cluster-ycor random-pycor
    let new-cluster-coord list new-cluster-xcor new-cluster-ycor

    set resource-cluster-center-coords lput new-cluster-coord resource-cluster-center-coords
  ]

  ; Tag the patch for each center
  ask patches [
    ; Start all patches with a binary assignment of not being in a cluster
    set cluster-patch? 0

  ; Below, determine the center of clusters and create higher quality patches
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
  let cluster-coords resource-cluster-center-coords
  let num-coords length cluster-coords           ; Number of cluster-makers to create

  create-cluster-makers num-coords [
    let coord item who cluster-coords            ; Assign each cluster-maker to a coordinate from the list
    setxy first coord last coord                 ; Set the cluster-maker's position based on the coordinates generated in previous procedure (setup-resource-cluster-center-coords)
  ]

  ; Now have each cluster-maker perform the random walk and modify patches
  ask cluster-makers [
    repeat (max-per-cluster-area) [                                  ; Size of one cluster, normalized so that the patch the turtle is currently on won't overshoot the max cluster area
      ask min-one-of patches with [ cluster-patch? = 0 ] [ distance myself ] [
        set cluster-patch? 1
        set resource-CN-ratio (target-cluster-avg-resource-CN-ratio * (1 + (random-float (2 * overall-resource-patchiness) - overall-resource-patchiness)))
      ]
      rt random 360
      fd (0.5 + resource-cluster-scatter)
    ]
    die                                           ; Remove the cluster-maker after it finishes
  ]

  ; It is possible that the clusters may not be exactly the resource cluster area. For simplicity and computational efficiency, this is corrected after the cluster formation:
  set total-resource-cluster-patches count patches with [cluster-patch? = 1]                                                     ; Sum all patches in a resource cluster
  let discrepancy-resource-cluster-area (max-cluster-resource-area - total-resource-cluster-patches)                             ; Calculate the difference between the actual resource cluster patches and intended

  let resource-cluster-patches patches with [cluster-patch? = 1]                                                                 ; Temp assignment to any cluster patches
  let resource-cluster-border-patches resource-cluster-patches with [count neighbors with [cluster-patch? = 1] < 8]              ; Find cluster patches with fewer than 8 neighbors that are also in a cluster
  let resource-cluster-edge-patches resource-cluster-border-patches with [count neighbors with [cluster-patch? = 0] >= 3]        ; Find the cluster patches with at least one neighbor that isn't in a cluster, and find patches that have at least 3 neighbors that are in a cluster to mimic morerealistic border

  let resource-noncluster-patches patches with [cluster-patch? = 0]                                                              ; Temp assignment to any noncluster patches
  let resource-noncluster-border-patches resource-noncluster-patches with [count neighbors with [cluster-patch? = 1] > 0]        ; Finds noncluster patches with at least one neighbor that is in a cluster
  let resource-noncluster-edge-patches resource-noncluster-border-patches with [count neighbors with [cluster-patch? = 1] >= 3]  ; Find the cluster patches with at least one neighbor that isn't in a cluster, and find patches that have at least 3 neighbors that are in a cluster to mimic morerealistic border

  if total-resource-cluster-patches > max-cluster-resource-area [
    ask n-of (abs discrepancy-resource-cluster-area) resource-cluster-edge-patches [
      set resource-CN-ratio target-noncluster-avg-resource-CN-ratio
      set cluster-patch? 0
      set total-resource-cluster-patches count patches with [cluster-patch? = 1]
    ]
  ]

  if total-resource-cluster-patches < max-cluster-resource-area[
    ask n-of (abs discrepancy-resource-cluster-area) resource-noncluster-edge-patches [
      set resource-CN-ratio target-cluster-avg-resource-CN-ratio
      set cluster-patch? 1
      set total-resource-cluster-patches count patches with [cluster-patch? = 1]
    ]
  ]

end

;;========================== SETUP LOWER QUALITY NON-CLUSTER RESOURCES ==========================;;
;;===============================================================================================;;

to setup-resouce-nonclusters
 let noncluster-patches patches with [cluster-patch? = 0]

  ; Set the resource-NC-ratio to be the target noncluster average, with randomness
  ask noncluster-patches [
    set resource-CN-ratio (target-noncluster-avg-resource-CN-ratio * (1 + (random-float (2 * overall-resource-patchiness) -  overall-resource-patchiness)))
  ]
end


;=========================== SETUP SOILS RESOURCES =================================
; ==================================================================================
to setup-soil
  ;====workflow====
  ; setup initial soil N pool from bar
  ; add % proportion of total resource N patch
  ask patches [
    set soil-N-quantity (resource-quantity-Nmol * patch-soil-N-multiplier)
  ]

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
  ]

  ; Assign values for DEB equations
  initialize-DEB-parameters
  initialize-DEB-variables

end


to initialize-DEB-parameters

 ; SETUP POPULATION-CONSTANT PARAMETERS
  ask turtles[

    ; Stoichiometry parameters
    set q-F_C [resource-CN-ratio] of patch-here
    set q-F_N 1
    set q-Ec_C 1
    set q-Ec_N 0
    set q-En_C 4
    set q-En_N 1
    set q-V_C consumer-body-CN-ratio
    set q-V_N 1
    set q-R_C q-V_C
    set q-R_N 1

    ; Moles of state variables to calculate yields
    set n_V 1
    set n_R 1
    set n-V_En 1
    set n-R_En 1
    set n-V_Ec  (n_V * (q-V_C  - (q-En_C * n-V_En) / q-Ec_C))
    set n-R_Ec (n_R * (q-R_C - (q-En_C * n-R_En) / q-Ec_C))

    ; Yield coefficients
    set Y_EcF 0.7
    set Y_EnF 0.7
    set Y_VS  1.1
    set Y_HRD 1.1
    set Y_EcV ((n-V_Ec * q-Ec_C) / (n_V * q-V_C))
    set Y_EnV ((n-V_En * q-En_C) / (n_V * q-V_C))
    set Y_EcR ((n-R_Ec * q-Ec_C) / (n_R * q-R_C))
    set Y_EnR ((n-R_En * q-En_C) / (n_R * q-R_C))

    ; Ingestion/Assimilation parameters
    set i_max_CrlBrack 0.00461
    set F_h 0.00071
    set Food [resource-quantity] of patch-here
    set f (Food / (Food + F_h))

    ; Movement parameters
    set sigma movement-somatic-maintenance-cost-multiplier

    ; SU preference parameters
    set rho_S_N 0.1
    set rho_S_V 0.001
    set rho_D_N 0.1
    set rho_D_H 0.001

    ; Somatic growth parameters
    set M_V_SqBrack 0.0125523
    set k_G 1
    set M_V_max	4.365501
    set M_V_baby 0.1895243

    ; Reserve parameters
    set nu 0.08996
    set mdens_tot_max 6.183382
    set mdens_Ec_max (0.6 * mdens_tot_max)
    set mdens_En_max (0.4 * mdens_tot_max)
    set mdens_Ec_max_SqBrack (mdens_Ec_max * M_V_SqBrack)
    set mdens_En_max_SqBrack (mdens_En_max * M_V_SqBrack)
    set kappa 0.90653
    set kappaG_Ec 0.8
    set kappaR_Ec 0.8
    set kappaG_En 0.9
    set kappaR_En 0.9



    ; Maintenance parameters
    set k_M 0.06313964
    set k_J 0.002

    ; Development and reproduction parameters
    set k_R 1
    set M_H_p 11.76545
    set M_H_baby 	0.04925455
    set M_R_bufferbaby 1.71828

    set is-adult? 0
  ]

; SETUP HETEROGENOUS POPULATION TRAITS

  ; Determine the stage structure of population
  ask n-of (initial-population * initial-fraction-adults) turtles [
    set is-adult? 1
  ]


end


to initialize-DEB-variables
  ask turtles [
    ifelse is-adult? = 0 [
      ; Set state variables for juvenilles
      set M_V (M_V_baby + (random-float ((3 * M_V_baby) - M_V_baby)))              ; Size
      set M_H (M_H_baby + (random-float ((3 * M_H_baby) - M_H_baby)))               ; Development
      set M_R 0                  ; Reproduction buffer
    ][
      set M_V ((M_V_max / 3) + (random-float (M_V_max - (M_V_max / 3))))
      set M_H M_H_p    ; Development
      set M_R (0 + (random-float(M_R_bufferbaby - 0)))                ; Reproduction buffer
    ]

    set V (M_V / M_V_SqBrack)

    set scaled_Mdens_Ec initial-mean-carbon-reserve-filled-fraction + random-float (2 * initial-error-carbon-reserve-filled-fraction) - initial-error-carbon-reserve-filled-fraction

    if scaled_Mdens_Ec < 0.1 [
      set scaled_Mdens_Ec 0.1
    ]
    if scaled_Mdens_Ec > 1 [
      set scaled_Mdens_Ec 0.99
    ]

    set Mdens_Ec (scaled_Mdens_Ec * mdens_Ec_max)
    set M_Ec (Mdens_Ec * M_V)

    set scaled_Mdens_En initial-mean-carbon-reserve-filled-fraction + random-float (2 * initial-error-nitrogen-reserve-filled-fraction) - initial-error-nitrogen-reserve-filled-fraction
    if scaled_Mdens_En < 0.1 [
      set scaled_Mdens_En 0.1
    ]
    if scaled_Mdens_En > 1 [
      set scaled_Mdens_En 0.99
    ]

    set Mdens_En (scaled_Mdens_Ec * mdens_Ec_max)
    set M_En (Mdens_En * M_V)
  ]
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
  ]
end


;; ==================== TURTLE MOVEMENT PROCEDURES ====================;;
;;=====================================================================;;

to move
  ask turtles[
    if scaled_Mdens_En <= mismatch-threshold
    [
      ;; Seek high nitrogen food when body nitrogen (C:N) is above the threshold
      let best-patch max-one-of patches in-radius ((sqrt world-area) / 5) [resource-NC-ratio]  ;; Look for best nitrogen patch nearby
      face best-patch  ;; Face the patch with the highest nitrogen
      forward 1 ;; move forward
    ]
    if scaled_Mdens_En > mismatch-threshold
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




;;========================================================================================================================================================;;
;;===========================================================   GO/RUN SIMULATION PROCEDURES:  ===========================================================;;
;;===========================================================                                  ===========================================================;;
;;===========================================================             PATCHES              ===========================================================;;
;;========================================================================================================================================================;;

;;======================================================================================================================;;
;; ============================================ PATCH RESOURCE CALCULATIONS ============================================;;
;;======================================================================================================================;;

to update-patches

;;=== UPDATE RESOURCE QUANTITY AND QUALITY, UPDATE RESOURCE UPTAKE OF N FROM SOIL
  ask patches[
  ;; Establish values for parameter setting within different ifelse conditions
    set soil-N-crit 0
    set resource-quantity-min 0.01
    let min-kickstart-growth-Nmol (resource-quantity-min * (1 / target-max-CN-resource-ratio))
    let carryingcapacity-resource-Nmol (resource-carrying-capacity * (1 / target-min-CN-resource-ratio))
    let remaining-Nmol-to-capacity (carryingcapacity-resource-Nmol - resource-quantity-Nmol)

   ;; Calculate how much growth resources can do, and also set a hard check that this will not exceed carrying capacity due to rounding errors, etc
    set max-resource-grown (resource-quantity * resource-specific-growth-rate * (resource-carrying-capacity - resource-quantity) / resource-carrying-capacity)  ;Logistic growth equation
    if (max-resource-grown + resource-quantity) > resource-carrying-capacity[
      set max-resource-grown (resource-carrying-capacity - resource-quantity)
    ]

   ;; Based on above calculations, calculate the absolute amount of C and N moles that can be added to resource to maintain current resource C:N ratios
    set max-new-resource-Cmol (max-resource-grown)
    set max-new-resource-Nmol (resource-NC-ratio * max-new-resource-Cmol)

   ;; Determine how many resources are grown this tick, and how much C and N are added to the resource pool

      ; Check if resource-quantity is 0. If so, apply minimal growth to prevent being stuck at 0, so long as the soil N pool has at least that amount needed to maintain minimum C:N ratio of resoures
    (ifelse (resource-quantity = 0 and soil-N-quantity >= min-kickstart-growth-Nmol)[
      set new-resource-grown resource-quantity-min  ;; Small positive value to kickstart growth
      set resource-CN-ratio target-max-CN-resource-ratio
      set new-resource-quantity-Nmol min-kickstart-growth-Nmol
      set new-resource-quantity-Cmol (new-resource-quantity-Nmol * resource-CN-ratio)
      ]
      ; If resource-quantity is 0 and the soil N pool cannot establish a minimum amount of resources, no resources are grown and no C or N is added this tick
      (resource-quantity = 0 and soil-N-quantity < min-kickstart-growth-Nmol)[
        set new-resource-grown 0
        set new-resource-quantity-Nmol 0
        set new-resource-quantity-Cmol 0
      ]
      ; If resource quantity is at carrying capacity, and the soil N can supplement the resources to increase resource quality (to the minimum C:N ratio of resources specified), then N is added
      (resource-quantity = resource-carrying-capacity and soil-N-quantity >= remaining-Nmol-to-capacity)[
        set new-resource-grown 0
        set new-resource-quantity-Cmol 0
        set new-resource-quantity-Nmol remaining-Nmol-to-capacity
      ]
      ; If resource quantity is at carrying capacity, and there isn't enough N in the soil to supplement to the highest quality, the resources will simply drain the remaining N from soil to improve quality
      (resource-quantity = resource-carrying-capacity and soil-N-quantity < remaining-Nmol-to-capacity)[
        set new-resource-grown 0
        set new-resource-quantity-Cmol 0
        set new-resource-quantity-Nmol soil-N-quantity
      ]
      ; If resource quantity is below carrying capacity but there is not enough soil N to maintain the current C:N ratio of resources, the resources will still try to grow, and the quality of resources will decrease
      (resource-quantity < resource-carrying-capacity and soil-N-quantity < max-new-resource-Nmol)[
        set new-resource-quantity-Nmol soil-N-quantity
        set new-resource-grown max-resource-grown
        set new-resource-quantity-Cmol max-new-resource-Cmol

        ;; However, there is an absolute minimum quality resources must maintain (maximum specificied resource C:N ratio), and if soil N cannot supplement even to this amount, then the growth of resources is stunted such that the minimum C:N ratio is reached
        if ((max-new-resource-Cmol + resource-quantity-Cmol)/(new-resource-quantity-Nmol + resource-quantity-Cmol)) < target-max-CN-resource-ratio[
          set new-resource-quantity-Cmol (new-resource-quantity-Nmol * target-max-CN-resource-ratio)
          set new-resource-grown new-resource-quantity-Cmol
        ]
      ]
      ; If the resource is not at carry capacity, and the N in soil is abundant beyond the current projected growth need to maintain C:N ratio of resources, then they will grow and add C and N as projected above
      (resource-quantity < resource-carrying-capacity and soil-N-quantity >= max-new-resource-Nmol)[
        set new-resource-grown max-resource-grown
        set new-resource-quantity-Cmol max-new-resource-Cmol
        set new-resource-quantity-Nmol max-new-resource-Nmol
      ][
        ; There shouldn't be another condition..but in case, there is an error message thrown here
        write "something funky happened in the update-patches procedure, check end of ifelse() statement"
      ]
      )

    ;; Now finish calculating the new resource quantity, C and N in resources, and how much soil N is left
    set resource-quantity (new-resource-grown + resource-quantity)
    set resource-quantity-Cmol (new-resource-quantity-Cmol + resource-quantity-Cmol)
    set resource-quantity-Nmol (new-resource-quantity-Nmol + resource-quantity-Nmol)
    set soil-N-quantity (soil-N-quantity - new-resource-quantity-Nmol)

    ;; Re-calculate resource quality ratios
    if resource-quantity > 0 [
      set resource-CN-ratio (resource-quantity-Cmol / resource-quantity-Nmol)
      set resource-NC-ratio (1 / resource-CN-ratio)
    ]
  ]

end

to update-patch-misc
  ; Update informational monitors
  ifelse (num-resource-clusters = 0) [
    set realized-avg-world-resource-CN-ratio mean [resource-CN-ratio] of patches
    set realized-avg-world-resource-NC-ratio (1 / realized-avg-world-resource-CN-ratio)
  ]
  [
    set realized-avg-world-resource-CN-ratio mean [resource-CN-ratio] of patches
    set realized-avg-world-resource-NC-ratio (1 / realized-avg-world-resource-CN-ratio)

    set realized-noncluster-avg-resource-CN-ratio mean [resource-CN-ratio] of patches with [cluster-patch? = 0]
    set realized-noncluster-avg-resource-NC-ratio (1 / realized-noncluster-avg-resource-CN-ratio)

    set realized-cluster-avg-resource-CN-ratio mean [resource-CN-ratio] of patches with [cluster-patch? = 1]
    set realized-cluster-avg-resource-NC-ratio (1 / realized-cluster-avg-resource-CN-ratio)
    ]
  ; Update patch coloring
  color-patches

end

;;======================== DEB MODEL PROCEDURES ======================;;
;;====================================================================;;

;; ========== Calculate fluxes for tick ==========;;
;;================================================;;

to go-DEB-fluxes

  go-DEB-flux-ingest-assimilate-egest
  go-DEB-flux-rootsolve-fluxes
  go-DEB-flux-excrete


end

to go-DEB-flux-ingest-assimilate-egest

  ask turtles [
  ; Ingestion
  set J_FA ((f * i_max_CrlBrack * ((M_V / M_V_SqBrack) ^ (2 / 3))) * q-F_C)

;   Ensure J_FA does not exceed available resource-quantity
  if J_FA > [resource-quantity] of patch-here [
    set J_FA ([resource-quantity] of patch-here)
  ]

  ;== Remove eaten food from resource quantity on patch and remove the N and C total resource moles thus far eaten
  ask patch-here [
    let max-J_FA [J_FA] of myself  ; Capture J_FA value from the turtle
    ; Ensure resource-quantity never goes below zero
    if max-J_FA > resource-quantity [
      set max-J_FA (resource-quantity)
    ]

    ; Update the resource-quantity and mole values based on the actual amount eaten
    set resource-quantity (resource-quantity - max-J_FA)
    let Cmol-eaten-this-tick  max-J_FA
    let Nmol-eaten-this-tick (max-J_FA * resource-NC-ratio)

    set resource-quantity-Cmol (resource-quantity-Cmol - Cmol-eaten-this-tick)
    set resource-quantity-Nmol (resource-quantity-Nmol - Nmol-eaten-this-tick)
  ]

    ; Assimilation and egestion
    set q-F_C [resource-CN-ratio] of patch-here
    set q-F_N 1

    set J_AEn (Y_EnF  * J_FA * q-En_C * (q-F_N / q-F_C))
    set J_AEc ((Y_EcF  * J_FA) - J_AEn)
    set J_Wc ((1 - Y_EcF) * J_FA)
    set J_Wn ((1 - Y_EnF) * J_FA * (q-F_N / q-F_C))
  ]
end

to go-DEB-flux-rootsolve-fluxes
 ask turtles[
  ;;set parameters
    let parList []
    set parList lput M_V parList
    set parList lput M_H parList
    set parList lput M_Ec parList
    set parList lput M_En parList

    set parList lput k_G parList
    set parList lput M_V_SqBrack parList
    set parList lput k_M parList
    set parList lput kappa parList
    set parList lput sigma parList
    set parList lput nu parList
    set parList lput rho_S_N parList
    set parList lput rho_S_V parList
    set parList lput Y_EcV parList
    set parList lput Y_EnV parList

    set parList lput k_R parlist
    set parList lput rho_D_N parList
    set parList lput rho_D_H parList
    set parList lput Y_EcR parList
    set parList lput Y_EnR parList
    set parList lput k_J parList

    (sr:set "parList" parList)

    sr:run "x <- flux_calculations(parList)"
;  print sr:runresult "x"
    set J_EcC sr:runresult "x[1]"
    set J_EnC sr:runresult "x[2]"
    set J_GV sr:runresult "x[3]"
    set r_dot sr:runresult "x[4]"
    set J_EcS sr:runresult "x[5]"
    set J_EnS sr:runresult "x[6]"
    set J_VS sr:runresult "x[7]"
    set J_EcG sr:runresult "x[8]"
    set J_EnG sr:runresult "x[9]"
    set J_EcG_plus sr:runresult "x[10]"
    set J_EnG_plus sr:runresult "x[11]"
    set J_EcG_neg sr:runresult "x[12]"
    set J_EnG_neg sr:runresult "x[13]"
    set J_EcD sr:runresult "x[14]"
    set J_EnD sr:runresult "x[15]"
    set J_HD sr:runresult "x[16]"
    set J_EcR sr:runresult "x[17]"
    set J_EnR sr:runresult "x[18]"
    set J_RR sr:runresult "x[19]"
    set J_EcR_plus sr:runresult "x[20]"
    set J_EnR_plus sr:runresult "x[21]"
    set J_EcR_neg sr:runresult "x[22]"
    set J_EnR_neg sr:runresult "x[23]"
  ]

end

to go-DEB-flux-excrete
 ask turtles[
  ; Calculate the recycling fluxes to reserves from the growth and reproduction SUs
    set J_GREc ((J_EcG_plus * kappaG_Ec) + (J_EcR_plus * kappaR_Ec))
    set J_GREn ((J_EnG_plus * kappaG_En) + (J_EnR_plus * kappaR_En))

 ; Calculate the excretion from rejected fluxes from the growth and reproduction SUs
    set J_GRX_Ec ((J_EcG_plus * (1 - kappaG_Ec)) + (J_EcR_plus * (1 - kappaR_Ec)))
    set J_GRX_Ec ((J_EnG_plus * (1 - kappaG_En)) + (J_EnR_plus * (1 - kappaR_En)))
  ]
end

;; ========== Calculate change in state variables ==========;;
;;==========================================================;;

to  go-DEB-ODE
  ask turtles[
    ; Differential Equations/ Mass balances
    set dWc-dt (J_Wc)
    set dWn-dt (J_Wn)
    set dM_Ec-dt (J_AEc + J_GREc - J_EcC)
    set dM_En-dt (J_AEn + J_GREn - J_EnC)
    set dM_V-dt (J_GV)
    ifelse (M_H < M_H_p) [
      set dM_H-dt (J_RR - J_HD)
      set dM_R-dt 0
    ][
      set dM_H-dt 0
      set dM_R-dt J_RR
    ]
    ifelse (M_H < M_H_p) [
        set dXc-dt J_EcS + J_EnS + J_EcD + J_EnD + J_GRX_Ec + J_GRX_En + J_RR
    ][
        set dXc-dt J_EcS + J_EnS + J_EcD + J_EnD + J_GRX_Ec + J_GRX_En
    ]
    ifelse (M_H < M_H_p) [
        set dXn-dt (((q-En_N / q-En_C) * (J_EnS + J_EnD + J_GRX_En)) + ((q-R_N / q-R_C) * (J_RR)))
    ][
        set dXn-dt  ((q-En_N / q-En_C) * (J_EnS + J_EnD + J_GRX_En))
    ]
  ]
end

;; ========== Update turtle model values  ==========;;
;;===============================================;;

to update-turtle-variables

  ;; Update DEB model variables
  ask turtles[
    set M_V (M_V + dM_V-dt)
    set V (M_V / M_V_SqBrack)
    set M_Ec (M_Ec + dM_Ec-dt)
    set M_En (M_En + dM_En-dt)
    set Mdens_Ec (M_Ec / M_V)
    set Mdens_En (M_En / M_V)
    set M_H (M_H + dM_H-dt)
    set M_R (M_R + dM_R-dt)
  ]
end

;; ========== Reproduction submodel  ==========;;
;;=============================================;;

to have-baby

  ask turtles [
    if M_R >= M_R_bufferbaby [
      hatch 1 [
        ; Set important parameters
        set is-adult? 0
        set did-baby? 0
        set low-Nreserve-threshold 0.2
        set is-Nreserve-low? 0
        set starved-to-death? 0
        set demographic-stochasticity-death? 0
        set should-die? 0

        ; Below, same parameters that original consumer turtles intialized with
            ; Stoichiometry parameters
        set q-F_C [resource-CN-ratio] of patch-here
        set q-F_N 1
        set q-Ec_C 1
        set q-Ec_N 0
        set q-En_C 4
        set q-En_N 1
        set q-V_C consumer-body-CN-ratio
        set q-V_N 1
        set q-R_C q-V_C
        set q-R_N 1

        ; Moles of state variables to calculate yields
        set n_V 1
        set n_R 1
        set n-V_En 1
        set n-R_En 1
        set n-V_Ec  (n_V * (q-V_C  - (q-En_C * n-V_En) / q-Ec_C))
        set n-R_Ec (n_R * (q-R_C - (q-En_C * n-R_En) / q-Ec_C))

        ; Yield coefficients
        set Y_EcF 0.7
        set Y_EnF 0.7
        set Y_EcV ((n-V_Ec * q-Ec_C) / (n_V * q-V_C))
        set Y_EnV ((n-V_En * q-En_C) / (n_V * q-V_C))
        set Y_EcR ((n-R_Ec * q-Ec_C) / (n_R * q-R_C))
        set Y_EnR ((n-R_En * q-En_C) / (n_R * q-R_C))

        ; Ingestion/Assimilation parameters
        set i_max_CrlBrack 0.00461
        set F_h 0.00071
        set Food [resource-quantity] of patch-here
        set f (Food / (Food + F_h))

        ; Movement parameters
        set sigma movement-somatic-maintenance-cost-multiplier

        ; SU preference parameters
        set rho_S_N 0.1
        set rho_S_V 0.001
        set rho_D_N 0.1
        set rho_D_H 0.001

        ; Somatic growth parameters
        set M_V_SqBrack 0.0125523
        set k_G 1
        set M_V_max	4.365501
        set M_V_baby 0.1895243

        ; Reserve parameters
        set nu 0.08996
        set mdens_tot_max 6.183382
        set mdens_Ec_max (0.6 * mdens_tot_max)
        set mdens_En_max (0.4 * mdens_tot_max)
        set mdens_Ec_max_SqBrack (mdens_Ec_max * M_V_SqBrack)
        set mdens_En_max_SqBrack (mdens_En_max * M_V_SqBrack)
        set kappa 0.90653
        set kappaG_Ec 0.8
        set kappaR_Ec 0.8
        set kappaG_En 0.9
        set kappaR_En 0.9

        ; Maintenance parameters
        set k_M 0.06313964
        set k_J 0.002

        ; Development and reproduction parameters
        set k_R 1
        set M_H_p 11.76545
        set M_H_baby 	0.04925455
        set M_R_bufferbaby 1.71828

        set is-adult? 0

        ; Set state variables
        set M_V M_V_baby
        set V (M_V / M_V_SqBrack)
        set Mdens_Ec [Mdens_Ec] of myself
        set Mdens_En [Mdens_En] of myself
        set M_Ec (Mdens_Ec * M_V)
        set M_En (Mdens_En * M_V)
        set M_H M_H_baby
        set M_R 0
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

to reset-ODE-and-fluxes
  ask turtles[
    ; Reset all ODEs to 0 for the next tick
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

    ; Reset all fluxes to 0 for the next tick
    set J_FA 0
    set J_AEc 0
    set J_AEn 0
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
    set J_GV 0
    set J_VS 0
    set J_EcG_plus 0
    set J_EnG_plus 0
    set J_EcG_neg 0
    set J_EnG_neg 0
    set J_EcR 0
    set J_EnR 0
    set J_RR 0
    set J_HD 0
    set J_EcR_plus 0
    set J_EnR_plus 0
    set J_EcR_neg 0
    set J_EnR_neg 0
    set J_GREc 0
    set J_GREn 0
    set J_GRX_Ec 0
    set J_GRX_En 0
  ]
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
end

to plot-total-N
;; Calculate the sum of the patch variable
  set sum-resource-N-total sum [resource-quantity-Nmol] of patches
  set sum-consumer-N-total (sum [ ((q-V_N / q-V_C) * M_V) + ((q-En_N / q-En_C) * M_En) + ((q-R_N / q-R_C) * (M_H + M_R)) +  temp-N-eaten ] of turtles)
  set sum-soil-N-total (sum [soil-N-quantity] of patches)
  set sum-world-N-total (sum-resource-N-total + sum-consumer-N-total + sum-soil-N-total)
  set avg-resource-quantity mean [resource-quantity] of patches

  set-current-plot "Total N Dynamics"
end

to plot-resource-densities
  set total-resource-quantity sum [resource-quantity] of patches

  set-current-plot "Total Resources"
end

to plot-consumer-densities
  set total-consumers count turtles
  set total-adult-consumers count turtles with [is-adult? = 1]
  set total-juvenille-consumers count turtles with [is-adult? = 0]

  set-current-plot "Total Consumers"
end
@#$#@#$#@
GRAPHICS-WINDOW
12
10
633
632
-1
-1
6.07
1
10
1
1
1
0
1
1
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
976
257
1198
290
initial-population
initial-population
0
100
30.0
1
1
NIL
HORIZONTAL

BUTTON
654
52
718
86
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
642
102
797
162
consumer-body-CN-ratio
7.0
1
0
Number

SLIDER
976
330
1198
363
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
720
52
784
86
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
649
439
839
472
overall-resource-patchiness
overall-resource-patchiness
0
0.5
0.25
0.05
1
NIL
HORIZONTAL

SLIDER
976
601
1200
634
mismatch-threshold
mismatch-threshold
0.1
0.5
0.5
0.1
1
NIL
HORIZONTAL

CHOOSER
976
553
1200
598
movement-type
movement-type
"brownian-walk" "correlated-random-walk"
1

SLIDER
976
637
1200
670
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
976
293
1198
326
initial-fraction-adults
initial-fraction-adults
0
1
0.8
0.1
1
NIL
HORIZONTAL

SLIDER
648
331
947
364
resource-specific-growth-rate
resource-specific-growth-rate
0.001
0.002
0.001
0.005
1
NIL
HORIZONTAL

SLIDER
648
293
947
326
resource-carrying-capacity
resource-carrying-capacity
1
20
15.0
2
1
NIL
HORIZONTAL

MONITOR
1033
42
1270
87
Consumer Body N:C Ratio
consumer-body-NC-ratio
3
1
11

TEXTBOX
1200
18
1345
37
Informational Calculations
12
0.0
1

TEXTBOX
680
378
924
404
Additional Resource Quality Controls
12
0.0
1

SLIDER
648
257
946
290
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
976
517
1200
550
mortality-risk
mortality-risk
0
0.01
0.0056
0.0001
1
NIL
HORIZONTAL

PLOT
1217
187
1621
522
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
"World Total N" 1.0 0 -16777216 true "" "plot sum-world-N-total"
"Resource N" 1.0 0 -13840069 true "" "plot sum-resource-N-total"
"Consumer N" 1.0 0 -5825686 true "" "plot sum-consumer-N-total"
"Soil N" 1.0 0 -6459832 true "" "plot sum-soil-N-total"

INPUTBOX
805
160
982
220
target-min-CN-resource-ratio
15.0
1
0
Number

INPUTBOX
805
101
982
161
target-max-CN-resource-ratio
50.0
1
0
Number

INPUTBOX
805
41
982
101
target-world-avg-resource-CN-ratio
30.0
1
0
Number

SLIDER
645
581
836
614
patch-soil-N-multiplier
patch-soil-N-multiplier
0.1
2
0.6
0.1
1
NIL
HORIZONTAL

MONITOR
1274
42
1511
83
Actual World  Resource CN Ratio Mean
realized-avg-world-resource-CN-ratio
2
1
10

MONITOR
1274
130
1511
171
Actual High-Quality Resource Cluster CN Ratio Mean
realized-cluster-avg-resource-CN-ratio
2
1
10

MONITOR
1274
86
1511
127
Actual Low-Quality  Resource CN Ratio Mean
realized-noncluster-avg-resource-CN-ratio
2
1
10

SLIDER
649
402
839
435
num-resource-clusters
num-resource-clusters
0
50
40.0
1
1
NIL
HORIZONTAL

MONITOR
1033
89
1270
130
% of world patches in high-quality clusters
count patches with [cluster-patch? = 1] / world-area
2
1
10

SLIDER
649
475
839
508
resource-cluster-scatter
resource-cluster-scatter
0
10
0.5
0.5
1
NIL
HORIZONTAL

PLOT
1422
523
1622
673
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
1216
523
1416
673
Total Resources
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
"default" 1.0 0 -14439633 true "" "plot total-resource-quantity"

TEXTBOX
689
16
984
51
Required User Input/Simulation Controls
14
0.0
1

TEXTBOX
698
237
913
267
Additional Resource Growth Controls
12
0.0
1

TEXTBOX
1019
237
1169
255
Additional Turtle Controls
12
0.0
1

TEXTBOX
733
525
859
550
Additional Soil Controls
12
0.0
1

SLIDER
976
367
1199
400
initial-mean-carbon-reserve-filled-fraction
initial-mean-carbon-reserve-filled-fraction
0.3
1
0.8
0.1
1
NIL
HORIZONTAL

SLIDER
976
404
1200
437
initial-error-carbon-reserve-filled-fraction
initial-error-carbon-reserve-filled-fraction
0
0.4
0.05
0.05
1
NIL
HORIZONTAL

SLIDER
976
442
1201
475
initial-mean-nitrogen-reserve-filled-fraction
initial-mean-nitrogen-reserve-filled-fraction
0.3
1
0.7
0.1
1
NIL
HORIZONTAL

SLIDER
976
479
1201
512
initial-error-nitrogen-reserve-filled-fraction
initial-error-nitrogen-reserve-filled-fraction
0
0.3
0.05
0.05
1
NIL
HORIZONTAL

MONITOR
12
645
174
690
Mean scaled C-reserve density
mean [scaled_Mdens_Ec] of turtles
2
1
11

MONITOR
176
646
379
691
StDev of scaled C-reserve density
(standard-deviation [scaled_Mdens_Ec] of turtles)
2
1
11

MONITOR
355
646
514
691
Mean scaled N-reserve density
mean [scaled_Mdens_En] of turtles
2
1
11

MONITOR
515
646
718
691
StDev of scaled N-reserve density
standard-deviation [scaled_Mdens_En] of turtles
2
1
11

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
