(* Mathematica Package *)

(* Created by the Wolfram Workbench 11.05.2010 *)

BeginPackage["EscapeProbability`"]
(* Exported symbols added here with SymbolName::usage *) 
Options[iterate] = {Geometry -> "LVG"};
Geometry::usage="Is an options whose setting specifies which escape \
probability approximation will be applied. Possible choices are:\ 
\"Sphere\", \"LVG\", \"PlaneParallel\", \"Turbulent\".";
escapeProbability::optx = 
  "`1` is an unknown geometry option. Switching to default \"LVG\".";
EscapeProbabilityRun::conflopts = 
  "Conflicting options. `1` and `2` may not both be True. Setting `3` \
to False.";
EscapeProbabilityRun::denxmatch =
  "Too few densities of collisional partners provided, Densities for\
  `1` collision partners are needed. Replicating the provided\
  density(ies) to fill up the list. WARNING! Results is likely to be\
  unphysical.";
EscapeProbabilityRun::usage = 
  "EscapeProbabilityRun[LAMBDAfile, N, FWHM, n, T] performs a \
radiative transfer calculation using the escape probability \
approximation. The input parameters are the column density of the \
required species N in \!\(\*SuperscriptBox[\"cm\", 
RowBox[{\"-\", \"2\"}]]\), the FWHM line width in km/sec, the total \
hydrogen density in \!\(\*SuperscriptBox[\"cm\", 
RowBox[{\"-\", \"3\"}]]\), and the kinetic gas temperature T in K. \
The atomic and molecular data is taken from the LAMBDA Leiden Atomic \
and Molecular Database.
  ";
ShowCollisionPartners::combop = 
  "If ortho- and para-H2 is given seperately as collision partners \
  the default behavior is to combine both into a single collision \
partner. In that case only one \
  H2 density needs to be given.";
ShowCollisionPartners::usage =
  "ShowCollisionPartners[LAMBDA file] displays collisional partners\
  for which collision rates are given in the LAMBDA database files.";
ShowCollisionTemperatures::usage =
  "ShowCollisionTemperatures[LAMBDA file] displays for each collisional\
  partner the temperatures for which collision rates are listed in the \
  LAMBDA database files.";
ShowCollisionRates::usage =
  "ShowCollisionRates[LAMBDA file] displays for each collisional\
  partner the collision rates which are listed in the \
  LAMBDA database files.";
ShowEnergyLevels::usage =
  "ShowEnergyLevels[LAMBDA file] displays which energy levels are \
  listed in the LAMBDA database files.";
ShowTransitions::usage =
  "ShowTransitions[LAMBDA file] displays which levels transitions are \
  listed in the LAMBDA database files.";
Options[ShowCollisionRates] = {CombineOrthoPara -> True};
Options[ShowCollisionTemperatures] = {CombineOrthoPara -> True};
Options[lambdaInput2] = {CombineOrthoPara -> True};
CombineOrthoPara::usage="CombineOrthoPara is an option which defines \
how to handle seperate collision rates for ortho- and para-H2. The \
default behavior is to assume LTE and combine both into a single rate."; 
Options[EscapeProbabilityRun] = {
    Geometry -> "Sphere", 
    MaximumIterations -> 10000, 
    MinimumLevelPopulation -> 1. 10^-30,
    FullOutput -> True, 
    ShowLevelPopulation -> False, 
    ShowAntennaTemperature -> False, 
    ShowExcitationTemperature -> False,
    ShowEscapeProbability -> False,
    ShowOpticalDepth -> False, 
    FrequencyRange -> {0, 2000}, 
    TableOutput -> True, 
    BackgroundTemperature -> 2.73, 
    InitialPopulationTemperature -> 2.73,
    SameTest->Automatic,
    CombineOrthoPara -> True};  
MaximumIterations::usage="MaximumIterations is an option that specifies\
 the maximum number of iterations that should be tried in the escape probability\ 
 calculations. The default is 10000.";
MinimumLevelPopulation::usage="MinimumLevelPopulation is an option that specifies\
 the lowest level population. Default is 10^-30.";
FullOutput::usage="FullOutput is an option that specifies what output should be\ 
produced by EscapeProbabilityRun. FullOutput->True  is the default.";
ShowLevelPopulation::usage="ShowLevelPopulation is an option for EscapeProbabilityRun\ 
 that will set the output to the level population numbers.";
ShowAntennaTemperature::usage="ShowAntennaTemperature is an option for EscapeProbabilityRun.\
  It will give the antenna temperatures as output.";
ShowExcitationTemperature::usage="ShowExcitationTemperature is an option for EscapeProbabilityRun.\
  It will give the excitation temperatures as output.";
ShowOpticalDepth::usage="ShowOpticalDepth is an option for EscapeProbabilityRun.\
  It will give the optical depths as output.";
ShowEscapeProbability::usage="ShowEscapeProbability is an option for EscapeProbabilityRun.\
  It will give the final escape probabilities as output.";
FrequencyRange::usage="FrequencyRange is an option for EscapeProbabilityRun which specifies\ 
  which frequency range should be included in the output.";
TableOutput::usage="TableOutput is an option for EscapeProbabilityRun which specifies\ 
  whether the output should be as formated table/grid.";  
BackgroundTemperature::usage="BackgroundTemperature is an option for EscapeProbabilityRun which specifies\ 
  the temperature of the cosmic microwave background. Default is 2.73 K.";   
InitialPopulationTemperature::usage="InitialPopulationTemperature is an option for EscapeProbabilityRun which specifies\ 
  which LTE temperature should be used for the initial level population. Default is 2.73 K.";
 
Begin["`Private`"]
(* Implementation of the package *)
 (* function to iterate in FixedPoint. Takes level population as input and output *)
  iterate[besetzungsvektor_, opts : OptionsPattern[]] :=
      Block[ {coeffMatrix, betavektor, beta, normierteBesetzungszahlen, 
         bnuTex, exr, radiativeRateMatrix,tau,geom,tex},
          geom = Geometry/. {opts} /. Options[EscapeProbabilityRun];
        (*normierteBesetzungszahlen=Normalize[besetzungsvektor];*)
          normierteBesetzungszahlen = besetzungsvektor;
          Quiet[tex = calcTexCompiled[besetzungsvektor, levels[[All, 3]], transitions[[All, 3]], 
          transitions[[All, 2]], levelenergiesKelvin,  minpop];
                {CompiledFunction::cfn}];
          (* calculate source function *)
(*        bnuTex = calcSource[besetzungsvektor];*)
          bnuTex = calcSourceCompiled[tex,xnu,xnut];
          (*bnuTex = calcSourceUncompiled[tex];*)
          (* calculate optical depths per transition *)
          tau = column/\[CapitalDelta]v tauListCompiled[besetzungsvektor, transitions[[All, 
             4]], levels[[All, 3]],transitions[[All, 3]], transitions[[All, 2]],xnut];
          (* calculate escape probability from the tau's *)
          betavektor = escapeProbability[geom][#]&/@tau;
(*          beta[n_Integer] :=
              betavektor[[n]];*)
          (* replaced xnu^3 with xnut *)
          (* calculate local radiation field available for level popukation *)
          exr = (background betavektor + (1. - betavektor) bnuTex)/(thc xnut);
          (* populate coefficient matrix of the rate equation with radiative rates *)
          radiativeRateMatrix = SparseArray[populateRadiativeMatrix[exr]];
		  (* inrement iteration counter *)
          iterCount++;
         (* 
            combine collisional and radiative rates. replace last matrix row with
            conservation equation for population numbers, i.e. sum of all populations
            equals 1 
          *)
         (* density multiplied to collision matrix outside of iterate now *)
          coeffMatrix = 
           ReplacePart[-collisionMatrix + radiativeRateMatrix, 
            Length@levels -> ConstantArray[1., Length@levels]];
          (*
            solve linear system and replace all level populations smaller than minpop
            with minpop
          *)
          Max[#, minpop] & /@ 
           LinearSolve[N[coeffMatrix, 30], 
            N[Append[ConstantArray[0., Length@levels - 1], 1.], 30], 
            FilterRules[{opts}, Options[LinearSolve]](*,Method->"Krylov"*)]
      ];  
    
   
  populateRadiativeMatrix[exr_] :=
      Block[ {(*levels,transitions,colltemps,collrates,wavenumber,hnu,xnu,*)
        ruleUpUp, ruleLowLow, ruleUpLow, ruleLowUp, allrules(*,gstat*)},
       (*{levels,transitions,colltemps,collrates}=onionFileInput;
       gstat=levels[[All,3]];*)
(*          ruleUpUp = 
           MapIndexed[
            Rule[{#1[[2]], #1[[2]]}, #[[4]] (1.0 + exr[[First@#2]])] &, 
            transitions];*)
          (* using Thread instead of MapIndexed is ca. 4 times faster *)
          ruleUpUp =   
           Thread[transitions[[All, {2, 2}]] -> 
           transitions[[All, 4]] (1.0 + exr)];
          ruleLowLow =
           Thread[transitions[[All, {3, 3}]] -> tabLowLow exr];
(*          ruleLowLow = 
           MapIndexed[
            Rule[{#1[[3]], #1[[3]]}, #[[4]] gstat[[#1[[2]]]] exr[[First@#2]]/
                gstat[[#1[[3]]]]] &, transitions];*)
          ruleUpLow = 
          Thread[transitions[[All, {2, 3}]] -> -tabLowLow exr];
(*          ruleUpLow = 
           MapIndexed[
            Rule[{#1[[2]], #1[[3]]}, -#[[4]] gstat[[#1[[2]]]] exr[[First@#2]]/
                gstat[[#1[[3]]]]] &, transitions];*)
          ruleLowUp = 
           Thread[transitions[[All, {3, 2}]] -> 
           -transitions[[All, 4]] (1.0 + exr)];
          (*ruleLowUp = 
           MapIndexed[
            Rule[{#1[[3]], #1[[2]]}, -#[[4]] (1.0 + exr[[First@#2]])] &, 
            transitions];*)
          (* combine all rules into a single list *)
          allrules = 
          Flatten[{ruleUpUp, ruleLowLow, ruleUpLow, 
           ruleLowUp}];
(*          allrules = 
           Join[ruleUpUp, ruleLowLow, ruleUpLow, 
            ruleLowUp];(*//.sumSameElements*)*)
          (*Dispatch@*)
          (* sum up all rates that belong to the same {i,j}->... *)
          (* Gather is quite slow, but I haven't found a better solution yet *)
          Rule[First[#][[1]], #[[All, 2]] // Total] & /@ 
          Gather[allrules, (#1[[1]] === #2[[1]] )&]
      ]
totalBackground[T_?NumberQ] :=
    Block[ {hnu},
    	(* calculate the planck function for a given temperature for 
    	all transition frequencies *)
        hnu = fk xnu/T; (* photon energy *)
        Table[If[ hnu[[i]] > 160.,
                  0.,
                  (* changed xnu[[i]]^3 to xnut[[i]]*)
                  (thc  xnut[[i]])/(
                  Exp[fk xnu[[i]]/T] - 1.0)
              ], {i, 1, Length@transitions}]
    ]
calcSource[besetzung_] :=
    Module[ {tex},
        (* calculate the source function, i.e. the planck function for the
        LTE temperature equivalent of the present population ratio of the
        upper and lower level of each transition  *)
        (* this is an old version that takes populations as input and 
        calculate Tex itself *)
        tex = calcTex[besetzung];
        (* changed xnu[[i]]^3 to xnut[[i]]*)
        (* changed h cl/kB to fk*)
        Table[If[ fk xnu[[i]] > 160.,
                  0.,
                  (thc  xnut[[i]])/(
                  Exp[fk xnu[[i]]/tex[[i]]] - 1.0)
              ], {i, 1, Length@transitions}]
    ]
    
calcSourceCompiled = 
 Compile[{{tex, _Real, 1}, {xnu, _Real, 1}, {xnut, _Real, 1}},
 	(* compiled version of the source function calculation*)
 	(* 
 	sometimes we have an error indicating that the input is not a
 	proper real. Haven't found the reason yet. Might also be in 
 	the excitation temp routine. 
 	*)
  Module[ {thc, fk},
      thc = 3.9728909742875315`*^-16;
      fk = 1.438775150569446`;
      Table[If[ fk xnu[[i]] > 160.,
                0.,
                (thc xnut[[i]])/(Exp[fk xnu[[i]]/tex[[i]]] - 1.0)
            ], {i, 1, 
        Length@xnu}]
  ]]
  
(* calcSourceUncompiled[tex_] :=
     Module[ {},
         Table[If[ fk xnu[[i]] > 160.,
                   0.,
                   (thc xnut[[i]])/(Exp[fk xnu[[i]]/tex[[i]]] - 1.0)
               ], {i, 1, 
           Length@xnu}]
     ]*)

calcTexCompiled = 
 Compile[{{besetzung, _Real, 1}, 
           {statisticalweights, _Real,1}, 
           {lowerlevelindex, _Integer, 1}, 
           {upperlevelindex, _Integer, 1}, 
           {levelenergiesKelvin, _Real, 1}, 
           {minpop, _Real}},
           (* calculate excitation temperature, i.e. the LTE temperature
           that produces the given population ratio of upper to lower level
           for each transition *)
           (* Has to be compiled because it will be executed very often.*)
            Table[
                If[ besetzung[[upperlevelindex[[i]]]] == 
   				  	besetzung[[lowerlevelindex[[i]]]] == minpop,
   				  	(* don't know if this makes sense.*)
                    (-(levelenergiesKelvin[[upperlevelindex[[i]]]] - 
                         levelenergiesKelvin[[lowerlevelindex[[i]]]]))/
                     Log[statisticalweights[[lowerlevelindex[[i]]]]/
                       statisticalweights[[upperlevelindex[[i]]]] (
                       besetzung[[upperlevelindex[[i]]]] + minpop/10.)/
                       besetzung[[lowerlevelindex[[i]]]]],
                           \
                    (-(levelenergiesKelvin[[upperlevelindex[[i]]]] - 
                    levelenergiesKelvin[[lowerlevelindex[[i]]]]))/
                    Log[statisticalweights[[lowerlevelindex[[i]]]]/
                    statisticalweights[[upperlevelindex[[i]]]] \
					besetzung[[upperlevelindex[[i]]]]/
                    besetzung[[lowerlevelindex[[i]]]]]
               	 	], 
               	 {i, 1, Length@lowerlevelindex}]
              ]
 
calcTex[besetzung_] :=
    Module[ {},
  (* uncompiled version*)
        Table[If[ 	besetzung[[transitions[[i, 2]]]] == 
          			besetzung[[transitions[[i, 3]]]] == minpop,
                  (-(levelenergiesKelvin[[transitions[[i, 2]]]] - 
                       levelenergiesKelvin[[transitions[[i, 3]]]]))/
                   Log[levels[[transitions[[i, 3]], 3]]/
                     levels[[transitions[[i, 2]], 3]] (
                     besetzung[[transitions[[i, 2]]]] + minpop/10.)/
                     besetzung[[transitions[[i, 3]]]]],
                  (-(levelenergiesKelvin[[transitions[[i, 2]]]] - 
                  levelenergiesKelvin[[transitions[[i, 3]]]]))/
                  Log[levels[[transitions[[i, 3]], 3]]/
                  levels[[transitions[[i, 2]], 3]] besetzung[[
                  transitions[[i, 2]]]]/besetzung[[transitions[[i, 3]]]]]
              	], 
              {i, 1, Length@transitions}]
    ]

 
tauListe[besetzungszahlen_] :=
    Module[ {},
    	(* calculate all optical depths *)
        Table[transitions[[i, 
             4]]/(8. N[\[Pi]] wavenumber[
              transitions[[i, {2, 3}]]]^3) (levels[[transitions[[i, 2]], 
              3]]/
             levels[[transitions[[i, 3]], 
              3]] besetzungszahlen[[transitions[[i, 3]]]] - 
            besetzungszahlen[[transitions[[i, 2]]]]), {i, 1, 
          Length@transitions}]
    ]



(*tauListe3[besetzung_,avalues_,statisticalweights_,lowerlevelindex_,upperlevelindex_,xnut_]:=
       Block[{}, Table[
            avalues[[i]]/(8. N[\[Pi]] xnut[[i]]) (statisticalweights[[upperlevelindex[[i]]]]/
             statisticalweights[[lowerlevelindex[[i]]]] besetzung[[lowerlevelindex[[i]]]] - 
            besetzung[[upperlevelindex[[i]]]]), 
                {i, 1, Length@avalues}]]*)


tauListCompiled = Compile[
    {
    {besetzung, _Real, 1}, 
    {avalues, _Real, 1},
    {statisticalweights, _Real, 1}, 
    {lowerlevelindex, _Integer, 1}, 
    {upperlevelindex, _Integer, 1}, 
    {xnut, _Real, 1}
    },
    (* compiled version because it gets executed a lot *)
        Table[
            avalues[[i]]/(8. N[\[Pi]] xnut[[i]]) (statisticalweights[[upperlevelindex[[i]]]]/
             statisticalweights[[lowerlevelindex[[i]]]] besetzung[[lowerlevelindex[[i]]]] - 
            besetzung[[upperlevelindex[[i]]]]), 
                {i, 1, Length@avalues}]
]

populateLTE[T_?NumberQ] :=
    Module[ {partitionFunction},
    	(* populate the energy levels according to a given LTE temperature *)
        (* calculate the partition function *)
        partitionFunction = 
         Sum[levels[[i, 3]]*Exp[-(levelenergiesKelvin[[i]])/T], {i, 1, 
           Length@levels}];
        (* calculate the actual LTE level populations *)
        Table[levels[[i, 3]]*
          Exp[-levelenergiesKelvin[[i]]/T]/partitionFunction, {i, 1, 
          Length@levels}]
    ]
 
SetUpCollisionMatrix[temperature_?NumberQ] :=
    Module[ {T = temperature,deltaEs, 
      statweightratio, interpolateCollRates, collratesTemp, 
      upwardCollisionRules, upCollratesTemp, downwardCollisionRules, 
       tmpMatrix, diagonalRules},
      (* set up the collision matrix *)
      (* energy difference between upper and lower level in Kelvin *)
        deltaEs = (levelenergiesKelvin[[#[[1]]]] - 
             levelenergiesKelvin[[#[[2]]]]) & /@ collrates[[All, 1 ;; 2]];
       (* ratio of the statistical weights to compute the upward rates 
        from statistical equilibrium 
        *)
        statweightratio = (levels[[#[[1]], 3]]/levels[[#[[2]], 3]] &) /@ 
          collrates[[All, 1 ;; 2]];
        (* interpolate the provided collision rates to the given temperature *)
        interpolateCollRates[temp_, list_, colltemps_] :=
            ListInterpolation[list, {colltemps}][temp];
        collratesTemp = 
         Map[interpolateCollRates[T, #[[3 ;; Length@colltemps + 2]], 
            colltemps] &, collrates];
        (* calculate the upward rates *)
        upCollratesTemp = collratesTemp statweightratio*Exp[-deltaEs/T];
        (* convert rates into rules to populate the sparse array *)
        upwardCollisionRules = 
         Rule[#[[1]], #[[2]]] & /@ 
          Transpose[{collrates[[All, 1 ;; 2]], upCollratesTemp}];
        downwardCollisionRules = 
         Rule[Reverse@#[[1]], #[[2]]] & /@ 
          Transpose[{collrates[[All, 1 ;; 2]], collratesTemp}];
        tmpMatrix = 
         SparseArray[Join[upwardCollisionRules, downwardCollisionRules]];
        diagonalRules = 
         Table[Rule[{i, i}, -Total[Normal[tmpMatrix[[All, i]]]]], {i, 1, 
           Length@levels}];
        tmpMatrix + 
         SparseArray[diagonalRules, {Length@levels, Length@levels}]
    ] 
    
(* same but accounting for multiple collision partners. Output is a list of coll matrices 
which have to be multiplied with their rerspective densities*)
SetUpCollisionMatrix2[temperature_?NumberQ] := 
 Module[{T = temperature, deltaEs, statweightratio, 
   interpolateCollRates, collratesTemp, upwardCollisionRules, 
   upCollratesTemp, downwardCollisionRules, tmpMatrix, diagonalRules, 
   mat},(*set up the collision matrix*)(*energy difference between \
upper and lower level in Kelvin*)
  interpolateCollRates[temp_, list_, colltemps_] := 
   ListInterpolation[list, {colltemps}][temp];
  mat = Table[
    deltaEs = (levelenergiesKelvin[[#[[1]]]] - 
         levelenergiesKelvin[[#[[2]]]]) & /@ 
      collrates[[i, All, 1 ;; 2]];
    (*ratio of the statistical weights to compute the upward rates \
from statistical equilibrium*)
    statweightratio = (levels[[#[[1]], 3]]/levels[[#[[2]], 3]] &) /@ 
      collrates[[i, All, 1 ;; 2]];
    (*interpolate the provided collision rates to the given \
temperature*)
    collratesTemp = 
     Map[interpolateCollRates[T, #[[3 ;; Length@colltemps[[i]] + 2]], 
        colltemps[[i]]] &, collrates[[i]]];
    (*calculate the upward rates*)
    upCollratesTemp = collratesTemp statweightratio*Exp[-deltaEs/T];
    (*convert rates into rules to populate the sparse array*)
    upwardCollisionRules = 
     Rule[#[[1]], #[[2]]] & /@ 
      Transpose[{collrates[[i, All, 1 ;; 2]], upCollratesTemp}];
    downwardCollisionRules = 
     Rule[Reverse@#[[1]], #[[2]]] & /@ 
      Transpose[{collrates[[i, All, 1 ;; 2]], collratesTemp}];
    tmpMatrix = 
     SparseArray[Join[upwardCollisionRules, downwardCollisionRules],{Length@levels, Length@levels}];
    diagonalRules = 
     Table[Rule[{i, i}, -Total[Normal[tmpMatrix[[All, i]]]]], {i, 1, 
       Length@levels}];
    tmpMatrix + 
     SparseArray[diagonalRules, {Length@levels, Length@levels}],
    {i, numberofCollisionPartners}]]


(* homogeneous, uniform sphere (Osterbrock, 1974, Astrophysics of Gaseous Nebulae) *)    
escapeProbability["Sphere"] = 
  Compile[{{tau, _Real}}, 
   If[ tau > 0.2 (* somewhat arbitrary limit, depending on how many terms 
       are included in the series below, this number can be changed *),
       1.5/tau (1. - 2./tau^2. + (2./tau + 2./tau^2.) Exp[-tau]),
       (* taylor series of the above form for small tau, because the above diverges for tau->0 *)
       1.` - 0.375` tau + 0.1` tau^2 - 0.020833333333333332` tau^3 + 
        0.0035714285714285718` tau^4 - 0.0005208333333333333` tau^5 + 
        0.00006613756613756614` tau^6 - 7.4404761904761905`*^-6 tau^7 + 
        7.515632515632515`*^-7 tau^8 - 6.889329805996473`*^-8 tau^9 + 
        5.781255781255781`*^-9 tau^10 - 4.473590783114593`*^-10 tau^11
   ]];
(* radially expanding sphere, aka LVG, aka Sobolov approximation 
(Elitzur, 1992, Astronomical Masers) *)
escapeProbability["LVG"] = 
  Compile[{{tau, _Real}}, If[ tau > 0,
                              (1. - Exp[-tau])/tau,
                              1. (* case tau==0, i.e. escape probability == 1 *)
                          ]];
(* homogeneous slab, i.e. semi infinite plane parallel geometry *) 
escapeProbability["PlaneParallel"] = 
  Compile[{{tau, _Real}}, 
   If[ tau > 0,
       (1. - Exp[-3. tau])/(3. tau),
       1.
   ]];
(* turbulent medium *)
(*escapeProbability["Turbulent"] = 
  Compile[{{tau, _Real}}, 1./(tau Sqrt[N[\[Pi]] Log[tau/2.]])];*)
escapeProbability[] = escapeProbability["Sphere"];
escapeProbability[_String] = escapeProbability["Sphere"];

(*escapeProbability[\[Tau]_, opts : OptionsPattern[]] :=
    Block[ {$MaxPrecision = 40, $MinPrecision = 40, geom, 
      geomList, \[Beta]List, tau},
        geomList = {"LVG", "PlaneParallel", "Turbulent", "Sphere"};
        \[Beta]List = {(1. - Exp[-tau])/tau,
          (1. - Exp[-3. tau])/(3. tau),
          1./(tau Sqrt[N[\[Pi]] Log[tau/2.]]),
          1.5/tau (1. - 2./tau^2 + (2./tau + 2./tau^2) Exp[-tau]),
          (*1.`-0.375` tau+0.1` tau^2-0.020833333333333332` \
        tau^3+0.0035714285714285718` tau^4,*)
          1.` - 0.375` tau + 0.1` tau^2 - 0.020833333333333332` tau^3 + 
           0.0035714285714285718` tau^4 - 0.0005208333333333333` tau^5 + 
           0.00006613756613756614` tau^6 - 7.4404761904761905`*^-6 tau^7 + 
           7.515632515632515`*^-7 tau^8 - 6.889329805996473`*^-8 tau^9 + 
           5.781255781255781`*^-9 tau^10 - 4.473590783114593`*^-10 tau^11,
          1.5/tau};
        geom = OptionValue[Geometry];
        If[ ! MemberQ[geomList, geom],
            Message[escapeProbability::optx, geom];
            geom = "LVG"
        ];
        Which[
         geom == geomList[[1]] && \[Tau] < 10^-8, 1.,
         geom == geomList[[1]], \[Beta]List[[1]] /. tau -> \[Tau],
         geom == geomList[[2]] && \[Tau] < 10^-8, 1.,
         geom == geomList[[2]], \[Beta]List[[2]] /. tau -> \[Tau],
         geom == geomList[[3]], \[Beta]List[[3]] /. tau -> \[Tau],
         geom == geomList[[4]] && \[Tau] < .001, \[Beta]List[[5]] /. 
          tau -> \[Tau],
         geom == geomList[[4]] && \[Tau] > 50., \[Beta]List[[6]] /. 
          tau -> \[Tau],
         geom == geomList[[4]], \[Beta]List[[4]] /. tau -> \[Tau]]
    ]*)
combineOrthoPara[orthoRate_, paraRate_, temperature_] :=
    Module[ {op1, op2, T1 = 18.73, T2 = 155.2},
    	(*  calculate a combined collision rates for cases where
    	ortho and para-H2 collision rates are given*) 
        Which[
         temperature < T1, op1 = 10^-3,
         T1 <= temperature <= T2, op1 = 9.0 Exp[-170.5/temperature],
         True, op1 = 3.0];
        op1 = op1/(1 + op1);
        op2 = 1.0 - op1;
        op2 paraRate + op1 orthoRate
    ]

lambdaInput[filename_String] :=
    Block[ {inp, pos, name, numLevels, levels, numTrans, numColPartners, 
      colType, numColTrans, numColTemps, colTemps, colTrans, clrat, cinp,
       trans, increment, orthoParaColTrans},
       (* import LAMBDA files that contain the atomic constants and the 
       collision rates *)
       (* see http://www.strw.leidenuniv.nl/~moldata/ *)
       
       (* only works for one or two collision partnerds! all files
       that give more partner will ignore all except first two *)
        inp = Import[filename, "Table"];

        pos = {5, 5 + 3 + inp[[6]], 
           5 + 3 + inp[[6]] + 3 + inp[[5 + 3 + inp[[6]] + 1]], 
           5 + 3 + inp[[6]] + 3 + inp[[5 + 3 + inp[[6]] + 1]] + 4, 
           5 + 3 + inp[[6]] + 3 + inp[[5 + 3 + inp[[6]] + 1]] + 6, 
           5 + 3 + inp[[6]] + 3 + inp[[5 + 3 + inp[[6]] + 1]] + 8} // 
    Flatten;
        (* molecule name line 2*)
        name = inp[[2]];
        (* numer of energy levels *)
        numLevels = inp[[pos[[1]] + 1]];
        (* energy levels from pos[[1]]+1;; pos[[1]]+1+numLevels-1*)
        levels = inp[[pos[[1]] + 3 ;; pos[[1]] + 3 + First@numLevels - 1]];
        (* numer of transitions *)
        numTrans = inp[[pos[[2]] + 1]];
        (* Transitions: !TRANS+UP+LOW+EINSTEINA (s^-1)+FREQ (GHz)+E_u (K)*)
        trans = inp[[pos[[2]] + 3 ;; pos[[2]] + 3 + First@numTrans - 1]];
        (* number of collision partners *)
        numColPartners = inp[[pos[[3]] + 1]];
        If[ numColPartners[[1]] == 1,
         (* type of collision *)
            colType = inp[[pos[[3]] + 3]];
            (*number of collisional transitions*)
            numColTrans = inp[[pos[[4]] + 1]];
            (* number of collision tempreature points*)
            numColTemps = inp[[pos[[5]] + 1]];
            (* collision temperature poinst*)
            colTemps = inp[[pos[[6]] + 1]];
            (* collisional transistions *)
            colTrans = 
             inp[[pos[[6]] + 3 ;; pos[[6]] + 3 + First@numColTrans - 1]];,
            colType = {};
            numColTemps = {};
            numColTrans = {};
            colTemps = {};
            colTrans = {{}};
            increment = 0;
            For[i = 1, i <= numColPartners[[1]], i++,
             {colType = {colType, inp[[pos[[3]] + 3 + increment]]};
              colType = Flatten /@ colType;
              (*number of collisional transitions*)
              numColTrans = {numColTrans, inp[[pos[[4]] + 1 + increment]]};
              numColTrans = Flatten[numColTrans];
              (* number of collision tempreature points*)
              numColTemps = {numColTemps, inp[[pos[[5]] + 1 + increment]]} // 
       Flatten;
              colTemps = {colTemps, inp[[pos[[6]] + 1 + increment]]};
              colTemps = Flatten /@ colTemps;
              (* collisional transistions *)
              colTrans = 
               Append[colTrans, 
                inp[[pos[[6]] + 3 + increment ;; 
                   pos[[6]] + 3 + Last[numColTrans] - 1 + increment]]];
              If[ i < First@numColPartners,
                  increment = First@inp[[pos[[4]] + 1 + increment]] + 9
              ];
              }];
            colTrans = Drop[colTrans, 1];
            If[ colTemps[[1]] == colTemps[[2]],
                orthoParaColTrans = 
                 Table[Flatten[
                   Prepend[Table[
                     combineOrthoPara[colTrans[[1, j]][[i + 3]], 
                      colTrans[[2, j]][[i + 3]], colTemps[[1, i]]], {i, 1, 
                      Length@colTemps[[1]]}], colTrans[[1, j]][[{1, 2, 3}]]]], {j,
                    1, Length@colTrans[[1]]}],
                Return[$Failed];
                Break
            ];
            colTemps = colTemps[[1]];
            colTrans = orthoParaColTrans;
        ];
        {levels, trans, colTemps, colTrans}
    ]
fileRules={"13co" -> "13co.dat", "13cs" -> "13cs@xpol.dat", 
 "a-ch3oh" -> "a-ch3oh.dat", "c17o" -> "c17o.dat", 
 "c18o" -> "c18o.dat", "c" -> "catom.dat", "c+" -> "c+.dat", 
 "co" -> "co.dat", "cs" -> "cs@xpol.dat", "e-ch3oh" -> "e-ch3oh.dat", 
 "h13cn" -> "h13cn@xpol.dat", "h13co+" -> "h13co+@xpol.dat", 
 "hc18o+" -> "hc18o+@xpol.dat", "hc3n" -> "hc3n.dat", 
 "hcl" -> "hcl.dat", "hcl@hfs" -> "hcl@hfs.dat", 
 "hcn" -> "hcn@xpol.dat", "hco+" -> "hco+@xpol.dat", 
 "hcs+" -> "hcs+@xpol.dat", "hnc" -> "hnc@xpol.dat", 
 "n2h+_hfs" -> "n2h+_hfs.dat", "no" -> "no.dat", "o2" -> "o2.dat", 
 "o" -> "oatom.dat", "o-c3h2" -> "o-c3h2.dat", 
 "ocs" -> "ocs@xpol.dat", "o-h2co" -> "o-h2co.dat", 
 "o-h2o" -> "o-h2o.dat", "oh2o-h2" -> "oh2o-h2.dat", 
 "o-h2o@lowT" -> "o-h2o@lowT.dat", "o-h3o+" -> "o-h3o+.dat", 
 "oh" -> "oh@hfs.dat", "o-nh3" -> "o-nh3.dat", 
 "o-sic2" -> "o-sic2.dat", "p-c3h2" -> "p-c3h2.dat", 
 "p-h2co" -> "p-h2co.dat", "p-h2o" -> "p-h2o.dat", 
 "ph2o-h2" -> "ph2o-h2.dat", "p-h2o@lowT" -> "p-h2o@lowT.dat", 
 "p-h3o+" -> "p-h3o+.dat", "p-nh3" -> "p-nh3.dat", "sio" -> "sio.dat",
  "sis" -> "sis@xpol.dat", "so2" -> "so2@xpol.dat", "so" -> "so.dat", 
 "13CO" -> "13co.dat", "13CS" -> "13cs@xpol.dat", 
 "a-CH3OH" -> "a-ch3oh.dat", "C17O" -> "c17o.dat", 
 "C18O" -> "c18o.dat", "C" -> "catom.dat", "C+" -> "c+.dat", 
 "CO" -> "co.dat", "CS" -> "cs@xpol.dat", "e-CH3OH" -> "e-ch3oh.dat", 
 "H13CN" -> "h13cn@xpol.dat", "H13CO+" -> "h13co+@xpol.dat", 
 "HC18O+" -> "hc18o+@xpol.dat", "HC3N" -> "hc3n.dat", 
 "HCl" -> "hcl.dat", "HCl@hfs" -> "hcl@hfs.dat", 
 "HCN" -> "hcn@xpol.dat", "HCO+" -> "hco+@xpol.dat", 
 "HCS+" -> "hcs+@xpol.dat", "HNC" -> "hnc@xpol.dat", 
 "N2H+_hfs" -> "n2h+_hfs.dat", "NO" -> "no.dat", "O2" -> "o2.dat", 
 "O" -> "oatom.dat", "o-C3H2" -> "o-c3h2.dat", 
 "OCS" -> "ocs@xpol.dat", "o-H2CO" -> "o-h2co.dat", 
 "o-H2O" -> "o-h2o.dat", "oH2O-H2" -> "oh2o-h2.dat", 
 "o-H2O@lowT" -> "o-h2o@lowT.dat", "o-H3O+" -> "o-h3o+.dat", 
 "OH" -> "oh@hfs.dat", "o-NH3" -> "o-nh3.dat", 
 "o-SiC2" -> "o-sic2.dat", "p-C3H2" -> "p-c3h2.dat", 
 "p-H2CO" -> "p-h2co.dat", "p-H2O" -> "p-h2o.dat", 
 "pH2O-H2" -> "ph2o-h2.dat", "p-H2O@lowT" -> "p-h2o@lowT.dat", 
 "p-H3O+" -> "p-h3o+.dat", "p-NH3" -> "p-nh3.dat", "SiO" -> "sio.dat",
  "SiS" -> "sis@xpol.dat", "SO2" -> "so2@xpol.dat", "SO" -> "so.dat"};
    
lambdaInput2[filename_String, opts : OptionsPattern[]] := 
 Block[{inp, pos, name, numLevels, levels, numTrans, numColPartners, 
   colType, numColTrans, numColTemps, colTemps, colTrans, clrat, cinp,
    trans, increment, orthoParaColTrans, combineorthopara, opos, ppos,file},
  combineorthopara = 
   CombineOrthoPara /. {opts} /. Options[lambdaInput2];
  (*import LAMBDA files that contain the atomic constants and the \
collision rates*)(*see http://www.strw.leidenuniv.nl/~
  moldata/*)(*only works for one or two collision partnerds! all \
files that give more partner will ignore all except first two*)
  file=filename/.fileRules;
 (* TODO: check if file exists *)
  inp = Import[file, "Table"];
  pos = {5, 5 + 3 + inp[[6]], 
     5 + 3 + inp[[6]] + 3 + inp[[5 + 3 + inp[[6]] + 1]], 
     5 + 3 + inp[[6]] + 3 + inp[[5 + 3 + inp[[6]] + 1]] + 4, 
     5 + 3 + inp[[6]] + 3 + inp[[5 + 3 + inp[[6]] + 1]] + 6, 
     5 + 3 + inp[[6]] + 3 + inp[[5 + 3 + inp[[6]] + 1]] + 8} // 
    Flatten;
  (*molecule name line 2*)name = inp[[2]];
  (*numer of energy levels*)numLevels = inp[[pos[[1]] + 1]];
  (*energy levels from pos[[1]]+1;;pos[[1]]+1+numLevels-1*)
  levels = inp[[pos[[1]] + 3 ;; pos[[1]] + 3 + First@numLevels - 1]];
  (*numer of transitions*)numTrans = inp[[pos[[2]] + 1]];
  (*Transitions:!TRANS+UP+LOW+EINSTEINA (s^-1)+FREQ (GHz)+E_u (K)*)
  trans = inp[[pos[[2]] + 3 ;; pos[[2]] + 3 + First@numTrans - 1]];
  (*number of collision partners*)numColPartners = inp[[pos[[3]] + 1]];
  Which[numColPartners[[1]] == 1,
   (*type of collision*)colType = inp[[pos[[3]] + 3]];
   (*number of collisional transitions*)
   numColTrans = inp[[pos[[4]] + 1]];
   (*number of collision tempreature points*)
   numColTemps = inp[[pos[[5]] + 1]];
   (*collision temperature poinst*)colTemps = inp[[pos[[6]] + 1]];
   (*collisional transistions*)
   colTrans = 
    inp[[pos[[6]] + 3 ;; pos[[6]] + 3 + First@numColTrans - 1]];,
   numColPartners[[1]] == 2 && combineorthopara,
   colType = {};
   numColTemps = {};
   numColTrans = {};
   colTemps = {};
   colTrans = {{}};
   increment = 0;
   For[i = 1, i <= numColPartners[[1]], 
    i++, {colType = {colType, inp[[pos[[3]] + 3 + increment]]};
     colType = Flatten /@ colType;
     (*number of collisional transitions*)
     numColTrans = {numColTrans, inp[[pos[[4]] + 1 + increment]]};
     numColTrans = Flatten[numColTrans];
     (*number of collision tempreature points*)
     numColTemps = {numColTemps, inp[[pos[[5]] + 1 + increment]]} // 
       Flatten;
     colTemps = {colTemps, inp[[pos[[6]] + 1 + increment]]};
     colTemps = Flatten /@ colTemps;
     (*collisional transistions*)
     colTrans = 
      Append[colTrans, 
       inp[[pos[[6]] + 3 + increment ;; 
          pos[[6]] + 3 + Last[numColTrans] - 1 + increment]]];
     If[i < First@numColPartners, 
      increment = First@inp[[pos[[4]] + 1 + increment]] + 9];}];
   colTrans = Drop[colTrans, 1];
   If[colTemps[[1]] == colTemps[[2]], 
    orthoParaColTrans = 
     Table[Flatten[
       Prepend[Table[
         combineOrthoPara[colTrans[[1, j]][[i + 3]], 
          colTrans[[2, j]][[i + 3]], colTemps[[1, i]]], {i, 1, 
          Length@colTemps[[1]]}], colTrans[[1, j]][[{1, 2, 3}]]]], {j,
        1, Length@colTrans[[1]]}], Return[$Failed];
    Break];
   colTemps = colTemps[[1]];
   colTrans = orthoParaColTrans;
   numColPartners[[1]] = 1,
   True,
   (*Print["More than 2 collision partners!"];*)
   colType = {};
   numColTemps = {};
   numColTrans = {};
   colTemps = {};
   colTrans = {};
   increment = 0;
   For[i = 1, i <= numColPartners[[1]], 
    i++, {colType = Append[colType, inp[[pos[[3]] + 3 + increment]]];
     (*colType=Flatten/@colType;*)
     (*number of collisional transitions*)
     numColTrans = {numColTrans, inp[[pos[[4]] + 1 + increment]]};
     numColTrans = Flatten[numColTrans];
     (*number of collision tempreature points*)
     numColTemps = {numColTemps, inp[[pos[[5]] + 1 + increment]]} // 
       Flatten;
     (*colTemps={colTemps,inp[[pos[[6]]+1+increment]]};*)
     colTemps = Append[colTemps, inp[[pos[[6]] + 1 + increment]]];
     (*colTemps=Flatten/@colTemps;*)
     (*collisional transistions*)
     colTrans = 
      Append[colTrans, 
       inp[[pos[[6]] + 3 + increment ;; 
          pos[[6]] + 3 + Last[numColTrans] - 1 + increment]]];
     If[i < First@numColPartners, 
      increment = 
       increment + First@inp[[pos[[4]] + 1 + increment]] + 9]; 
     (*Print[increment];*)}];
   ppos = 
    Position[
      Pick[Select[#, StringQ], 
         StringFreeQ[Select[#, StringQ], 
          "+" | "-" | "" ~~ ___ ~~ "p" ~~ ___ ~~ "h2", 
          IgnoreCase -> True], False] & /@ colType, 
      x_List /; x =!= {}] // Flatten;
   opos = 
    Position[
      Pick[Select[#, StringQ], 
         StringFreeQ[Select[#, StringQ], 
          "+" | "-" | "" ~~ ___ ~~ "o" ~~ ___ ~~ "h2", 
          IgnoreCase -> True], False] & /@ colType, 
      x_List /; x =!= {}] // Flatten;
  (* Print[{opos, ppos}];
   Print[colTemps[[opos[[1]]]]];*)
   If[ppos =!= {} && opos =!= {} && ppos =!= opos &&combineorthopara,
    If[colTemps[[opos[[1]]]] == colTemps[[ppos[[1]]]],
     (*Print["HERE"];*)
     orthoParaColTrans = 
      Table[Flatten[
        Prepend[Table[
          combineOrthoPara[colTrans[[opos[[1]], j]][[i + 3]], 
           colTrans[[ppos[[1]], j]][[i + 3]], 
           colTemps[[opos[[1]], i]]], {i, 1, 
           Length@colTemps[[opos[[1]]]]}], 
         colTrans[[opos[[1]], j]][[{1, 2, 3}]]]], {j, 1, 
        Length@colTrans[[opos[[1]]]]}], Return[$Failed];
     Break];
    (*Print[orthoParaColTrans];*)
    (* remove one set of collision temperatures*)
    colTemps = Delete[colTemps, opos[[1]]];
    (*Print[colTrans];*)
    (* replace both o/p H2 rates with the combined rate *)
    colTrans = 
     Delete[ReplacePart[colTrans, ppos -> orthoParaColTrans], 
      opos[[1]]];
    numColPartners[[1]] = numColPartners[[1]] - 1;
    ];
   
   (*Print[colType];*)
   ];
  {levels, trans, colTemps, colTrans, numColPartners[[1]],colType}]
  

ShowCollisionPartners[filename_String] := 
 Block[{levels, transitions, colltemps, collrates, 
   numberofCollisionPartners, colType, ppos, opos},
  {levels, transitions, colltemps, collrates, 
    numberofCollisionPartners, colType} = lambdaInput2[filename/.fileRules];
  If[VectorQ[colType], colType = {colType}];
  ppos = Position[
     Pick[Select[#, StringQ], 
        StringFreeQ[Select[#, StringQ], 
         "+" | "-" | "" ~~ ___ ~~ "p" ~~ ___ ~~ "h2", 
         IgnoreCase -> True], False] & /@ colType, 
     x_List /; x =!= {}] // Flatten;
  opos = Position[
     Pick[Select[#, StringQ], 
        StringFreeQ[Select[#, StringQ], 
         "+" | "-" | "" ~~ ___ ~~ "o" ~~ ___ ~~ "h2", 
         IgnoreCase -> True], False] & /@ colType, 
     x_List /; x =!= {}] // Flatten;
  If[opos =!= {} && ppos =!= {}, 
   Message[ShowCollisionPartners::combop]];
  colType]
  
ShowCollisionTemperatures[filename_String, opts : OptionsPattern[]] := 
 Block[{levels, transitions, colltemps, collrates, 
   numberofCollisionPartners, colType, ppos, opos},
  {levels, transitions, colltemps, collrates, 
    numberofCollisionPartners, colType} = lambdaInput2[filename/.fileRules,opts];
  colltemps]


ShowCollisionRates[filename_String, opts : OptionsPattern[]] := 
 Block[{levels, transitions, colltemps, collrates, 
   numberofCollisionPartners, colType, ppos, opos},
  {levels, transitions, colltemps, collrates, 
    numberofCollisionPartners, colType} = lambdaInput2[filename/.fileRules,opts];
  collrates]
ShowEnergyLevels[filename_String] := 
 Block[{levels, transitions, colltemps, collrates, 
   numberofCollisionPartners, colType, ppos, opos},
  {levels, transitions, colltemps, collrates, 
    numberofCollisionPartners, colType} = lambdaInput2[filename/.fileRules];
  levels]
ShowTransitions[filename_String] := 
 Block[{levels, transitions, colltemps, collrates, 
   numberofCollisionPartners, colType, ppos, opos},
  {levels, transitions, colltemps, collrates, 
    numberofCollisionPartners, colType} = lambdaInput2[filename/.fileRules];
  transitions]  

(* EscapeProbabilityRun has to be a Block to allow dynamic scoping,
 i.e., all functions called from a Block[] can access its local data*)
 
EscapeProbabilityRun[lambdaFilename_String, columnDensity_, 
  deltavkmsec_, den_, temperature_, opts : OptionsPattern[]] :=
    Block[ {h, cl, kB, einsteinARules, einsteinARulesDiagonal, 
      einsteinAMatrix, BRules3, BRules4, BRules2, BRules1, bmatrix, 
      coeffMatrix, betavektor, beta, background3k, levels, transitions, 
      colltemps, collrates, bvalues21, bvalues12, 
      normierteBesetzungszahlen, background, bnuTex, exr, 
      radiativeRateMatrix, fk, thc, fgaus, xnu, minpop, \[CapitalDelta]v,
       startbesetzung, collisionMatrix, in, gstat, hnu, 
      levelenergiesKelvin, wavenumber, column, temp, density, maxiter, 
      finalPopulation, finalTau, finalLineBrightness, 
      finalAntennaTemperature, finalRadiationTemperature, 
      finalSourceFunction, finalTotalBackground, 
      finalTotalLineBrightness, finalKkms, finalergs, iterateOpts, out, 
      outTable, tableHeader, tableHeader2, freqRange, showPop, grid, 
      showTA, fullOut, tableOut, showTEX, showTAU, tbg,DEBUG,iterCount,xnut,
      sametest,tabLowLow,inipoptemp,avalues,levelEnergy,numberofCollisionPartners,
      colType,rangePos,showBETA,geom,comborpa},
        DEBUG = False;
        h = 6.62606896 10^-27;
        cl = 2.99792456 10^10;
        kB = 1.3806504 10^-16;
        (*iterateOpts=FilterRules[{opts},Options[iterate]];*)
        geom = 
         Geometry /. {opts} /. Options[EscapeProbabilityRun];
        comborpa = 
         CombineOrthoPara /. {opts} /. Options[EscapeProbabilityRun];
        minpop = 
         MinimumLevelPopulation /. {opts} /. Options[EscapeProbabilityRun];
        maxiter = 
         MaximumIterations /. {opts} /. Options[EscapeProbabilityRun];
        freqRange = 
         FrequencyRange /. {opts} /. Options[EscapeProbabilityRun];
        showPop = 
         ShowLevelPopulation /. {opts} /. Options[EscapeProbabilityRun];
        showTA = 
         ShowAntennaTemperature /. {opts} /. Options[EscapeProbabilityRun];
        fullOut = FullOutput /. {opts} /. Options[EscapeProbabilityRun];
        tableOut = TableOutput /. {opts} /. Options[EscapeProbabilityRun];
        showTEX = 
         ShowExcitationTemperature /. {opts} /. Options[EscapeProbabilityRun];
        showTAU = 
         ShowOpticalDepth /. {opts} /. Options[EscapeProbabilityRun];
        showBETA = 
         ShowEscapeProbability /. {opts} /. Options[EscapeProbabilityRun];
        tbg = BackgroundTemperature /. {opts} /. 
          Options[EscapeProbabilityRun];
        sametest = SameTest/. {opts} /. Options[EscapeProbabilityRun];
        inipoptemp = InitialPopulationTemperature/. {opts} /. Options[EscapeProbabilityRun];
        If[ fullOut && Or[showPop, showTA, showTEX, showTAU, showBETA],
            fullOut = False
        ];
        If[ showPop && showTA,
            Message[EscapeProbabilityRun::conflopts, ShowLevelPopulation, 
             ShowAntennaTemperature, ShowAntennaTemperature];
            showTA = False
        ];
        If[ showPop && showTEX,
            Message[EscapeProbabilityRun::conflopts, ShowLevelPopulation, 
             ShowExcitationTemperature, ShowExcitationTemperature];
            showTEX = False
        ];
        If[ showPop && showTAU,
            Message[EscapeProbabilityRun::conflopts, ShowLevelPopulation, 
             ShowOpticalDepth, ShowOpticalDepth];
            showTAU = False
        ];
         If[ showPop && showBETA,
            Message[EscapeProbabilityRun::conflopts, ShowLevelPopulation, 
             ShowEscapeProbability, ShowEscapeProbability];
            showBETA = False
        ];
        (* make all input parameters available in all functions being called within this Block *)
        column = columnDensity;
        density = den;
        temp = temperature;
        \[CapitalDelta]v = deltavkmsec 10^5;
        fk = (h cl)/kB;
        (*minpop=1. 10^-30;*)
        thc = 2. h cl;
        fgaus = 1.0645 8.0 N[Pi];
        (* Import LAMBDA file *)
        {levels, transitions, colltemps, collrates,numberofCollisionPartners, colType} = 
         lambdaInput2[lambdaFilename,CombineOrthoPara->comborpa];
        (* find positions that are in the given frequency range *)
         rangePos=Position[transitions[[All, 5]], 
 a_ /; IntervalMemberQ[Interval@freqRange, a]];
        (* remove running numbers from list of transitions*)
        If[numberofCollisionPartners==1, 
        	collrates = {Rest /@ collrates}; colltemps = {colltemps};,
       	    collrates = Map[Rest /@#&,collrates]
       	  ];     	         	
        wavenumber[{upper_, lower_}] :=
            levels[[upper, 2]] - levels[[lower, 2]];
        (* assign to variables *)
        (* xnu = transition frequency in wavenumbers [cm^-1]*)
         xnu = Table[
          wavenumber[transitions[[i]][[{2, 3}]]], {i, 1, 
           Length@transitions}];
        (* frequency ^3 [cm^-3]*)        
         xnut = xnu^3;
        (* statistical weights of each level *)
        gstat = levels[[All, 3]];
        (* level enegy in [cm^-1] *)
        levelEnergy=levels[[All,2]];
        (* level energy in [K] *)
        levelenergiesKelvin = levels[[All, 2]] (h cl)/kB;
        (* einstein A values [s^-1]*)
        avalues=transitions[[All,4]];
        (* factor for radiative rates : A_21*g_2/g_1 *)
        tabLowLow = Table[transitions[[i, 4]] gstat[[transitions[[i, 2]]]]/
                   gstat[[transitions[[i, 3]]]], {i, Length@transitions}];
        
        in = {levels, transitions, colltemps, collrates};
        (* calculate total background intensity *)
        background = totalBackground[tbg];
        (* calculate intial level population *)
        (* changed initial population to T->background*)
        startbesetzung = populateLTE[inipoptemp];
        (* populate collision matrix - this has to be done only once*)
        (* Test whether number of collision partners and number of
        densities is the same*)
        If[(*numberofCollisionPartners==1&&*)density[[0]]=!=List,density={density}];
        If[numberofCollisionPartners=!=Length[density],
        	Message[EscapeProbabilityRun::denxmatch,numberofCollisionPartners];
        	density=PadRight[density,numberofCollisionPartners,density]];
        (*density*SetUpCollisionMatrix[temp];*)
        (* multiply each coll matrix with its density and add all up*)
        collisionMatrix = Apply[Plus,density*SetUpCollisionMatrix2[temp]];
        If[ DEBUG,
            Print["Initial Population:\n", startbesetzung]
        ];
        (* initialize iteration counter *)
        iterCount = 0;
        (* This is the actual iteration *)
        Quiet[
         finalPopulation = 
           Normalize@FixedPoint[
            iterate[N[#, 30](*,density,column,\[CapitalDelta]v,
              collisionMatrix,in*), opts] &, startbesetzung, maxiter,
            SameTest->sametest], {Power::"infy"}];
        (* calculate all output quantities from the final population *)
        Quiet[out = calculateOutput[finalPopulation, in, opts, 
        	Geometry->geom], {Power::"infy"}];
        Print["Number of Iterations:",iterCount];
        (* prepare table of final output *)
        outTable = Table[
          {levels[[transitions[[i, 2]], 4]] (*quantum numbers of upper level *),
           levels[[transitions[[i, 3]], 4]] (*quantum numbers of lower level *),
           transitions[[i, 6]] (* Energy of upper level in [K] *),
           transitions[[i, 5]] (* transition frequency in [GHz]*),
           cl/transitions[[i, 5]]/10^5 (* transition wavelength in [micro-m]*),
           out[[1, i]] (* ectitation temperature in [K]*),
           out[[2, i]] (* optical depth t *),
           out[[3, i]] (* antenna temperature against background [K] *),
           out[[4, i]] (* velocity integrated flux in [K km/s]*),
           out[[5, i]] (* velocity integrated flux in [erg/s/cm^2/sr]*)
           }, {i, 1, Length@transitions}];
        (* select all lines that are part of the wanted frequency range *)
        If[ ListQ[freqRange] && Length[freqRange] == 2,
            outTable = 
             Select[outTable, freqRange[[1]] <= #[[4]] <= freqRange[[2]] &]
        ];
        (* set up table headings *)
        tableHeader = {"Line", SpanFromLeft, " E_UP", "Frequency", 
          " Wavelength", "T_EX", "Tau", "T_R", "Flux", "Flux"};
        tableHeader2 = {"", "", "[K]", "[GHz]", "[\[Micro]m]", "[K]" , "", 
          "[K]", "[K km/s]", "[erg/cm2/s]"};
        (* prepare table form *)
        grid = Grid[Prepend[Prepend[outTable, tableHeader2], tableHeader], 
          Background -> {None, {Lighter[Yellow, .9], 
             Lighter[Yellow, .9], {White, 
              Lighter[Blend[{Blue, Green}], .8]}}}, 
          Dividers -> {{Darker[Gray, .6], {Lighter[Gray, .5]}, 
             Darker[Gray, .6]}, {Darker[Gray, .6], Darker[Gray, .6], 
             Darker[Gray, .6], {False}, Darker[Gray, .6]}}, 
          Alignment -> {{Center, Center, Left, Left, 
             Right, {Left}}},(*ItemSize->{{10,3,5,5}},*)
          Frame -> Darker[Gray, .6],(* ItemStyle -> 14, *)
          Spacings -> {Automatic, .8}];
        
        Which[
         showPop,(finalPopulation/(Total@finalPopulation))[[Extract[transitions, rangePos][[All, 2 ;; 3]] // Flatten // Union]] (* only give out final level population *),
         showTA && showTEX && showTAU, Extract[out[[{1, 2, 3}]],rangePos] (* print Tex, tau and TA*),
         showTA && showTEX, Extract[out[[{1, 3}]],rangePos] (* print Tex and TA*),
         showTA && showTAU, Extract[out[[{2, 3}]],rangePos] (*print tau and TA*),
         showTEX && showTAU, Extract[out[[{1, 2}]],rangePos] (*print TEX and tau *),
         showTA, Extract[out[[3]],rangePos] (*print TA*),
         showTAU, Extract[out[[2]],rangePos] (*print TAU*),
         showTEX, Extract[out[[1]],rangePos] (* print Tex *),
         showBETA, Extract[out[[6]],rangePos] (* print beta *),
         tableOut, grid (* print formatted table output *),
         True, outTable (* print unformatted output *)]
    ]
  
 calculateOutput[besetzungsvektor_, onionFileInput_, 
  opts : OptionsPattern[]] :=
     Block[ {ftau, toti, background, tbl, tback, wh, levels, transitions, 
       colltemps, collrates, h, cl, kB, fk, thc, fgaus, wavenumber, xnu, 
       hnu, ta, tr, finalTau, finalSourceFunction, bnu, beta,geom},
         {levels, transitions, colltemps, collrates} = onionFileInput;
         h = 6.62606896 10^-27;
         cl = 2.99792456 10^10;
         kB = 1.3806504 10^-16;
         fk = (h cl)/kB;
         (*minpop=1. 10^-30;*)
         thc = 2. h cl;
         fgaus = 1.0645 8.0 N[Pi];
         wavenumber[{upper_, lower_}] :=
             levels[[upper, 2]] - levels[[lower, 2]];
         xnu = Table[
           wavenumber[transitions[[i]][[{2, 3}]]], {i, 1, 
            Length@transitions}];
         finalTau = (column/\[CapitalDelta]v tauListe[besetzungsvektor]);
         ftau = If[ # <= 300.,
                    Exp[-#],
                    0.
                ] & /@ finalTau;
         background = totalBackground[tbg];
         (*Calculate source function*)
         finalSourceFunction = calcSource[besetzungsvektor];
         toti = background*ftau + finalSourceFunction (1 - ftau);
         wh = (thc xnu^3)/toti + 1.;
         geom = Geometry/.{opts};
         beta = escapeProbability[geom][#]&/@finalTau;
        (* beta = escapeProbability[#, 
             FilterRules[{opts}, Options[escapeProbability]]] & /@ finalTau;*)
         (*Calculate line brightness in excess of background*)
         tbl = Table[
           If[ toti[[i]] == 0,
               0.,
               If[ wh[[i]] <= 0.,
                   toti[[i]]/(thc xnu[[i]] xnu[[i]]/fk),
                   fk xnu[[i]]/Log[10, wh[[i]]]
               ]
           ], {i, 1, Length@toti}];
         tback = 
          Table[If[ background[[i]] == 0,
                    0.,
                    fk xnu[[i]]/Log[10, (thc xnu[[i]]^3)/background[[i]] + 1.]
                ], {i, 
            1, Length@background}];
         tbl = tbl - tback;
         hnu = fk xnu;
         (*Calculate antenna temperature*)
         ta = Table[
           If[ Abs[tback[[i]]/hnu[[i]]] <= 0.02,
               toti[[i]],
               toti[[i]] - background[[i]]
           ], {i, 1, Length@background}];
         ta = ta/(thc xnu^2/fk);
         bnu = (background beta + (1 - beta) finalSourceFunction);
         (*Calculate radiation temperature*)
         tr = Table[
           If[ bnu[[i]] == 0,
               background[[i]],
               If[ wh[[i]] <= 0,
                   bnu[[i]]/(thc xnu[[i]]^2/fk),
                   fk xnu[[i]]/Log[10, wh[[i]]]
               ]
           ], {i, 1, Length@background}];
         {calcTex[besetzungsvektor], finalTau, ta, 
          1.0645 \[CapitalDelta]v ta/10^5, 
          fgaus kB \[CapitalDelta]v ta xnu^3,
          beta}
     ]
     

End[]

EndPackage[]

