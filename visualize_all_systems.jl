using ChaoticDynamicalSystemLibrary, OrdinaryDiffEq, Plots

ALL_SYSTEMS = (ChaoticDynamicalSystemLibrary.Aizawa,
    ChaoticDynamicalSystemLibrary.AnishchenkoAstakhov,
    ChaoticDynamicalSystemLibrary.Arneodo,
    ChaoticDynamicalSystemLibrary.ArnoldBeltramiChildress,
    ChaoticDynamicalSystemLibrary.ArnoldWeb,
    ChaoticDynamicalSystemLibrary.AtmosphericRegime,
    ChaoticDynamicalSystemLibrary.BeerRNN,
    ChaoticDynamicalSystemLibrary.BelousovZhabotinsky,
    ChaoticDynamicalSystemLibrary.BickleyJet,
    ChaoticDynamicalSystemLibrary.Blasius,
    ChaoticDynamicalSystemLibrary.BlinkingRotlet,
    ChaoticDynamicalSystemLibrary.BlinkingVortex,
    ChaoticDynamicalSystemLibrary.Bouali,
    ChaoticDynamicalSystemLibrary.Bouali2,
    ChaoticDynamicalSystemLibrary.BurkeShaw,
    ChaoticDynamicalSystemLibrary.CaTwoPlus,
    ChaoticDynamicalSystemLibrary.CaTwoPlusQuasiperiodic,
    ChaoticDynamicalSystemLibrary.CellCycle,
    ChaoticDynamicalSystemLibrary.CellularNeuralNetwork,
    ChaoticDynamicalSystemLibrary.Chen,
    ChaoticDynamicalSystemLibrary.ChenLee,
    ChaoticDynamicalSystemLibrary.Chua,
    ChaoticDynamicalSystemLibrary.CircadianRhythm,
    ChaoticDynamicalSystemLibrary.CoevolvingPredatorPrey,
    ChaoticDynamicalSystemLibrary.Colpitts,
    ChaoticDynamicalSystemLibrary.Coullet,
    ChaoticDynamicalSystemLibrary.Dadras,
    ChaoticDynamicalSystemLibrary.DequanLi,
    ChaoticDynamicalSystemLibrary.DoubleGyre,
    ChaoticDynamicalSystemLibrary.DoublePendulum,
    ChaoticDynamicalSystemLibrary.Duffing,
    ChaoticDynamicalSystemLibrary.ExcitableCell,
    ChaoticDynamicalSystemLibrary.Finance,
    ChaoticDynamicalSystemLibrary.FluidTrampoline,
    ChaoticDynamicalSystemLibrary.ForcedBrusselator,
    ChaoticDynamicalSystemLibrary.ForcedFitzHughNagumo,
    ChaoticDynamicalSystemLibrary.ForcedVanDerPol,
    ChaoticDynamicalSystemLibrary.GenesioTesi,
    ChaoticDynamicalSystemLibrary.GuckenheimerHolmes,
    ChaoticDynamicalSystemLibrary.Hadley,
    ChaoticDynamicalSystemLibrary.Halvorsen,
    ChaoticDynamicalSystemLibrary.HastingsPowell,
    ChaoticDynamicalSystemLibrary.HenonHeiles,
    ChaoticDynamicalSystemLibrary.HindmarshRose,
    ChaoticDynamicalSystemLibrary.Hopfield,
    ChaoticDynamicalSystemLibrary.HyperBao,
    ChaoticDynamicalSystemLibrary.HyperCai,
    ChaoticDynamicalSystemLibrary.HyperJha,
    ChaoticDynamicalSystemLibrary.HyperLorenz,
    ChaoticDynamicalSystemLibrary.HyperLu,
    ChaoticDynamicalSystemLibrary.HyperPang,
    ChaoticDynamicalSystemLibrary.HyperQi,
    ChaoticDynamicalSystemLibrary.HyperRossler,
    ChaoticDynamicalSystemLibrary.HyperWang,
    ChaoticDynamicalSystemLibrary.HyperXu,
    ChaoticDynamicalSystemLibrary.HyperYan,
    ChaoticDynamicalSystemLibrary.HyperYangChen,
    ChaoticDynamicalSystemLibrary.IsothermalChemical,
    ChaoticDynamicalSystemLibrary.ItikBanksTumor,
    ChaoticDynamicalSystemLibrary.JerkCircuit,
    ChaoticDynamicalSystemLibrary.KawczynskiStrizhak,
    ChaoticDynamicalSystemLibrary.Laser,
    ChaoticDynamicalSystemLibrary.LidDrivenCavityFlow,
    ChaoticDynamicalSystemLibrary.Lorenz,
    ChaoticDynamicalSystemLibrary.Lorenz84,
    ChaoticDynamicalSystemLibrary.Lorenz96,
    ChaoticDynamicalSystemLibrary.LorenzBounded,
    ChaoticDynamicalSystemLibrary.LorenzCoupled,
    ChaoticDynamicalSystemLibrary.LorenzStenflo,
    ChaoticDynamicalSystemLibrary.LiuChen,
    ChaoticDynamicalSystemLibrary.LuChen,
    ChaoticDynamicalSystemLibrary.LuChenCheng,
    ChaoticDynamicalSystemLibrary.MacArthur,
    ChaoticDynamicalSystemLibrary.MooreSpiegel,
    ChaoticDynamicalSystemLibrary.MultiChua,
    ChaoticDynamicalSystemLibrary.NewtonLiepnik,
    ChaoticDynamicalSystemLibrary.NoseHoover,
    ChaoticDynamicalSystemLibrary.NuclearQuadrupole,
    ChaoticDynamicalSystemLibrary.OscillatingFlow,
    ChaoticDynamicalSystemLibrary.PanXuZhou,
    ChaoticDynamicalSystemLibrary.PehlivanWei,
    ChaoticDynamicalSystemLibrary.Qi,
    ChaoticDynamicalSystemLibrary.QiChen,
    ChaoticDynamicalSystemLibrary.RabinovichFabrikant,
    ChaoticDynamicalSystemLibrary.RayleighBenard,
    ChaoticDynamicalSystemLibrary.RikitakeDynamo,
    ChaoticDynamicalSystemLibrary.Rossler,
    ChaoticDynamicalSystemLibrary.Rucklidge,
    ChaoticDynamicalSystemLibrary.Sakarya,
    ChaoticDynamicalSystemLibrary.SaltonSea,
    ChaoticDynamicalSystemLibrary.SanUmSrisuchinwong,
    ChaoticDynamicalSystemLibrary.ShimizuMorioka,
    ChaoticDynamicalSystemLibrary.SprottA,
    ChaoticDynamicalSystemLibrary.SprottB,
    ChaoticDynamicalSystemLibrary.SprottC,
    ChaoticDynamicalSystemLibrary.SprottD,
    ChaoticDynamicalSystemLibrary.SprottE,
    ChaoticDynamicalSystemLibrary.SprottF,
    ChaoticDynamicalSystemLibrary.SprottG,
    ChaoticDynamicalSystemLibrary.SprottH,
    ChaoticDynamicalSystemLibrary.SprottI,
    ChaoticDynamicalSystemLibrary.SprottJ,
    ChaoticDynamicalSystemLibrary.SprottJerk,
    ChaoticDynamicalSystemLibrary.SprottK,
    ChaoticDynamicalSystemLibrary.SprottL,
    ChaoticDynamicalSystemLibrary.SprottM,
    ChaoticDynamicalSystemLibrary.SprottMore,
    ChaoticDynamicalSystemLibrary.SprottN,
    ChaoticDynamicalSystemLibrary.SprottO,
    ChaoticDynamicalSystemLibrary.SprottP,
    ChaoticDynamicalSystemLibrary.SprottQ,
    ChaoticDynamicalSystemLibrary.SprottR,
    ChaoticDynamicalSystemLibrary.SprottS,
    ChaoticDynamicalSystemLibrary.SprottTorus,
    ChaoticDynamicalSystemLibrary.StickSlipOscillator,
    ChaoticDynamicalSystemLibrary.SwingingAtwood,
    ChaoticDynamicalSystemLibrary.Thomas,
    ChaoticDynamicalSystemLibrary.Torus,
    ChaoticDynamicalSystemLibrary.Tsucs2,
    ChaoticDynamicalSystemLibrary.TurchinHanski,
    ChaoticDynamicalSystemLibrary.VallisElNino,
    ChaoticDynamicalSystemLibrary.WangSun,
    ChaoticDynamicalSystemLibrary.WindmiReduced,
    ChaoticDynamicalSystemLibrary.YuWang,
    ChaoticDynamicalSystemLibrary.YuWang2,
    ChaoticDynamicalSystemLibrary.ZhouChen)
wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)
for system in ALL_SYSTEMS
    @info "$system"
    prob = system()
    # sol = solve(prob, Tsit5(), tspan=(0, 1000))
    sol = solve(prob, Tsit5(), abstol = 1e-7, reltol = 1e-5)
    # p = plot(sol)
    p = plot(sol, layout = (length(sol.u[1]), 1), xlabel = "", xticks = false,
        color = collect(1:length(sol.u[1]))', legend = false,
        title = [i == 1 ? "$system" : "" for i in 1:length(sol.u[1])] |> permutedims)
    display(p)
    wait_for_key("Press any key to continue")
end
