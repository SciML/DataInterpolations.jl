import StableRNGs: StableRNG
using RegularizationTools

# create scattered data
npts = 50
xmin = 0.0
xspan = 3 / 2 * π
x = collect(range(xmin, xmin + xspan, length = npts))
rng = StableRNG(655)
x = x + xspan / npts * (rand(rng, npts) .- 0.5)
# select a subset randomly
idx = unique(rand(rng, collect(eachindex(x)), 20))
t = x[unique(idx)]
npts = length(t)
ut = sin.(t)
stdev = 1e-1 * maximum(ut)
u = ut + stdev * randn(rng, npts)
# data must be ordered if t̂ is not provided
idx = sortperm(t)
tₒ = t[idx]
uₒ = u[idx]

tolerance = 1e-3

@testset "Direct smoothing" begin
    # fixed with default λ = 1.0
    A = RegularizationSmooth(uₒ, tₒ; alg = :fixed)
    ans = [0.6456173647252937 0.663974701324226 0.7631218523665086 0.778654700697601 0.7489958320589535 0.7319087707475104 0.6807082599508811 0.6372557895089508 0.5832859790765743 0.5021274805916013 0.3065928203396211 0.1353332321156384 -0.3260000640060584 -0.6557906092739154 -0.9204882447932498]'
    @test isapprox(A.û, ans, rtol = tolerance)
    @test isapprox(A.(tₒ), ans, rtol = tolerance)
    # non-default d and λ
    A = RegularizationSmooth(uₒ, tₒ, 4; λ = 1e-2, alg = :fixed)
    ans = [0.19865190868740357 0.2885349151737291 0.6756699442978945 0.9165887141895426 0.9936113717653254 1.0042825002191034 0.9768118192829827 0.9184595331808411 0.8214983284892922 0.6538356458824783 0.28295521578898 0.018060767871253963 -0.5301723647977373 -0.8349855890541111 -1.1085048455468356]'
    @test isapprox(A.û, ans, rtol = tolerance)
    # GCV (default) to determine λ
    A = RegularizationSmooth(uₒ, tₒ)
    @test isapprox(A.λ, 0.12788440382063268, rtol = tolerance)
    ans = [0.21974931848164914 0.2973284508009968 0.6908546278415386 0.9300465474303226 0.9741453042418977 0.9767572556868123 0.9432951659303452 0.8889834700087442 0.804842790047182 0.6603217445567791 0.30341652659101737 0.05924456463634589 -0.5239939779242144 -0.8421768233191822 -1.107517099580091]'
    @test isapprox(A.û, ans, rtol = tolerance)
    # L-curve to determine λ
    A = RegularizationSmooth(uₒ, tₒ; alg = :L_curve)
    @test isapprox(A.λ, 0.9536286111306728, rtol = tolerance)
    ans = [
        0.6261657429321232,
        0.6470204841904836,
        0.7599270022828396,
        0.7835725197598925,
        0.7567105872094757,
        0.7404815750685363,
        0.6906841961987067,
        0.647628105931872,
        0.5937273796308717,
        0.512087780658067,
        0.3136272387739983,
        0.1392761732695201,
        -0.3312498167413961,
        -0.6673268474631847,
        -0.9370342562716745,
    ]
    @test isapprox(A.û, ans, rtol = tolerance)
    @test isapprox(A.(tₒ), ans, rtol = tolerance)
end

@testset "Smoothing with weights" begin
    # midpoint rule integration
    A = RegularizationSmooth(uₒ, tₒ, nothing, :midpoint)
    @test isapprox(A.λ, 0.10787235405005478, rtol = tolerance)
    ans = [
        0.3068904607028622,
        0.3637388879266782,
        0.6654462500501238,
        0.9056440536733456,
        0.9738150157541853,
        0.9821315604309402,
        0.9502526946446999,
        0.8953643918063283,
        0.8024431779821514,
        0.6415812230114304,
        0.2834706832220367,
        0.05281575111822609,
        -0.5333542714497277,
        -0.8406745098604134,
        -1.0983391396173634,
    ]
    @test isapprox(A.û, ans, rtol = tolerance)
    @test isapprox(A.(tₒ), ans, rtol = tolerance)
    # arbitrary weights for wls (and fixed λ, GCV not working well for some of these)
    A = RegularizationSmooth(uₒ, tₒ, nothing, collect(1:npts); λ = 1e-1, alg = :fixed)
    ans = [
        0.24640196218427968,
        0.3212059975226125,
        0.6557626475144205,
        0.9222911426465459,
        0.9913331910731215,
        1.0072241662103494,
        0.9757899817730779,
        0.935880516370941,
        0.8381074902073471,
        0.6475589703422522,
        0.2094170714475404,
        0.09102085384961625,
        -0.5640882848240228,
        -0.810519277110118,
        -1.1159124134900906,
    ]
    @test isapprox(A.û, ans, rtol = tolerance)
    @test isapprox(A.(tₒ), ans, rtol = tolerance)
    # arbitrary weights for wls and wr
    nhalf = Int(floor(npts / 2))
    wls = vcat(ones(nhalf), 10 * ones(npts - nhalf))
    wr = collect(1:(npts - 2))
    A = RegularizationSmooth(uₒ, tₒ, nothing, wls, wr; λ = 1e-1, alg = :fixed)
    ans = [
        0.21878709713242372,
        0.3118480645325099,
        0.7669822464946172,
        1.0232343854914931,
        1.0526513115274412,
        1.0469579284244412,
        0.9962426294084775,
        0.9254407155702626,
        0.8204764044515936,
        0.6514510142804217,
        0.27796896299068763,
        0.04756024028728636,
        -0.5301034620974782,
        -0.8408107101140526,
        -1.1058428573417736,
    ]
    @test isapprox(A.û, ans, rtol = tolerance)
    @test isapprox(A.(tₒ), ans, rtol = tolerance)
end

@testset "Smoothing with t̂ provided" begin
    N̂ = 20
    t̂ = collect(range(xmin, xmin + xspan, length = N̂))
    # with t̂, no weights
    A = RegularizationSmooth(u, t, t̂)
    @test isapprox(A.λ, 0.138273889585313, rtol = tolerance)
    ans = [0.21626377852882872 0.39235926952322575 0.5573848799950002 0.7072474496656729 0.8361906119247042 0.9313473799797176 0.9809844353757837 0.9750833208625507 0.9096038940899813 0.7816929736202427 0.6052694276527628 0.4015903497629387 0.1913719025253403 -0.01979786871512895 -0.23400354001942947 -0.44481229967011127 -0.6457913359497256 -0.8405146928672158 -1.0367229293434395 -1.2334090099343238]'
    @test isapprox(A.û, ans, rtol = tolerance)
    @test isapprox(A.(t̂), ans, rtol = tolerance)
    # t̂ and wls
    A = RegularizationSmooth(u, t, t̂, collect(1:npts))
    @test isapprox(A.λ, 0.26746430253489195, rtol = tolerance)
    ans = [0.3118247878815087 0.44275860852897864 0.5705834985506882 0.6979119448253899 0.8234189540704866 0.9289458273102476 0.9970803470992273 1.0071205506077525 0.9443157518324818 0.7954860908242515 0.5847385548859145 0.34813493129868633 0.1237494751337505 -0.0823517516424196 -0.28265170846635246 -0.4760833187699964 -0.6615795059853024 -0.844779821396189 -1.0341162283806349 -1.225270266213379]'
    @test isapprox(A.û, ans, rtol = tolerance)
    @test isapprox(A.(t̂), ans, rtol = tolerance)
    # t̂, wls, and wr
    nhalf = Int(floor(npts / 2))
    wls = vcat(ones(nhalf), 10 * ones(npts - nhalf))
    wr = collect(1:(N̂ - 2))
    A = RegularizationSmooth(u, t, t̂, wls, wr)
    @test isapprox(A.λ, 0.04555080890920959, rtol = tolerance)
    ans = [
        0.2799800686914433,
        0.4627548444527547,
        0.5611922868318674,
        0.6647761469309206,
        0.7910803348948329,
        0.9096001134420562,
        1.0067644979677808,
        1.0541868144785513,
        0.9889720466386331,
        0.8088479651943575,
        0.5677592185997403,
        0.31309698432269184,
        0.08587106716465115,
        -0.11476265128730469,
        -0.30749376694236485,
        -0.4942769809562725,
        -0.676806367664006,
        -0.8587832527770329,
        -1.0443430843364814,
        -1.2309001260104093,
    ]
    @test isapprox(A.û, ans, rtol = tolerance)
    @test isapprox(A.(t̂), ans, rtol = tolerance)
end
