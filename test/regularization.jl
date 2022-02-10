import StableRNGs: StableRNG

# create scattered data
npts = 50
xmin = 0.0
xspan = 3/2*π
x = collect(range(xmin, xmin+xspan, length=npts))
rng = StableRNG(655)
x = x + xspan/npts*(rand(rng,npts) .- 0.5)
# select a subset randomly
idx = unique(rand(rng, collect(eachindex(x)), 20))
t = x[unique(idx)]
npts = length(t)
ut = sin.(t)
stdev = 1e-1*maximum(ut)
u = ut + stdev*randn(rng, npts)
# data must be ordered if t̂ is not provided
idx = sortperm(t)
tₒ = t[idx]
uₒ = u[idx]

tolerance = 1e-6

@testset "Direct smoothing" begin
    # fixed with default λ = 1.0
    A = RegularizationSmooth(uₒ,tₒ; alg=:fixed)
    ans = [0.6456173647252937 0.663974701324226 0.7631218523665086 0.778654700697601 0.7489958320589535 0.7319087707475104 0.6807082599508811 0.6372557895089508 0.5832859790765743 0.5021274805916013 0.3065928203396211 0.1353332321156384 -0.3260000640060584 -0.6557906092739154 -0.9204882447932498]'
    @test isapprox(A.û, ans, rtol=tolerance)
    # non-default d and λ
    A = RegularizationSmooth(uₒ,tₒ, 4; λ=1e-2, alg=:fixed)
    ans = [0.19865190868740357 0.2885349151737291 0.6756699442978945 0.9165887141895426 0.9936113717653254 1.0042825002191034 0.9768118192829827 0.9184595331808411 0.8214983284892922 0.6538356458824783 0.28295521578898 0.018060767871253963 -0.5301723647977373 -0.8349855890541111 -1.1085048455468356]'
    @test isapprox(A.û, ans, rtol=tolerance)
    # GCV (default) to determine λ
    A = RegularizationSmooth(uₒ,tₒ)
    @test isapprox(A.λ, 0.12788440382063268, rtol=tolerance)
    ans = [0.21974931848164914 0.2973284508009968 0.6908546278415386 0.9300465474303226 0.9741453042418977 0.9767572556868123 0.9432951659303452 0.8889834700087442 0.804842790047182 0.6603217445567791 0.30341652659101737 0.05924456463634589 -0.5239939779242144 -0.8421768233191822 -1.107517099580091]'
    @test isapprox(A.û, ans, rtol=tolerance)
    # L-curve to determine λ
    A = RegularizationSmooth(uₒ,tₒ; alg=:L_curve)
    @test isapprox(A.λ, 0.3339750913481082, rtol=tolerance)
    ans = [0.3038273346447746 0.3695574797289306 0.7261740581164446 0.8945169297047363 0.9100773683859941 0.9048101771171335 0.8663722295718835 0.8188893880450858 0.7520030529566242 0.641893306140736 0.35912131844587486 0.12553463033078033 -0.45446129364457777 -0.8246966349034558 -1.118321479210823]'
    @test isapprox(A.û, ans, rtol=tolerance)
end

@testset "Smoothing with weights" begin
    # midpoint rule integration
    A = RegularizationSmooth(uₒ,tₒ, nothing, :midpoint)
    @test isapprox(A.λ, 0.10768759624567165, rtol=tolerance)
    ans = [0.30685734533678666 0.3637106676306693 0.6653926790689593 0.9056721139798669 0.9738743511381565 0.9822029085679189 0.9503276381232628 0.8954395464903986 0.8024964523151614 0.6415834595513483 0.2833951659794013 0.052823783087081515 -0.5334066052968335 -0.8406485796257146 -1.0983153546407178]'
    @test isapprox(A.û, ans, rtol=tolerance)
    # arbitrary weights for wls
    A = RegularizationSmooth(uₒ,tₒ, nothing, collect(1:npts))
    @test isapprox(A.λ, 0.1447077651169621, rtol=tolerance)
    ans = [0.25185799867326253 0.32130333169345704 0.6679088195845279 0.927260505998722 0.9914594632495953 1.0028090695917116 0.9749355187531624 0.9254679051053563 0.8291064785367862 0.6493316082309926 0.2324329561930743 0.07635651700236243 -0.5564454824473795 -0.817911933888329 -1.1130193817730065]'
    @test isapprox(A.û, ans, rtol=tolerance)
    # arbitrary weights for wls and wr
    nhalf = Int(floor(npts/2))
    wls = vcat(ones(nhalf),10*ones(npts-nhalf))
    wr = collect(1:npts-2)
    A = RegularizationSmooth(uₒ,tₒ, nothing, wls, wr)
    @test isapprox(A.λ, 0.01920769986569919, rtol=tolerance)
    ans = [0.16161677980889166 0.34630920002600707 0.6469896059918103 0.9219255811754361 1.0001382485019457 1.0260666135098766 1.0278670301547874 0.9823839714037528 0.8512856706983483 0.6416080323126064 0.1892473613483678 0.10264529505906245 -0.5680977857976282 -0.8075007645180948 -1.1168524120952026]'
    @test isapprox(A.û, ans, rtol=tolerance)
end

@testset "Smoothing with t̂ provided" begin  ##### TBD
    N̂ = 20
    t̂ = collect(range(xmin, xmin+xspan, length=N̂))
    # with t̂, no weights
    A = RegularizationSmooth(u,t, t̂)
    @test isapprox(A.λ, 0.138273889585313, rtol=tolerance)
    ans = [0.21626377852882872 0.39235926952322575 0.5573848799950002 0.7072474496656729 0.8361906119247042 0.9313473799797176 0.9809844353757837 0.9750833208625507 0.9096038940899813 0.7816929736202427 0.6052694276527628 0.4015903497629387 0.1913719025253403 -0.01979786871512895 -0.23400354001942947 -0.44481229967011127 -0.6457913359497256 -0.8405146928672158 -1.0367229293434395 -1.2334090099343238]'
    @test isapprox(A.û, ans, rtol=tolerance)
    # t̂ and wls
    A = RegularizationSmooth(u,t, t̂, collect(1:npts))
    @test isapprox(A.λ, 0.26746430253489195, rtol=tolerance)
    ans = [0.3118247878815087 0.44275860852897864 0.5705834985506882 0.6979119448253899 0.8234189540704866 0.9289458273102476 0.9970803470992273 1.0071205506077525 0.9443157518324818 0.7954860908242515 0.5847385548859145 0.34813493129868633 0.1237494751337505 -0.0823517516424196 -0.28265170846635246 -0.4760833187699964 -0.6615795059853024 -0.844779821396189 -1.0341162283806349 -1.225270266213379]'
    @test isapprox(A.û, ans, rtol=tolerance)
    # t̂, wls, and wr
    nhalf = Int(floor(npts/2))
    wls = vcat(ones(nhalf),10*ones(npts-nhalf))
    wr = collect(1:N̂-2)
    A = RegularizationSmooth(u,t, t̂, wls, wr)
    @test isapprox(A.λ, 0.04542202002453003, rtol=tolerance)
    ans = [0.27982685119190176 0.4629485181120683 0.5612727224498559 0.6647076034312198 0.790945652904372 0.9094680579558058 1.0067349357353885 1.0543162288524837 0.9890994274368918 0.8088681588493614 0.5676793701053269 0.3129632437767546 0.08578111296361303 -0.11476994629546472 -0.30745805411723526 -0.49422228232466525 -0.6767431598309912 -0.858727390770621 -1.0443219211010253 -1.2309211555227055]'
    @test isapprox(A.û, ans, rtol=tolerance)
end