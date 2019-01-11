# SigmoidFit
function SigmoidFit(u,t,p=zeros(2))
  m = (x, p) -> @. p[1]/(1+exp(x-p[2]))
  ####
end

# HillFit
function HillFit(u,t,p=zeros(2),n=nothing)
  if n == nothing
    push!(p,rand(1)[1])
    m = (x, p) -> @. p[1] * inv(1.0 + p[2]/(x^p[3]))
  else
    m = (x, p) -> @. p[1] * inv(1.0 + p[2]/(x^n))
  end
  ###
end

# WeibullFit
function WeibullFit(u,t,p=zeros(3))
  m = weibull_fun
  ###
end

# OneCompartmentPKFit
function OneCompartmentPKFit(u,t,p=zeros(5))
  m = (x,p) -> @. p[1] * p[2] * p[3] * (exp(-p[4]*x) - exp(-p[1]*x)) / (p[5] * (p[1] - p[4]))
  ###
end
