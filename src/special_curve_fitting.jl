# SigmoidFit
function SigmoidFit(u,t,p=zeros(2))
  m = (x, p) -> @. p[1]/(1+exp(x-p[2]))
  c_f = curve_fit(m,t,u,p)
  Curvefit{true}(u,t,m,c_f)
end

# HillFit
function HillFit(u,t,p=zeros(2),n=nothing)
  if n == nothing
    push!(p,rand(1)[1])
    m = (x, p) -> @. p[1] * inv(1.0 + p[2]/(x^p[3]))
  else
    m = (x, p) -> @. p[1] * inv(1.0 + p[2]/(x^n))
  end
  c_f = curve_fit(m,t,u,p)
  Curvefit{true}(u,t,m,c_f)
end

# WeibullFit
function WeibullFit(u,t,p=zeros(3))
  m = weibull_fun
  c_f = curve_fit(m,t,u,p)
  Curvefit{true}(u,t,m,c_f)
end

# OneCompartmentPKFit
function OneCompartmentPKFit(u,t,p=zeros(5))
  m = (x,p) -> @. p[1] * p[2] * p[3] * (exp(-p[4]*x) - exp(-p[1]*x)) / (p[5] * (p[1] - p[4]))
  c_f = curve_fit(m,t,u,p)
  Curvefit{true}(u,t,m,c_f)
end
