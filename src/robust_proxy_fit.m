function [thetas] = robust_proxy_fit(PV,proxies,opt_pars,init)

if nargin<4
    init = [];
end

n_proxies = size(proxies,2);
lb = zeros(n_proxies,1);
ub = ones(n_proxies,1)*max(PV);

if isempty(init)
    theta0 = 1e-2*ones(n_proxies,1);
else
    theta0 = init;
end

fmincon_options = optimoptions('fmincon','PlotFcn',{@optimplotfval,@optimplotx});

if strcmp(opt_pars.loss_function,'epsiloninsensitive') || strcmp(opt_pars.loss_function,'eps')
    thetas = fmincon(@(thetas) epsiloninsensitive(PV,proxies,thetas,opt_pars),theta0,[],[],[],[],lb,ub,[],fmincon_options);
elseif strcmp(opt_pars.loss_function,'gaussian') || strcmp(opt_pars.loss_function,'ols')
    thetas = fmincon(@(thetas) ols(PV,proxies,thetas,opt_pars),theta0,[],[],[],[],lb,ub,[],fmincon_options);
elseif strcmp(opt_pars.loss_function,'huber')
    thetas = fmincon(@(thetas) huber(PV,proxies,thetas,opt_pars),theta0,[],[],[],[],lb,ub,[],fmincon_options);
elseif strcmp(opt_pars.loss_function,'biweight') || strcmp(opt_pars.loss_function,'tukey')
    thetas = fmincon(@(thetas) biweight(PV,proxies,thetas,opt_pars),theta0,[],[],[],[],lb,ub,[],fmincon_options);
end


end


function OF = epsiloninsensitive(PV,proxies,thetas,opt_pars)
lambda = opt_pars.lambda;
epsilon = opt_pars.epsilon;
err = (proxies*thetas -PV);
if opt_pars.rescale
    err = abs(err)./(PV.^0.5);
end
OF = mean(max([zeros(size(err)),abs(err)-epsilon],[],2).^2) + lambda*(sum(abs(thetas)));
end

function OF = ols(PV,proxies,thetas,opt_pars)
lambda = opt_pars.lambda;
err = (proxies*thetas -PV);
if opt_pars.rescale
    err = abs(err)./(PV.^0.5);
end
OF = 0.5*mean((err).^2) + lambda*(sum(abs(thetas)));
end

function OF = huber(PV,proxies,thetas,opt_pars)
lambda = opt_pars.lambda;
epsilon = opt_pars.epsilon;
err = (proxies*thetas -PV);
%err(err>0)=err(err>0)*0.25;
if opt_pars.rescale
    err = abs(err)./(PV.^0.5);
end
huber_loss = (0.5*err.^2).*(abs(err)<epsilon) + (epsilon*abs(err)-0.5*epsilon^2).*(abs(err)>epsilon);
OF = mean(huber_loss) + lambda*(sum(abs(thetas)));
end

function OF = biweight(PV,proxies,thetas,opt_pars)
lambda = opt_pars.lambda;
err = (proxies*thetas -PV);
%err(err>0)=err(err>0)*0.25;
if opt_pars.rescale
    err = abs(err)./(PV.^0.5);
end
epsilon = opt_pars.epsilon;
%     epsilon = 4.6851;
one_vect = ones(size(err));
biweight_loss = ((epsilon^2)/6)*((one_vect - (one_vect-(err/epsilon).^2).^3).*(abs(err)<epsilon) +(abs(err)>epsilon));
OF = mean(biweight_loss) + lambda*(sum(abs(thetas)));
end

