function par0 = getParameterGuesses(parameters,ll_fun,n,lb,ub,varargin)


if nargin > 5
    options = varargin{1};
end

ind_failed_old = 1:n;
counter = 1;

par0 = nan(numel(lb),n);


while ~isempty(ind_failed_old) && counter < 5e2
    temp_n = numel(ind_failed_old);
    
    par0(:,ind_failed_old) = bsxfun(@plus,lb,bsxfun(@times,ub - lb,...
            rand(parameters.number,numel(ind_failed_old))));

    if nargin > 5
        par0(options.fixedParameters,ind_failed_old) = options.fixedParameterValues(:) * ones(1,numel(ind_failed_old));
    end
    ind_failed_new = [];
    for i = 1:temp_n
        try
            %warning off
            [ll1,dll] = ll_fun(par0(:,ind_failed_old(i)));
            [ll2] = ll_fun(par0(:,ind_failed_old(i)));
            warning on
            if isinf(ll1) || isinf(ll2)
                ind_failed_new = [ind_failed_new,ind_failed_old(i)];
            end
            if isnan(ll1) || isnan(ll2)
                ind_failed_new = [ind_failed_new,ind_failed_old(i)];
            end
        catch
            ind_failed_new  = [ind_failed_new,ind_failed_old(i)];
        end
    end
    counter = counter+1
    numel(ind_failed_new)
    ind_failed_new
    ind_failed_old = ind_failed_new;
end
end