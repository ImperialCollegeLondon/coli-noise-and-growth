function mu_dblgs_pr_hr = alpha2mu ( alpha_per_min )

mu_dblgs_pr_hr = alpha_per_min .* 60 ./ log(2) ;

end

