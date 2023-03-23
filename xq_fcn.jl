
###############################################
# vf for bidding
###############################################

function value_fcn_calcu(lambda, seg_num, Ne, T, c, P, eta, ed, ef, sigma)
    # generate value function samples
    vEnd = zeros(Ne,1) ;
    # use 1000 as the penalty for final discharge level
    println("  Ne = ", Ne)
    println("  ef = ", ef)
    println("  sigma = ", sigma)
    T_final = Int64(floor(ef*Ne));
    if T_final != 0
        vEnd[1:T_final] = 1000;
    end
    # initialize the value function series
    v = zeros(Ne, T+1);
    # v(1,1) is the marginal value of 0# SoC at the beginning of day 1
    # V(Ne, T) is the maringal value of 100% SoC at the beginning of the last operating day
    # update final value function
    v[:,end] = vEnd; 
    # process index
    es = (0:ed:1)';
    Ne = length(es);
    # calculate soc after charge vC = (v_t(e+P*eta))
    eC = es .+ P.*eta;
    iC = zeros(Int64, Ne);
    # round to the nearest sample
    for i = 1:Ne
        iC[i] = Int64(ceil(eC[i]./ed)+1); 
        if  iC[i] > Ne+1
            iC[i] = Ne+2;
        end
        if iC[i] < 2
            iC[i] = 1;
        end
    end
    
    # calculate soc after discharge vC = (v_t(e-P/eta))
    eD = es .- P./eta;
    iD = zeros(Int64, Ne);
    # round to the nearest sample
    for i = 1:Ne
        iD[i] = Int64(floor(eD[i]./ed)+1); 
        if  iD[i] > Ne+1
            iD[i] = Ne+2;
        end
        if iD[i] < 2
            iD[i] = 1;
        end
    end

    # Recalculate ES value fcn
    for t = T:-1:1 #     
        vi = v[:,t+1];
        vo = calc_valueF(lambda[t], sigma, c, eta, vi, iC, iD);
        v[:,t] = vo;
    end
    for i = 1:Ne
        for t = 1:T
            if v[i,t] <0
                v[i,t] = 0;
            end
        end
    end
    vp = v./eta .+ c;
    vb = v.*eta;
    ## modify here if need non-uniform segments
    vpAvgt = zeros(seg_num,T+1);
    vbAvgt = zeros(seg_num,T+1);
    NN = Int64(round((Ne-1)/seg_num));
    @printf("Type: %s \n", typeof(NN))

    for i = 1:seg_num
        id_st = Int64(round((i-1)*NN + (1)));
        id_ed = Int64(round((i-1)*NN + (NN+1)))
        for t = 1:T
        vpAvgt[i,t] = mean(vp[id_st:id_ed,t]);
        vbAvgt[i,t] = mean(vb[id_st:id_ed,t]);
        end
    end
    # projection from T slots to H hours ? I did not have it.
    # vpAvg = zeros(seg_num,H);
    # vbAvg = zeros(seg_num,H);
    # NH = T/H;
    # for j = 1:H
    #     id_seg = Int64(round((j-1)*NH + (1:(NH+1))));
    #     vpAvg[:,j] = mean(vpAvgt(:,(j-1)*NH + (1:(NH+1))),2);
    #     vbAvg[:,j] = mean(vbAvgt(:,(j-1)*NH + (1:(NH+1))),2);
    # end
    return vpAvgt, vbAvgt
end



###############################################
# VALUATION FUNCTION
###############################################
function calc_valueF(price, sigma, m_c, eta, vi, e_chg_index, e_dch_index)
    # define penalty for limit soc values
    penalty = 100000
    # expand value function, mc and eta vectors to include penalty for limit soc values
    v_ext = vcat(penalty, vi, -penalty)
    mc_ext = vcat(1,m_c,1)
    eta_ext = vcat(1,eta,1)
    
    # map charging index vector to corresponding VF
    v_chg = v_ext[e_chg_index]
    # map discharging index vector to corresponding VF
    v_dch = v_ext[e_dch_index]
    
    # calculate event vectors i.e. CDF and PDF corresponding to the deterministic prices case
    if sigma ==0 # storage participants do not consider price uncertainty
        FtEC = vi.*eta .> price
        FtCC = v_chg.*eta .> price
        FtED = (vi./eta .+ m_c).*(vi./eta .+ m_c .> 0) .> price
        FtDD = (v_dch./eta .+ m_c).*(v_dch./eta .+ m_c .> 0) .> price
        # calculate terms of next value function q_t 
        Term1 = v_chg .* FtCC
        Term2 = price .* (FtEC - FtCC) ./ eta
        Term3 = vi .* (FtED - FtEC)
        Term4 = (price .- m_c) .* eta .* (FtDD - FtED)
        Term5 = v_dch .* (1 .- FtDD)

        out = Term1 + Term2 + Term3 + Term4 + Term5;
    else    # storage participants consider price uncertainty
        mu = price; # price is the mu here
        FtEC = normcdf(vi*eta, mu, sigma); # F_t(v_t(e)*eta)
        FtCC = normcdf(v_chg*eta, mu, sigma); # F_t(v_t(e+P*eta)*eta)
        FtED = normcdf((vi/eta + m_c).*((vi/eta + m_c) > 0), mu, sigma); # F_t(v_t(e)/eta + m_c) 
        FtDD = normcdf((v_dch/eta + m_c).*((v_dch/eta + m_c) > 0), mu, sigma); # F_t(v_t(e-P/eta)/eta + m_c) 

        ftEC = normpdf(vi*eta, mu, sigma); # f_t(v_t(e)*eta)
        ftCC = normpdf(v_chg*eta, mu, sigma); # f_t(v_t(e+P*eta)*eta)
        ftED = normpdf((vi/eta + m_c).*((vi/eta + m_c) > 0), mu, sigma).*((vi/eta + m_c) > 0); # f_t(v_t(e)/eta + m_c) 
        ftDD = normpdf((v_dch/eta + m_c).*((v_dch/eta + m_c) > 0), mu, sigma).*((v_dch/eta + m_c) > 0); # f_t(v_t(e-P/eta)/eta + m_c) 
        
        Term1 = v_chg .* FtCC;
        Term2 = ( mu * (FtEC - FtCC) - sigma^2 * (ftEC - ftCC) ) / eta;
        Term3 = vi .* (FtED - FtEC);
        Term4 = ( mu * (FtDD - FtED) - sigma^2 * (ftDD - ftED) ) * eta;
        Term5 = - m_c * eta * (FtDD - FtED);
        Term6 = v_dch .* (1-FtDD);

        out = Term1 + Term2 + Term3 + Term4 + Term5 + Term6;
    end

    
    return out
end


function normcdf(x, mu, sigma)
    d = Normal(mu,sigma); # normal distribution with mu and sigma
    out = zeros(length(x))
    for i = 1: length(x)
        out[i] = cdf(d,x[i])
    end
    return out
end