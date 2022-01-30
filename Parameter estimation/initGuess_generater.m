function para_initGuess_modified  = initGuess_generater(init_para,population_size,LB,UB,para_totalnum, Largerpertubation, LargerParanum)
%INITGUESS_GENERATER 
%  
%para_lb = 0.75;
%para_ub = 1.25; 
coeff_lower_bounds = LB*init_para;
coeff_upper_bounds = UB*init_para;
 

initial_population = zeros(population_size, para_totalnum);%matrix
for  index_population  =1:population_size
    initial_population(index_population,:)=(coeff_upper_bounds-coeff_lower_bounds).*rand(1, para_totalnum)+coeff_lower_bounds;
    %initial_population(index_population,:)=(coeff_upper_bounds-init_para).*rand(1, para_totalnum)+(coeff_lower_bounds-init_para).*rand(1, para_totalnum);
    
end

if Largerpertubation ==1
    for  index_population  =1:population_size
        coeff_lower_bounds_larger = LB*init_para(1:LargerParanum);%*0.5
        coeff_upper_bounds_larger = UB*init_para(1:LargerParanum);%*1.5
        initial_population( index_population, 1:LargerParanum )=(  coeff_upper_bounds_larger...
            -  coeff_lower_bounds_larger).*rand(1, LargerParanum)+coeff_lower_bounds_larger;
    end
end
para_initGuess_modified = initial_population;
end
