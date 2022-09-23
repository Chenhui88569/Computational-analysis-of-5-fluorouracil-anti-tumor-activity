clear
close all
T =  5000;

%init_para = [2.5101e-4,1.1623e1,4.4318e-5 ,0.2,0.2];
%init_para = [1253.66895622927,0.692746650820480,0.781733880493165,0.576082730498741,0.000225278706464815,2.36269456175158,3.16948420762964e-06, 0.2,0.2];
%init_para_col = [581.367899173836,1.41924572740766,23.9608321851255,161.434899824332,0.000428266684476464,1.42037832167050,2.79732808262959e-05, 0.2 , 0.2];
init_para_col = [497.378960533090,7.62618240795310,0.673867033252338,5.14201932190539,7.66432919815531e-05,4.84137784103185,1.15640279846837e-05;581.367899173836,1.41924572740766,23.9608321851255,161.434899824332,0.000428266684476464,1.42037832167050,2.79732808262959e-05;893.445064075148,0.00376576201660318,0.175114274349874,0.238015233076624,0.000212517367465980,1.24679407958498,5.05055761679821e-05;893.445064075148,0.00376576201660318,0.175114274349874,0.238015233076624,0.000251011503209797,1.72303673413627,4.43176610699279e-05;1253.66895622927,0.692746650820480,0.781733880493165,0.576082730498741,0.000225278706464815,2.36269456175158,3.16948420762964e-06;893.445064075148,0.00376576201660318,0.175114274349874,0.238015233076624,0.000283709997955965,11.6233335469169,3.42893935816365e-05;893.445064075148,0.00376576201660318,0.175114274349874,0.238015233076624,0.000150078500719029,2.64192539865509,1.77069663982568e-05;1115.96765386021,87.8831577334009,0.506721148813672,62.8754391498571,0.000231632570338714,8.98484143967960,4.03604272204016e-06;1300.24459206718,72.2011253506358,1.52495680944258,49.6560305403402,0.000273700742242568,7.56442192772004,5.51313396085374e-06];
Model_Index = 7;  LoadData = 1;
if ismember( Model_Index, [4,6,7])
    para_label = [  "\lambda_{g}",   "P_{max}",   "\lambda_{d}",  "\sigma_{treated}^2","\sigma_{control}^2" ];
    init_para = [init_para_col(Model_Index,5:7),0.2, 0.2];
else
    para_label = [ "T_{d,Tv}", "IC_{50}",   "E_{max , damage }",   "EC_{50,damage}",   "\lambda_{g}",   "P_{max}",   "\lambda_{d}",  "\sigma_{treated}^2","\sigma_{control}^2" ];
      init_para = [init_para_col(Model_Index,:),0.2,0.2];

end
init_para = init_para';

if LoadData  == 1
    num_para = length(para_label);
    load(['MCMC_TGI_' num2str(Model_Index) '/' 'Para_TGI_col_uni_' num2str(Model_Index) '.mat'])
    Para_col_hist = Para_col;
    Para_col = zeros( T,    num_para ) ;
    Para_col(1,:) =Para_col_hist(end,1:num_para);%Para_col_hist(end,1:num_para); %init_para ;% Para_col_hist(end,1:num_para
    
    S_curr  = S_curr( 1: num_para, 1:num_para );
    myfun = @Target_fun_LE_TGI_Uniform;
    outputfile_dir = ['MCMC_TGI_'  num2str(Model_Index)  '/'];
    [Para_col,  acc_rate_col , S_curr] = MyRAM(  T, num_para, Para_col , S_curr,myfun,para_label, outputfile_dir ,Model_Index );
    Para_col = [Para_col_hist; Para_col(2:end,:)];
    save( [ outputfile_dir   'Para_TGI_col_uni_' num2str(Model_Index) '.mat']  ,'Para_col'  , 'acc_rate_col',  'S_curr' ,'-v7.3');
else
    num_para = length(para_label);
    Para_col = zeros( T,    num_para ) ;
    Para_col(1,:) = init_para;%Para_col_hist(end,1:num_para); %init_para ;% Para_col_hist(end,1:num_para);
    
    if  ~ismember( Model_Index, [4,6,7])
        S_curr = zeros(num_para);
        for i = 1: num_para
            S_curr(i,1:i) =  abs(randn)*0.01;
            if  i == 1
                S_curr(i,1:i) =  abs(randn);
            end
            if  i == 5 ||  i == 7
                S_curr(i,1:i) =  abs(randn)*1e-4;
            end
%             if  i == 2
%                 S_curr(i,1:i) =  abs(randn);
%             end
            if  Model_Index== 3 && i ==2
                S_curr(i,1:i) =  abs(randn)*0.001;
            end
        end
    else
        S_curr = zeros(num_para);
        for i = 1: num_para
            S_curr(i,1:i) =  abs(randn)*0.001;
            if  i == 1 ||  i == 3
                S_curr(i,1:i) =  abs(randn)*1e-5;
            end
        end
    end
    S_curr  = S_curr( 1: num_para, 1:num_para );
    myfun = @Target_fun_LE_TGI_Uniform;
    outputfile_dir = ['MCMC_TGI_'  num2str(Model_Index)  '/'];
    [Para_col,  acc_rate_col , S_curr] = MyRAM(  T, num_para, Para_col , S_curr,myfun,para_label, outputfile_dir ,Model_Index );
    save( [ outputfile_dir   'Para_TGI_col_uni_' num2str(Model_Index) '.mat']  ,'Para_col'  , 'acc_rate_col',  'S_curr' ,'-v7.3');
end
