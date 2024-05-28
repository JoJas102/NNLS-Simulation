function [simResults,statistic] = evaluation(resultNNLS, resultNLLS, resultNLLSraw, dAUC, fAUC, fdInput, np, snrIn, m, fd_offset)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Preperations and calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % Read variables
    fIn        = fdInput(1,1:np+1)';
    dIn        = fdInput(2,1:np+1)';
    dNNLS(:,:) = resultNNLS(:,2,:) - fd_offset(1,:)';
    fNNLS(:,:) = resultNNLS(:,4,:) - fd_offset(2,:)';
    dNLLS(:,:) = resultNLLS(1:np+1,:) - fd_offset(1,:)';
    fNLLS(:,:) = resultNLLS(np+2:end-1,:) - fd_offset(2,:)';
    dNLLSraw(:,:) = resultNLLSraw(1:np+1,:) - fd_offset(1,:)';
    fNLLSraw(:,:) = resultNLLSraw(np+2:end-1,:) - fd_offset(2,:)';
    iter       = length(dNNLS(1,:));
    
    % Deviation to ground truth
    dDiffNNLS     = dNNLS-dIn;
    dDiffNLLS     = dNLLS-dIn;
    dDiffNLLSraw  = dNLLSraw-dIn;
    dDiffNNLSp    = abs(mean(dNNLS,2)-dIn)./dIn*100;
    dDiffNLLSp    = abs(mean(dNLLS,2)-dIn)./dIn*100;
    fDiffNNLS     = (fNNLS-fIn);
    fDiffNLLS     = (fNLLS-fIn);
    fDiffNLLSraw  = (fNLLSraw-fIn);
    fDiffNNLSmean = abs(mean(fDiffNNLS,2));
    fDiffNLLSmean = abs(mean(fDiffNLLS,2));

    % Simulation results table NEW
    statistic   = cat(3,dDiffNNLS, fDiffNNLS, dDiffNLLS, fDiffNLLS);
    methods     = ["NNLS" "NNLS_AUC" "NLLS*" "NLLS"];
    iter_simu   = repelem(1:iter,length(methods))';
    method      = repmat(methods,1,iter)'; 
    peaks       = repelem(sum(dNNLS~=0,1),4)';
    dAUC(:,4)   = zeros(iter,1);
    dAUC        = dAUC';
    dDiffAUC    = dAUC-dIn;
    fAUC(:,4)   = zeros(iter,1);
    fAUC        = fAUC';
    fDiffAUC    = fAUC-fIn;

    l = length(iter_simu);
    dResults(1:4:l,:) = dDiffNNLS';
    dResults(2:4:l,:) = dDiffAUC';
    dResults(3:4:l,:) = dDiffNLLS';
    dResults(4:4:l,:) = dDiffNLLSraw';
    fResults(1:4:l,:) = fDiffNNLS';
    fResults(2:4:l,:) = fDiffAUC';
    fResults(3:4:l,:) = fDiffNLLS';
    fResults(4:4:l,:) = fDiffNLLSraw';

    simResults  = table(method, iter_simu, peaks, dResults(:,1), dResults(:,2), dResults(:,3), dResults(:,4), fResults(:,1), fResults(:,2), fResults(:,3), fResults(:,4), 'VariableNames',["method","iter","peaks","d_slow","d_inter","d_fast","d_4","f_slow","f_inter","f_fast", "f_4"]);
    simResults(end+1,:) = {"gT" dIn(1) dIn(2) dIn(3) dIn(4) fIn(1) fIn(2) fIn(3) fIn(4) NaN NaN};
    simResults(end+1,:) = {"SNR, bins, iter" snrIn m iter NaN NaN NaN NaN NaN NaN NaN};
    
    % Farben
    HHU = [0/256,106/256,179/256];
    HHUeisblau = [181/256,203/256,214/256];
    HHUhellblau = [89/256,143/256,179/256];
    HHUdunkelblau = [0,57/256,100/256];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting statistic boxcharts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Statistics as table struct
    method = [repelem("NNLS",(np+1)*iter) repelem("NLLS",(np+1)*iter)]';
    method = categorical(method);
    compartment = repmat([repelem(1:np+1,iter)]',2,1);
    dValue = [reshape(dNNLS.',1,[]), reshape(dNLLS.',1,[])]';
    fValue = [reshape(fNNLS.',1,[]), reshape(fNLLS.',1,[])]';
    iterc = repmat(1:iter,1,2*(np+1))';
    results = table(method,compartment,dValue,fValue,iterc);
    idx = results.compartment<=3; % 3 peaks only 
    
    % diff coeff boxplot
    figure(2)
    b = boxchart(results.compartment(idx),results.dValue(idx),'GroupByColor',results.method(idx),'linew',1);
    
    % Color and design
    b(1).BoxFaceColor = HHUeisblau;
    b(1).BoxLineColor = HHUdunkelblau;
    b(1).WhiskerLineColor = HHUdunkelblau;
    b(1).MarkerColor = HHUhellblau;
    b(1).BoxFaceAlpha = 0.4;
    b(2).BoxFaceColor = HHUdunkelblau;
    b(2).MarkerColor = HHUdunkelblau;
    b(2).BoxFaceAlpha = 0.6;
    set(gca,'FontSize',14)
    set(gca,'YScale','log')
    set(gca,'xtick',[1,2,3],'xticklabel',{'D_{slow}';'D_{inter}';'D_{fast}'})
    ylabel('D (mm^2/s)');
    title('Diffusion coefficients','FontSize',16)
    xlim([.2,3.5])
    ylim([5e-4 4e-1])

    % Display gT values 
    hold on
    plot([0,1.35], [dIn(1),dIn(1)], ':', 'Color', 'k');
    text(.25,dIn(1)+0.0005,'D^{in}_{slow}','FontSize',13)
    plot([0,2.35], [dIn(2),dIn(2)], ':','Color', 'k');
    text(.25,dIn(2)+0.003,'D^{in}_{inter}','FontSize',13)
    if fIn(3) ~= 0 %! adjust manually if < 3 compartments including fast component are simulated
       plot([0,3.35], [dIn(3),dIn(3)], ':','Color', 'k');
        text(.25,dIn(3)+0.09,'D^{in}_{fast}','FontSize',13)
    end
        legend(b,'NLLS*','NNLS','Location','southeast')
    hold off

    % Vol frac boxplot
    figure(3)
    b2 = boxchart(results.compartment(idx),results.fValue(idx),'GroupByColor',results.method(idx),'linew',1);
    
    % Color and design
    b2(1).BoxFaceColor = HHUeisblau;
    b2(1).BoxLineColor = HHUdunkelblau;
    b2(1).WhiskerLineColor = HHUdunkelblau;
    b2(1).MarkerColor = HHUhellblau;
    b2(1).BoxFaceAlpha = 0.4;
    b2(2).BoxFaceColor = HHUdunkelblau;
    b2(2).MarkerColor = HHUdunkelblau;
    b2(2).BoxFaceAlpha = 0.6;
    set(gca,'FontSize',14)
    set(gca,'xtick',[1,2,3],'xticklabel',{'f_{slow}';'f_{inter}';'f_{fast}'})
    ylabel('Volume fraction (%)');
    title('Volume fractions','FontSize',16)
    xlim([.2,3.5])
    ylim([0 80])

    % Display gT values 
    hold on
    plot([0,1.35], [fIn(1),fIn(1)], ':', 'Color', 'k');
    text(.25,fIn(1)+5,'D^{in}_{slow}','FontSize',13)
    plot([0,2.35], [fIn(2),fIn(2)], ':','Color', 'k');
    text(.25,fIn(2)+5,'D^{in}_{inter}','FontSize',13)
    if fIn(3) ~= 0 %! adjust manually if < 3 compartments including fast component are simulated
       plot([0,3.35], [fIn(3),fIn(3)], ':','Color', 'k');
        text(.25,fIn(3)+5,'D^{in}_{fast}','FontSize',13)
    end
        legend(b2,'NLLS*','NNLS','Location','northeast')
    hold off

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting statistic bar graphs 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    figure(4)

    % Avg vol frac bar graph
    subplot(2,1,1);
    barsMean = [mean(fNLLS(1,:)) mean(fNNLS(1,:)); mean(fNLLS(2,:)) mean(fNNLS(2,:)); mean(fNLLS(3,:)) mean(fNNLS(3,:))];
    barsDiff = [fDiffNLLSmean(1) fDiffNNLSmean(1); fDiffNLLSmean(2) fDiffNNLSmean(2); fDiffNLLSmean(3) fDiffNNLSmean(3)];
    
    yyaxis left
    plot([0,1.25], [fIn(1),fIn(1)], ':', 'Color', 'k');
    text(0.05,fIn(1)+2,'f^{in}_{slow}','FontSize',12)
    hold on
    plot([0,2.25], [fIn(2),fIn(2)], ':','Color', 'k');
    text(0.05,fIn(2)+2,'f^{in}_{inter}','FontSize',12)
    if fIn(3) ~= 0
        plot([0,3.25], [fIn(3),fIn(3)], ':','Color', 'k');
        text(0.05,fIn(3)+2,'f^{in}_{fast}','FontSize',12)
    end
    b = bar(barsMean,'FaceColor','flat');
    set(gca,'FontSize',14)
    set(gca,'YColor','k')
    b(1).CData = HHUhellblau;
    b(2).CData = HHUdunkelblau;
    ylim([0 70]);
    ylabel('Volume fraction f (%)');
    hold off
    
    yyaxis right
    b2 = bar(barsDiff);
    set(gca,'FontSize',14)
    ylim([0 70]);
    set(gca,'xtick',[1:3],'xticklabel',{'f_{slow}';'f_{inter}';'f_{fast}'})
    ylabel('Deviation');
    title('(a) Average volume fractions','FontSize',16)
    legend([b b2(1)],'NLLS*','NNLS','Deviation','Location', 'northeast')
    
    
    % Avg diff coeff bar graph 
    subplot(2,1,2); 
    barsMean = [mean(dNLLS(1,:))      mean(dNNLS(1,:))     ; mean(dNLLS(2,:))      mean(dNNLS(2,:))     ; mean(dNLLS(3,:))      mean(dNNLS(3,:))     ];
    barsDiff = [mean(dDiffNLLSp(1,:)) mean(dDiffNNLSp(1,:)); mean(dDiffNLLSp(2,:)) mean(dDiffNNLSp(2,:)); mean(dDiffNLLSp(3,:)) mean(dDiffNNLSp(3,:))];
    
    yyaxis left
    semilogy([0,1.25], [dIn(1),dIn(1)], ':', 'Color', 'k');
    text(0.05,dIn(1)+0.0002,'D^{in}_{slow}','FontSize',12)
    hold on
    semilogy([0,2.25], [dIn(2),dIn(2)], ':','Color', 'k');
    text(0.05,dIn(2)+0.001,'D^{in}_{inter}','FontSize',12)
    if fIn(3) ~= 0
        semilogy([0,3.25], [dIn(3),dIn(3)], ':','Color', 'k');
        text(0.05,dIn(3)+0.03,'D^{in}_{fast}','FontSize',12)
    end
    b = bar(barsMean,'FaceColor','flat');
    set(gca,'YScale','log')
    set(gca,'YColor','k')
    set(gca,'FontSize',14)
    b(1).CData = HHUhellblau;
    b(2).CData = HHUdunkelblau;
    ylabel('D (mm^2/s)');
    ylim([1e-4 1e-0]);
    hold off
    
    yyaxis right
    bar(barsDiff);
    set(gca,'FontSize',14)
    yticks([50 100])
    yticklabels({'0.5','1'})
    ylim([0 400]);
    set(gca,'xtick',[1:3],'xticklabel',{'D_{slow}';'D_{inter}';'D_{fast}'})
    ylabel('Relative deviation');
    title('(b) Average diffusion coefficients','FontSize',16)

end
