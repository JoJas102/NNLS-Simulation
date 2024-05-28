function plotSimu(inputSimu, resultSimu, dNNLS, dNLLS, fNLLS, numberOfIter)

    A = size(inputSimu);
    b = inputSimu(1,:,1);
    dIn = inputSimu(3,1:3,1);
    signal(:,:) = inputSimu(5,:,:);
    DValues =  resultSimu(1,:,1);
    s(:,:) =  resultSimu(2,:,:);
    s2(:,:) =  resultSimu(3,:,:);
    DBasis = exp(-kron(b', DValues));
    sNLLS = exp(-kron(b, dNLLS(1)))*fNLLS(1) + exp(-kron(b, dNLLS(2)))*fNLLS(2) + exp(-kron(b, dNLLS(3)))*fNLLS(3);
    
    figure(1)
    if numberOfIter == 1
        singlePlot(signal, DBasis, b, DValues, dIn, s, dNNLS, s2, dNLLS, sNLLS);
    else
        multiPlot(A, s, s2, DValues, dIn);
    end
end

function singlePlot(signal, DBasis, b, DValues, dIn, s, dNNLS, s2, dNLLS, sNLLS)
    close all
    np = length(dIn);
    
    % compare signal and fitting result
    subplot(3,1,1);
    y_recon = DBasis*s';
    y_recon2 = DBasis*s2';
    plot( b , signal , 'ko' , b , y_recon , 'b-', b, y_recon2 , 'r-', b, sNLLS, 'm-') 
    ylim([0 101]);
    New_YTickLabel = get(gca,'ytick');
    Newer_YTickLabel = New_YTickLabel./100; 
    set(gca,'YTickLabel',Newer_YTickLabel);
    title('(a) Signal decay fitting');
    xlabel('b-values (s/mm^2)');
    ylabel('Signal amplitude');
    legend('data (a), D^{in}_i (c)','lsqnonneg','NNLSreg', 'NLLS fit')

    subplot(3,1,2);
    plot(b, signal-signal, 'k-', b, signal'-y_recon , 'b-', b, signal'-y_recon2 , 'r-', b, signal-sNLLS, 'm-')
    title('(b) Residual plot');
    a=[abs(mean(signal'-y_recon)) abs(mean(signal'-y_recon2)) abs(mean(signal-sNLLS))];
    xlabel(['Avg residual: lsqnonneg = ' num2str(a(1)) ' | NNLSreg = ' num2str(a(2)) ' | NLLS fit = ' num2str(a(3))]);
    
    subplot(3,1,3);
    yyaxis left
    semilogx(DValues,s);
    y=ylim;
    ytext=y(2)-2;
    text(0.9*1e-3,ytext,'D_{slow}')
    text(4*1e-3,ytext,'D_{inter}')
    text(90*1e-3,ytext,'D_{fast}')
    yr=y(2)-0.5;
    h(1)=rectangle('Position',[2*1e-3 0.2 8*1e-3 yr],'FaceColor',[.9 .9 .9],'EdgeColor','none'); % [Wong2019], [Stabinska2020], [MA E.Wilken] & [Periquito2021]
    uistack(h(1),'bottom')
    title('(c) Results');
    xlabel('D (mm^2/s) (\cdot 10^{-3})');
    ylabel('Amplitude');
    yyaxis right
    semilogx(DValues,s2,'r-', dIn, ones(np), 'ko', dNNLS(1:3), ones(np).*3/4, 'ro', dNLLS(1:3), ones(np)./2, 'mo');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'r';
end

function multiPlot(A, s, s2, DValues, dIn)
    close all
    
    hold on;
    view(3)                                                  % set dim=3 for 3-D plot
    surf(1:A(3),DValues,s2);                                 % s (reg) or s2 (noReg)
    plot3([0,A(3)+1], [dIn(1),dIn(1)] ,[0,0], 'Color', 'r'); % lines with ground truth values
    plot3([0,A(3)+1], [dIn(2),dIn(2)] ,[0,0], 'Color', 'r');
    plot3([0,A(3)+1], [dIn(3),dIn(3)] ,[0,0], 'Color', 'r');
    xlim([0 A(3)+1])
    ylim([0.0007 0.3])
    grid on
    hold off;
    set(gca,'yscale','log')
    if A(3) >= 100
        shading interp                                      % smoothing 3D plot for high number of simus
    end
    map = load('colormapHHU.mat');                          % use HHU colors
    colormap(map.ColormapHHU)
    
    xlabel('Number of simulation');
    ylabel('Diffusion coefficient (mm^2/s)');
    zlabel('Amplitude');    
end