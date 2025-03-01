%% show details 
tempProt = find(Protype<3);

clf
fig2 = 450;
fig3 = 1000;
colormapw = brewermap(fig2,'Blues'); % or fig3 as you want

selgoodrepeats = [md1,md2,md3,drawlist(15,1:2)]; % fig2
selgoodrepeats = [drawlist(27,[1,2,4,5]+30),drawlist(15,[1,2,12])]; %fig 3

ordertemp = bestorder(StrainSumProm(Protype<3,selgoodrepeats)');
selgoodrepeats = selgoodrepeats(ordertemp);

%% plot
   
data = StrainSumProm(Protype<3, selgoodrepeats); 
numVars = size(selgoodrepeats, 2);
varNames = strrep(strains.strain(selgoodrepeats),'_',' ');
% Compute correlation matrix
corrMatrix = corr(data);

% Set colormap
cmap = [repmat([1 1 1],1000-height(colormapw),1);colormapw]; % You can use other colormaps like 'parula', 'hot', etc.

% Normalize correlation values to fit colormap index (range [1, size(cmap,1)])
minCorr = 0; maxCorr = 1; % Correlation values range from 0 to 1
cn=0;
for i = 1:numVars
    for j = 1:numVars
        cn=cn+1;
        subplot(numVars,numVars,cn)

        % Label column names (X-axis)
        if i == 1
          title(varNames{j},'FontSize',6)
        end
        
        if i == j % between repeats
            % Show repeats correlation
                corrValue = strains.max_corr(selgoodrepeats(i));
        
                % Set background color based on correlation value
                hold on
                % plot repleats
                % select top 300 targets
                [~,tempsel] = maxk(JingsumProm(Protype<3,strains.bestrps{selgoodrepeats(i)}),300);
                scatter(JingsumProm(tempProt(unique(tempsel)),strains.bestrps{selgoodrepeats(i)}(1)), ...
                JingsumProm(tempProt(unique(tempsel)),strains.bestrps{selgoodrepeats(i)}(2)),10,'white','filled','MarkerFaceAlpha',0.5);
                bgColor = interp1(linspace(minCorr, maxCorr, size(cmap,1)), cmap, corrValue);
                getylim = ylim;
                getxlim = xlim;
                % background color: corr. values
                r=rectangle('Position', [0 0 getxlim(2) getylim(2)], 'FaceColor', bgColor, 'EdgeColor', 'k');
              
                text(0.2*getxlim(2), 0.9*getylim(2), sprintf('r=%.2f', corrValue), 'FontSize', 8, ...
                 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'w');

                temp = [max(getxlim(1),getylim(1)),min(getylim(2),getxlim(2))];
                plot(temp,temp,'--','color',[.75 .75 .75])
                xlim(getxlim)
                ylim(getylim)  
                uistack(r,'bottom')
                ax = gca;
                ax.YAxis.Exponent = 4;
                ax.XAxis.Exponent = 4;

        else % between mutants
            corrValue = corrMatrix(i, j); % Extract correlation
            getylim =[0, max(max(data(:,[i,j])))];
            % Set background color based on correlation value
            hold on;
            if i<j
            % Scatter plot for each pair of variables
                [~,tempsel] = maxk(data(:,[i,j]),300);
                scatter(data((unique(tempsel)), j), data((unique(tempsel)), i), ...
                    10, 'filled','MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',0.5, ...
                    'MarkerEdgeColor','k','LineWidth',0.25);                 
                getxlim =xlim;
                getylim =ylim;
                text(0.2*getxlim(2), 0.9*getylim(2), sprintf('r=%.2f', corrValue), 'FontSize', 8, ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'm');
                %             set(ax,'color','none')
                %                 text([0 getxlim(2)],[getylim(1) getylim(1)],num2str(getxlim'))
                temp = [max(getxlim(1),getylim(1)),min(getylim(2),getxlim(2))];
                plot(temp,temp,'--','color',[.75 .75 .75])
                
                bgColor = interp1(linspace(minCorr, maxCorr, size(cmap,1)), cmap, corrValue);
                r= rectangle('Position', [getxlim(1) getylim(1) getxlim(2)-getxlim(1) getylim(2)-getylim(1)], ...
                    'FaceColor', bgColor, 'EdgeColor', [.75 .75 .75]);
                uistack(r,'bottom')
                
                xlim(getxlim)
                ylim(getylim)
                ax = gca;
                ax.YAxis.Exponent = 4;
                ax.XAxis.Exponent = 4;
            end    
        end
    end
end

% Add a colorbar to indicate correlation strength
colormap(cmap);
c = colorbar;

c.Label.String = 'Correlation'; 
c.Position= 'out';

 exportgraphics(gcf, ['/home/labs/barkailab/jingliu/Documents/Finalfiles/Fig1_newversion/', ...
     'Fig2_FTop300.pdf'], 'ContentType', 'vector');
