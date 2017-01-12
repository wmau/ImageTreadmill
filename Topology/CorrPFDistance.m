function CorrPFDistance(md,fieldType)
%CorrPFDistance(md,fieldType)
%
%   Does pairwise correlations of place fields then looks for anatomical
%   distance between place cells. This initial version looks only at place
%   cells and does a Spearman correlation for all pairs. Then it only takes
%   the significant place field correlations (.05, Bonferroni-corrected)
%   and does a Pearson correlation with anatomical distance. 
%

%%
    DATA = CompileMultiSessionData(md,{'placecells',fieldType});
    PlaceCells = DATA.placecells{1};
    nPCs = length(PlaceCells);
    
%% 
    centroids = getNeuronCentroids(md,'neurons',PlaceCells);
    centroids = centroids(~isnan(centroids(:,2)),:);
    
    D = nan(nPCs);
    for n1=1:nPCs
        %Centroid for cell 1.
        x1 = centroids(n1,1);
        y1 = centroids(n1,2); 
        
        for n2=n1+1:nPCs
            %centroid for cell 2.
            x2 = centroids(n2,1);
            y2 = centroids(n2,2);

            %Anatomical distance. 
            D(n1,n2) = sqrt((x2-x1)^2 + (y2-y1)^2);
        end
    end
    
%% 
    [R,p] = deal(nan(nPCs));
    for n1=1:nPCs
        PF1 = DATA.(fieldType){1}{PlaceCells(n1)};
        
        for n2=n1+1:nPCs
            PF2 = DATA.(fieldType){1}{PlaceCells(n2)};
            
            [R(n1,n2),p(n1,n2)] = corr(PF1(:),PF2(:),'rows','complete',...
                'type','spearman');
        end
    end
    
    Dflat = D(:);
    Rflat = R(:);
    p = p(:);
    nComparisons = sum(~isnan(p));
    good = p < .05/nComparisons & Rflat > 0;
    
    [~,pval] = corr(Dflat(good),Rflat(good),'rows','complete');
    s = scatter(Dflat(good),Rflat(good),10,'filled');
    alpha(s,.1);
    title(['p = ',num2str(pval)]);
    xlabel('Anatomical distance [microns]');
    ylabel('Correlation Coefficient');
    
end