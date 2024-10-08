function [omegaMat,xiMat,omegaxiMat] = CMT2omegadc_xi0_inter(M,bfigure)
%CMT2OMEGADC_XI0 compute the omegadc and xi0 angles between two moment tensors
%
% INPUT
%   M           6 x n moment tensors: M = [M11 M22 M33 M12 M13 M23]
%
% OUTPUT
%   omegaMat    n x n matrix of omegaDC(Mi,Mj) for angle between Mi and Mj
%   xiMat       n x n matrix of xi0(Mi,Mj) for angle between Mi and Mj
%   omegaxiMat  n x n matrix with omegaDC on lower ttriangl, xi0 on upper triangle
%
% See CMT2omegadc_xi0.m for details.
% 
% Carl Tape, 2024-09-19
%

if nargin==1, bfigure=true; end
[~,n] = size(M);

omegaMat = zeros(n,n);
xiMat = zeros(n,n);
for ii=1:n
    for jj=1:n  % this will compute 2x what we need
        if ii==jj, continue; end
        [omegaMat(ii,jj),xiMat(ii,jj)] = CMT2omegadc_xi0(M(:,ii),M(:,jj),0,0);
        %if imag(xiMat(ii,jj))
        %    ii,jj,M(:,ii),M(:,jj),omegaMat(ii,jj),xiMat(ii,jj)
        %    error
        %end
    end
end

omegaMat
xiMat

triu(omegaMat)

% square matrix with lower triangle as omegaDC and upper triangle as xi0
omegaxiMat = omegaMat;
omegaxiMat(tril(true(size(omegaxiMat)))) = xiMat(tril(true(size(xiMat))));

if bfigure
    bgridlines = true;
    btextlabels = true;
    bticklabels = true;
    dotsize = 4;
    
    if n > 10, btextlabels = false; bticklabels = false; bgridlines = false; dotsize = 2; end
    
    Z = omegaxiMat;
    
    figure;
    imagesc(Z);
    colormap('hot'); colorbar; axis equal tight;

    [nRows, nCols] = size(Z);
    if btextlabels
        for i = 1:nRows
            for j = 1:nCols
                text(j, i, sprintf('%.1f', Z(i,j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.5 0.5 0.5]);
            end
        end
    end
    
    if bticklabels
        Mlabels = arrayfun(@(x) sprintf('M_{%d}', x), 1:nCols, 'UniformOutput', false);
        set(gca, 'XTick',1:nCols,'YTick',1:nRows,'XTickLabel',Mlabels,'YTickLabel',Mlabels,...
            'TickLength',[0 0],'XAxisLocation','top');
    else
        set(gca, 'XTick',[],'YTick',[]);
    end
    
    if bgridlines
        % custom gridlines
        hold on
        for i = 0.5:nRows, plot([0.5, nCols+0.5], [i, i], 'k', 'LineWidth', 1); end
        for j = 0.5:nCols, plot([j, j], [0.5, nRows+0.5], 'k', 'LineWidth', 1); end
        hold off
    end
    
    title({sprintf('\\xi_0(M_i, M_j) lower \\Delta, \\omega_{DC}(M_i, M_j) upper \\Delta'),sprintf('%i moment tensors',n)});
    
    % scatterplot
    domega = omegaMat(triu(true(size(omegaMat)),1));
    dxi    = xiMat(triu(true(size(xiMat)),1));
    figure; plot(domega,dxi,'ko','markersize',dotsize,'markerfacecolor','k');
    xlabel('\omega_{DC}, degrees');
    ylabel('\xi_0, degrees');
    grid on;
end
    
%==========================================================================
% EXAMPLES

if 0==1
    %M = load('~/Downloads/MT_4_Carl.txt')';
    n = 100;  % try 10 or 100
    M = uniformMT(n,0,0);
    [omegaMat,xiMat,omegaxiMat] = CMT2omegadc_xi0_inter(M);
    figure(1); print(gcf,'-dpng',sprintf('~/omegadc_xi0_matrix_n%i',n));
    figure(2); print(gcf,'-dpng',sprintf('~/omegadc_xi0_scatter_n%i',n));
end

%==========================================================================