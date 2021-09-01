function figo = PlotHeatMap(x,y,v,fig,range)

    try
        if(isnumeric(fig))
            fig = figure(fig);
        end
    catch
        fig = figure();
    end
    figo = fig;

    F = scatteredInterpolant(x,y,v);
    [xq,yq] = meshgrid(-3:0.1:3);
    F.Method = 'linear';
    vq = F(xq,yq);
    figure(fig);
    contourf(xq,yq,vq,150,'LineColor','none');
    colorbar;
    
    try 
        isempty(range);
    catch
        range = [min(v),max(v)];
    end
    
    if ((min(v)<range(1)) || (max(v)>range(2)) || (range(1)>range(2)))
        error('range is not set properly.');
    end
    
    ColorStepSize = length(fig.Colormap);
    ColorLower0 = [1,1,1];
    ColorUpper0 = [1,0,0];
    
    ColorLower = Affine(ColorLower0,ColorUpper0,(min(v)-range(1))/(range(2)-range(1)));
    ColorUpper = Affine(ColorLower0,ColorUpper0,(max(v)-range(1))/(range(2)-range(1)));
    GradRed     = linspace(ColorLower(1),ColorUpper(1),ColorStepSize)';
    GradGreen   = linspace(ColorLower(2),ColorUpper(2),ColorStepSize)';
    GradBlue    = linspace(ColorLower(3),ColorUpper(3),ColorStepSize)';
    fig.Colormap = [GradRed GradGreen GradBlue];
    %colormap([GradRed GradGreen GradBlue]);

%     figure;
%     mesh(xq,yq,vq);
%     hold on;
%     plot3(x,y,v,'.');
    
end

function m = Affine(x,y,r)
    m = x*(1-r) + y*r;
end